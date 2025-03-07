/*

Cell-ID is intended to identify cells in images and to calculate a
number of statistics, including statistics derived from corresponding
fluorescence images.
Copyright (C) 2005 Andrew Gordon

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

Andrew Gordon can be contacted at
agordon@molsci.org
Molecular Sciences Institute
2168 Shattuck Ave, 2nd Floor
Berkeley, CA  94704

**********************************************************************
Start-copyright-notice-for-libtiff
Libtiff software is used by Cell-ID for some of the reading in of
TIF image file data and also for creating new TIF files. Libtiff is
available at http://www.remotesensing.org/libtiff/. The libtiff software
was written by Sam Leffler while working for Silicon Graphics and
contains the following copyright notice:

   "Copyright (c) 1988-1997 Sam Leffler
    Copyright (c) 1991-1997 Silicon Graphics, Inc.

    Permission to use, copy, modify, distribute, and sell this software and
    its documentation for any purpose is hereby granted without fee, provided
    that (i) the above copyright notices and this permission notice appear in
    all copies of the software and related documentation, and (ii) the names
    of Sam Leffler and Silicon Graphics may not be used in any advertising or
    publicity relating to the software without the specific, prior written
    permission of Sam Leffler and Silicon Graphics.

    THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND,
    EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
    WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

    IN NO EVENT SHALL SAM LEFFLER OR SILICON GRAPHICS BE LIABLE FOR
    ANY SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND,
    OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
    WHETHER OR NOT ADVISED OF THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF
    LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
    OF THIS SOFTWARE."

End-copyright-notice-for-Libtiff
*********************************************




*/

// rcel22 imports
#include <unistd.h>
#include <getopt.h>
#include <libgen.h>           // for basename()

// original imports
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nums.h"
#include "tif_routines.h"
#include "segment.h"
#include "date_and_time.h"
#include "image_type.h"
#include "align_image.h"
#include "split_and_overlap.h"
#include "parameters.h"
#include "flatten.h"

//V1.2a
#include "file_type.h"

//V1.2a TEST
#include "oif.h"

int seconds_per_day=86400;
#define max_files 1000
#define max_strlen 500

float *flat_cors=NULL;

//new_phase tells if we have a new file to search for cells.  It's
//external since if there's not a new file, we don't need to try to
//match the cells again to their location in the list of known cells,
//when we calculate new fluorescence values.  We declare it here since
//it's extern since segment.c uses it to decide whether to do the match
//against the known cells or just to use the previous values.
int new_phase=0;

// rcell2 getopt declarations
extern char *optarg;
extern int optind, opterr, optopt;
int opt;

int main(int argc, char *argv[]){

  //The first argument should point to a file that contains a list of
  //phase contrast files.
  //The second argument should point to a file that contains a list of
  //fluorescence files.

  //The first list should be a series of phase contrast pictures,
  //and the second list should be a series of fluorescent images of
  //the _same_ pictures.  Ie, a 1-1 correspondence, the first picture in
  //the first list corresponds to the first pictures in the second list.

  //We're going to search the phase contrast picture for the cells and
  //their boundaries and then calculate how much fluorescence is in each
  //one with the second picture.

  //We'll loop over all of them, tracking the cells through the frames.

  //V1.2a When reading oif (Olympus Image File) the first arguments is the
  //name of the file ("myfile.oif") and the second argument is the output

  // Welcome message
  printf("\n\n*** Cell_ID Version 1.5.0 ***");
  printf("** 2022 redistribution. Now with awesome identified masks, and tailored to rcell2's workflow. **\n\n");

  int out_mask=0;        // mask_mod: optionally output interior/boundary coords to TSV, default disabled
  int label_cells=0;     // mask_mod: label cells in BF.out optionally, default disabled
  int mask_output=0;     // mask_mod: cellID-adjusted mask output option in BF.out TIFF, default disabled
  int fill_interior=0;   // mask_mod: fill cells in BF.out. overrides label_cells. default disabled
  int interior_offset=0; // mask_mod: offset interior and boundary pixel intensities

  FILE *fp_in;
  FILE *fp;

  int bf_index[max_files]; //V1.2a Index of the BF file used to find cells,
                           //for each fluorescent file
  FILE *bf_fl_file=NULL;
  char bf_fl_file_name[500];

  char *flat_files[max_files]; //For flattening image
  char *dark_files[max_files]; //For subtracting fixed background
  char *phase_files[max_files];
  char *fluor_files[max_files];
  char *third_files[max_files];
  int flag[max_files];
  int flag_bf[max_files];
  char c0,c1,c2;  // for flag/channel identification from file names

  int time_flag[max_files];
  int time_index[max_files];

  char *s;
  char *s2;
  char *s3;
  char *s4;
  char line[500];
  char line2[500];
  char line3[500];
  int i,j,k;
  int n_fluor,n_phase,n_third,n_dark,n_flat;
  int iloop,nloop;

  float *ftmp;

  int flat_d[max_files],flat_t[max_files],flat_dmin,flat_cur=-1;
  int fl_d[max_files],fl_t[max_files],fl_dmin;
  int ph_d[max_files],ph_t[max_files],ph_dmin;
  int third_d[max_files],third_t[max_files],third_dmin;
  int dtmp,t,d_cur,t_cur;
  unsigned short year;
  unsigned char month,day;
  int hours,minutes,seconds;
  int dt,dt_min,j_min;
  int j_cur=0;
  int third_cur=-999;
  int output_third_image=0;
  int first=0;
  int xmax_new;
  int ymax_new;

  //Images
  float *bf=NULL;
  float *fl=NULL;
  float *third_image=NULL;
  int *bf_fl_labels=NULL;
  int *third_labels=NULL;
  int xmax=-999;
  int ymax=-999;

  int *fret_label_array=NULL;

  int align_fl=0;

  int align_individual_cells=0;
  int align_individual_cells_boundary=0;
  int force_nucleus_in_center=0;

  int output_individual_cells=0;

  int fret_image;

  float xstage,ystage;
  float tmp,tmp1,tmp2;
  int itmp;
  double double_tmp;

  char append[2];
  int do_append;

  int nuclear_fret_lower;

  float *dark=NULL;
  float exposure[max_files];
  float t_exposure;

  double max_d_over_s_cut_t0,max_d_over_s_cut_save;
  double max_split_d_over_minor_t0,max_split_d_over_minor_save;


  #define n_recomb_cuts_max 100
  float recombination_cuts[n_recomb_cuts_max];
  int recombination_cuts_type[n_recomb_cuts_max];
  int recombination_cuts_flag[n_recomb_cuts_max];
  int do_recombination=0;
  int n_recomb_cuts=0;

  int treat_brightfield_as_fluorescence=0;

  int i_last_find_cell_call=0;

  int n_found_cells;
  struct point **found_boundary_points=NULL;
  struct point **found_interior_points=NULL;

  third_image=NULL;                //Default in case don't have third image type
  third_image_type=no_third_image; //default type, value defined in image_type.h as "#define no_third_image 4"
  nuclear_fret_lower=-1;
  recalculate_internal=1;          //First time through always do
  image_type=bright_field;         //default, image_type.h says: "#define bright_field 0"
  fret_image=0;                    //default to false (image_type for fret images also tells what _kind_ of fret image, so set this fret_image flag also).

  n_third=-1;
  overall_id_offset=-1; //flag value, will be changed to 0 later
  strcpy(append,"w");
  do_append=0;
  for(i=0;i<max_files;i++){
    flag_bf[i]=0; //for case of treat_bright_field_as_fluorescnece,
                  //this will mark which files were the brighr field files
  }

  //Check if there is a file with certain parameters to load
  //First set parameters to default values
  max_d_over_s_cut_t0=-999.0;
  max_split_d_over_minor_t0=-999.0;
  max_split_d_over_minor=0.5;
  max_d_over_s_cut=6.0;
  max_pixels_per_cell=1500;
  min_pixels_per_cell=75;
  background_reject_factor=3.0;
  I_over_U_for_match=0.2;
  nucleus_radii[0]=2; //V1.4.5 radii size initialization
  nucleus_radii[1]=3;
  nucleus_radii[2]=4;
  nucleus_radii[3]=5;
  nucleus_radii[4]=6;
  nucleus_radii[5]=7;


  file_type=list_file; //V1.2a
  bf_fl_mapping=bf_fl_mapping_time;//V1.2a
  int paw_output=0; //V1.3a
  int n_bf_as_fl=0; //V1.4.2



  //Command line parsing variables

  char *equal_sign = NULL;
  int help_flag = 0;
  int time_of_t0_flag = 0;  // Default disabled

  char *param_file = "parameters.txt";
  char *dark_list_file = "dark.txt";
  char *flat_list_file = "flat.txt";

  char *bright_list_file = NULL;
  char *fluor_list_file = NULL;
  char *third_list_file = NULL;
  char *output_basename = NULL;

  char *fret_bf = NULL;
  char *fret_nuclear = NULL;
  char *str_align_fl = NULL;

  //V1.4a
  char *file_basename = NULL; // rcell2: used to be gchar, unsure if this is ok
                                                            // for getting correct channel identification from filenames with paths

  char str_third_img_label[500];
  char *pnt_third_img_label = NULL;
  char str_image_type[500];
  char *pnt_image_type = NULL;

  // rcell2 chunk: getopts option parser
  opterr = 0;  // https://stackoverflow.com/a/24331449/11524079
  optind = 1;  // https://stackoverflow.com/a/25937743/11524079

  while((opt = getopt(argc, argv, "p:b:f:g:o:D:F:limtwzh")) != -1) {
    printf("Parsing getopt options\n");
    switch(opt) {
    case 'p':
       printf(" -p parameters: ");
       printf("%s\n", optarg);
      param_file=optarg;
      break;

    case 'b':
       printf(" -b brightfield: ");
       printf("%s\n", optarg);
      bright_list_file=optarg;
      break;

    case 'f':
       printf(" -f fluorescence: ");
       printf("%s\n", optarg);
      fluor_list_file=optarg;
      break;
      
    case 'g':
      printf(" -g third image (nucleus or vacuole): ");
      printf("%s\n", optarg);
      third_list_file=optarg;
      break;
      
    case 'o':
       printf(" -o output_prefix: ");
       printf("%s\n", optarg);
      output_basename=optarg;
      break;

    case 'D':
      printf(" -D dark file list: ");
      printf("%s\n", optarg);
      dark_list_file=optarg;
      break;
      
    case 'F':
      printf(" -F flat file list: ");
      printf("%s\n", optarg);
      flat_list_file=optarg;
      break;
      
    case 'l':
       printf(" -l Label cells in BF.\n");
       if(fill_interior==0) label_cells=1;
      break;

    case 'i':
       printf(" -i Fill interior pixels in BF.out\n");
       fill_interior=1;
       label_cells=0;
      break;

    case 'm':
       printf(" -m Replace BF.out with segmentation masks only, removing image data.\n");
       mask_output=1;
      break;

    case 't':
       printf(" -t Output cell boundary and interior coords to TSV file.\n");
       out_mask=1; // enable
      break;

    case 'w':
       printf(" -w Offset boundary and interior mask intensities.\n");
       interior_offset = 1; // enable
      break;

    case 'z':
      printf(" -z Time of t0 flag enabled.\n");
      time_of_t0_flag = 1;
      break;

    case 'h':
        printf(" -h Help flag enabled.\n");
        help_flag = 1;
      break;

    case ':':
       printf(" - option needs a value\n");
      break;

    case '?':
       printf(" - unknown option: ");
       printf("%c\n", optopt);
      break;
    }
  }

  // mask_mod: Check and set mask_output type
  if((fill_interior==1||label_cells==1)&&mask_output==0) mask_output=2;

  // Get all of the non-option arguments and print them
  //https://azrael.digipen.edu/~mmead/www/Courses/CS180/getopt.html
  if (optind < argc) {
    printf(" - Non-option args: ");
    while (optind < argc){
      printf("%s ", argv[optind++]);
    }
    printf("\n");
  }

  //Checking for parameters file option in command line manually
  for(i=1;i<argc;i++){
    if(strstr(argv[i],"-?")!=NULL||strstr(argv[i],"--help")!=NULL||argc==1){
      help_flag = 1;
      break;
    }
    if(strcmp(argv[i],"-p")==0 || strcmp(argv[i],"--param")==0){
      equal_sign=strstr(argv[i],"=");
      if(equal_sign==NULL){
        if(i+1<argc && argv[i+1][0]!='-'){
          param_file=argv[i+1];
        }else{
          printf("Filename required after -p or --param option");
          return 3;  // Exit code for bad parameters argument.
        }
      } else {
        param_file=++equal_sign;
      }
    }
  }

  if(help_flag==0 && (fp=fopen(param_file,"r"))!=NULL ){
    printf("Reading %s\n",param_file);
    //fscanf() only reads to next white space, not end of line
    while((fgets(line,450,fp))!=NULL){ //next line (while not EOF)
      if (line[0]!='#'){ //if not a comment
          if ((strstr(line,"max_split_over_minor"))!=NULL){
          if ((strstr(line,"_t0"))!=NULL){
              sscanf(line,"%s %le",line2,&double_tmp);
              max_split_d_over_minor_t0=double_tmp;
             }else{
              sscanf(line,"%s %le",line2,&double_tmp);
              max_split_d_over_minor=double_tmp;
            }
          }else if((strstr(line,"max_dist_over_waist"))!=NULL){
            if ((strstr(line,"_t0"))!=NULL){
              sscanf(line,"%s %le",line2,&double_tmp);
              max_d_over_s_cut_t0=double_tmp;
            }else{
               sscanf(line,"%s %le",line2,&double_tmp);
               max_d_over_s_cut=double_tmp;
            }
          }else if((strstr(line,"background_reject_factor"))!=NULL){
             sscanf(line,"%s %le",line2,&double_tmp);
             background_reject_factor=double_tmp;
          }else if((strstr(line,"tracking_comparison"))!=NULL){
             sscanf(line,"%s %le",line2,&double_tmp);
             I_over_U_for_match=double_tmp;
          }else if((strstr(line,"paw_output"))!=NULL){
             paw_output=1;
          }else if((strstr(line,"max_pixels_per_cell"))!=NULL){
             sscanf(line,"%s %i",line2,&itmp);
             max_pixels_per_cell=itmp;
          }else if((strstr(line,"min_pixels_per_cell"))!=NULL){
             sscanf(line,"%s %i",line2,&itmp);
             min_pixels_per_cell=itmp;
          }else if((strstr(line,"force_nucleus_in_center"))!=NULL){
             force_nucleus_in_center=1;
          }else if((strstr(line,"output_individual_cells"))!=NULL){
             output_individual_cells=1;
             //system("mkdir -p cells");

          //image types:
          }else if(strstr(line,"image_type")!=NULL){
             sscanf(line,"%s %s",line2,str_image_type);


        //V1.2a bright field to fluorescence mapping
        }else if(strstr(line,"bf_fl_mapping")!=NULL){
          sscanf(line,"%s %s",line2,line3);
          if(strstr(line3,"time")!=NULL){
            //printf("Mapping bright field and fluorescence images by time.\n");
            bf_fl_mapping=bf_fl_mapping_time;
          } else if(strstr(line3,"list")!=NULL){
            //printf("Mapping bright field and fluorescence images by list order.\n");
            //printf("Same number of elemtes required in list files.\n");
             bf_fl_mapping=bf_fl_mapping_list;
          } else {
            printf("-%s- is a invalid value for bf_fl_mapping in parameter.txt.\n",line3);
            printf("Using time mapping by default.\n");
            bf_fl_mapping=bf_fl_mapping_time;
          }

          }else if(strstr(line,"fret")!=NULL){
            sscanf(line,"%s %s",line2,line3);

            if (fret_image!=1){ //if we haven't been here before
              fret_image=1; //To indicate that we have a fret split image
              //printf("Searching a split FRET image.\n");
            }
            //Further check the fret line:
            if (strstr(line3,"bf_top_only")!=NULL){
              //Now check if the split-brightfield has an image on top only or
              //bottom only, or bottom and top (note that each of the argument
              //strings has "fret" in them, so they'll pass the above if
              //statement).
              image_type=fret_bf_top_only;
              //printf("Cells on top part of BF image only.\n");
            }else if((strstr(line3,"bf_bottom_only"))!=NULL){
              image_type=fret_bf_bottom_only;
              //printf("Cells on bottom part of BF image only.\n");
            }else if((strstr(line3,"bf_bottom_and_top"))!=NULL){
              image_type=fret_bf_bottom_and_top;
              //printf("Cells on bottom and top part of BF image.\n");
            }

            //Check if should use top or bottom of image that's used
            //for nuclear label. (If this isn't set, then it defaults
            //to bottom of image if there is a third image set and to
            //the top of the image otherwise--this is done below all this).
            if(strstr(line3,"nuclear_top")!=NULL){
              nuclear_fret_lower=1;
            }else if(strstr(line3,"nuclear_bottom")!=NULL){
              nuclear_fret_lower=0;
            }

          }else if(strstr(line,"align_fl_to_first")!=NULL){
            //Check whether we want to use the first fluorescence image to
            //align all the others
            //printf("Will align all FL files to first FL file.\n");
            align_fl=1;
          }else if(strstr(line,"align_fl_to_bf")!=NULL){
            //or whether to align fl images to bf
            //printf("Will align first FL files to brightfield.\n");
            align_fl=2;
            //Check for third list of images. The use of the
            //images depends on the image_type.
          }else if(strstr(line,"treat_brightfield_as_fluorescence_also")!=NULL){
            treat_brightfield_as_fluorescence=1;
            //printf("Adding BF image as additional fluorescence image.\n");
          }else if (strstr(line,"third_image")!=NULL){
            sscanf(line,"%s %s",line2,str_third_img_label);
            // Label image: nuclear or vacuole.
            // Here the "line" is read, the first value is stored in "line2",
            // and the second in "str_third_img_label", which is used later on.
            // As per "parameters_description.txt":
            //   ## If there is a third image type included, then this describes how
            //   ## to use it (nuclear_label means its a fluorescence tagged nucleus
            //   ## channel, vacuole_label means vacuole label)
            //   # third_image vacuole_label
            //   # third_image nuclear_label
          }else if(strstr(line,"align_individual_cells")!=NULL){
            align_individual_cells=1;
            //see if "boundary" is part of name
            if(strstr(line,"align_individual_cells_boundary")!=NULL){
              align_individual_cells_boundary=1;
              //printf("Will wiggle each cell around to re-align with BF");
              //printf(" using boundary.\n");
            }//else{
              //printf("Will wiggle each cell around to re-align with BF.\n");
            //}
          }else if(strstr(line,"do_recombination")!=NULL){
            do_recombination=1;
          }else if(strstr(line,"recombination_fl_per_a_nucleus")!=NULL){
            sscanf(line,"%s %i %e",line2,&itmp,&tmp);
            recombination_cuts_type[n_recomb_cuts]=0;
            recombination_cuts[n_recomb_cuts]=tmp;
            recombination_cuts_flag[n_recomb_cuts]=itmp;
            printf("Will try to recombine cells with nucleus fluorescence\n");
            printf("intensity below %e for fl images with flag=%i.\n",
                recombination_cuts[n_recomb_cuts],
                recombination_cuts_flag[n_recomb_cuts]);

          n_recomb_cuts++;
            if (n_recomb_cuts>=n_recomb_cuts_max){
              printf("Too many recombination cuts, dropping %s.\n",line);
              n_recomb_cuts--;
            }
          }else if(strstr(line,"recombination_is_a_cell_fl_per_a")!=NULL){
            sscanf(line,"%s %i %e",line2,&itmp,&tmp);
            recombination_cuts_type[n_recomb_cuts]=1;
            recombination_cuts[n_recomb_cuts]=tmp;
            recombination_cuts_flag[n_recomb_cuts]=itmp;
            printf("For recomb only, considering cells with total FL\n");
            printf("intensity above ");
            printf("%e for fluorescence images with flag=%i.\n",
                 recombination_cuts[n_recomb_cuts],
                 recombination_cuts_flag[n_recomb_cuts]);

          n_recomb_cuts++;
            if (n_recomb_cuts>=n_recomb_cuts_max){
              printf("Too many recombination cuts, dropping %s.\n",line);
              n_recomb_cuts--;
            }
          }else if(strstr(line,"recombination_all_new_cells")!=NULL){
            recombination_cuts_type[n_recomb_cuts]=2;
            recombination_cuts[n_recomb_cuts]=0.0;
            recombination_cuts_flag[n_recomb_cuts]=0;
            printf("All cells produced after i_t 0 will be recombined.\n");
            n_recomb_cuts++;
            if (n_recomb_cuts>=n_recomb_cuts_max){
              printf("Too many recombination cuts, dropping %s.\n",line);
              n_recomb_cuts--;
            }
          }else if(strstr(line,"append_output")!=NULL){
            sscanf(line,"%s %i",line2,&itmp);
            overall_id_offset=itmp;
          }else if(strstr(line,"nucleus_radius_1")!=NULL){
             sscanf(line,"%s %i",line2,&itmp);
                     nucleus_radii[0]=itmp;
          }else if(strstr(line,"nucleus_radius_2")!=NULL){
             sscanf(line,"%s %i",line2,&itmp);
                     nucleus_radii[1]=itmp;
          }else if(strstr(line,"nucleus_radius_3")!=NULL){
             sscanf(line,"%s %i",line2,&itmp);
                     nucleus_radii[2]=itmp;
          }else if(strstr(line,"nucleus_radius_4")!=NULL){
             sscanf(line,"%s %i",line2,&itmp);
                     nucleus_radii[3]=itmp;
          }else if(strstr(line,"nucleus_radius_5")!=NULL){
             sscanf(line,"%s %i",line2,&itmp);
                     nucleus_radii[4]=itmp;
          }else if(strstr(line,"nucleus_radius_6")!=NULL){
             sscanf(line,"%s %i",line2,&itmp);
                     nucleus_radii[5]=itmp;
          }
      } //End of check that first character wasn't a "#"
    } //End of while loop over parameters.txt
    fclose(fp);
  }else if (help_flag==0) {   //End of check that parameters.txt was open ok
    printf("\n%s not found. Using default parameters. \n",param_file);
  }

  if(pnt_third_img_label==NULL) pnt_third_img_label=&str_third_img_label[0];
  if(pnt_image_type==NULL) pnt_image_type=&str_image_type[0];

  // Cell-ID help message
  // Params: p:b:f:o:D:F:limtwzh
  if(help_flag==1){
    printf("\n\n");
    printf("CellID help:\n");
    printf("  -- Usage: cell -b BF/list/file.txt -f FL/list/file.txt -p parameters/list/file.txt -o output/path/prefix/Position4/out [-D <dark.tiff>] [-F <flat.tiff>] [-l] [-i] [-m] [-t] [-w] [-h] \n");
    printf("  -- For detailed usage, visit: https://github.com/darksideoftheshmoo/cellid-linux/tree/mask_mod\n");
    printf("  -- Parameters: \n");
    printf("     -p <path> Path to a text file containing a list of segmentation parameters.\n");
    printf("     -b <path> Path to a text file containing a list of paths to 'BF' images used for segmentation.\n");
    printf("     -f <path> Path to a text file containing a list of paths to 'FL' images used for quantification.\n");
    printf("     -o <path> Output path prefix; the relative or absolute path where output files will be saved, including the file prefix (tipically 'out').\n");
    printf("     -g <path> Path to a text file containing a list of paths to 'third' images used for estimating the position of the nucleus or vacuole.\n");
    printf("     -D <path> Path to a 'dark' image used for background correction (optional, correction skipped by default).\n");
    printf("     -F <path> Path to a 'flat' image used for flattening correction (optional, correction skipped by default).\n");
    printf("     -l        Set cell boundary pixel intensities proportional to each cellID, following the relationship 'cellID = 65535 - boundary_intensity - 1'. Also adds cellID numbers to the cells, with maximum pixel intensity (65535). (optional, enabled by default).\n");
    printf("     -i        Set cell boundary and interior pixel intensities proportional to each cellID (optional, disabled by default). This overrides cell labeling with '-l'.\n");
    printf("     -m        Replace contents of the BF.out TIFF files with segmentation masks only, removing image data (optional, disabled by default).\n");
    printf("     -t        Output cell boundary and interior coords to a compressed TSV file (optional, disabled by default).\n");
    printf("     -w        Offset boundary and interior pixel intensities by a calculated 'interior_offset' threshold. cellID will relate to interior pixel intensities with the relationship 'cellID = 65535 - boundary_intensity - interior_offset - 1'. The offset defaults to 5000, but may have a larger value for images or time series with more than 2500 cells.\n");
    printf("     -z        Write the absolute time of the first image to a text file (optional, disabled by default).\n");
    printf("     -h        Show the help message.\n");
    printf("\n\n");
    return 0;
  }

  if(output_individual_cells==1){
      system("mkdir -p cells");
  }

  //V1.4.5
  printf("Using nucleus radii %i %i %i %i %i %i px.\n",nucleus_radii[0],nucleus_radii[1]
              ,nucleus_radii[2],nucleus_radii[3],nucleus_radii[4],nucleus_radii[5]);

  //Reading fret parameters form command line and steping on parameters.txt
  //if required.
  if(fret_bf != NULL){
    fret_image=1;
    if (strstr(fret_bf,"top")!=NULL){
        //Now check if the split-brightfield has an image on top only or
        //bottom only, or bottom and top
      image_type=fret_bf_top_only;
      //printf("Cells on top part of BF image only.\n");
    }else if((strstr(fret_bf,"bottom"))!=NULL){
        image_type=fret_bf_bottom_only;
        //printf("Cells on bottom part of BF image only.\n");
      }else if((strstr(fret_bf,"both"))!=NULL){
        image_type=fret_bf_bottom_and_top;
        //printf("Cells on bottom and top part of BF image.\n");
      }
  }
    if(fret_nuclear!=NULL){
      fret_image=1;
      if(strstr(fret_nuclear,"top")!=NULL){
        nuclear_fret_lower=1;
      }else if(strstr(fret_nuclear,"bottom")!=NULL){
        nuclear_fret_lower=0;
      }
  }

  //Printing the final fret and image type parameters
  if (fret_image==1){
     printf("Searching a split FRET image.\n");
    if(image_type==fret_bf_top_only){
      printf("Cells on top part of BF image only.\n");
    }else if(image_type==fret_bf_bottom_only){
        printf("Cells on bottom part of BF image only.\n");
      }else if(image_type==fret_bf_bottom_and_top){
        printf("Cells on bottom and top part of BF image.\n");
      }
    //If nuclear_fret_lower wasn't set, set to defaults
    if (nuclear_fret_lower==(-1)){
      if (third_image_type==no_third_image){
          nuclear_fret_lower=1; //Upper part if third image
      }else{
          nuclear_fret_lower=0; //Lower part if no third image
      }
    }
    if (nuclear_fret_lower==1){
      printf("Using top part of nuclear-image for nucleus.\n");
    }else{
      printf("Using bottom part of nuclear-image for nucleus.\n");
    }

  } else {
    //image_type is either a fret type, or one of the following
    if (strstr(pnt_image_type,"bright")!=NULL){
        image_type=bright_field;                                         // "#define bright_field 0" at image_type.h
        printf("Searching brightfield image for cells.\n");
      }else if(strstr(pnt_image_type,"decon")!=NULL){
        image_type=metamorph_deconvolution;
        printf("Searching metamorph_deconvolution image for cells.\n");
      }else if(strstr(pnt_image_type,"hex")!=NULL){
        image_type=hexagonal_grid;
        printf("Searching hexagonal_grid image type for cells.\n");
      //}else if(strstr(str_image_type,"membrane_tagged_fl")!=NULL){
      //  image_type=membrane_tagged_fl;
      //  printf("Searching a membrane tagged fluor image for cells.\n");
      }else if(strstr(pnt_image_type,"confocal")!=NULL){
        image_type=confocal_transmission;
        printf("Searching a confocal transmission image for cells.\n");
      }
  }

  //Stepping parameters.txt values
  if(str_align_fl!=NULL){
    if(strstr(str_align_fl,"first")!=NULL){
      align_fl=1;
    }else if(strstr(str_align_fl,"bf")!=NULL){
      align_fl=2;
    }else{
      printf("'%s' is a invalid value for --align-fl, no aligment will be done.\n"
        ,str_align_fl);
    }
  }

  //Printing final values
  printf("\nConfiguration for this run:\n");
  if(align_fl==1){
    printf("    Will align all FL files to first FL file.\n");
  }else if (align_fl==2){
    printf("    Will align first FL files to brightfield.\n");
  }


  if(treat_brightfield_as_fluorescence==1){
      printf("    Adding BF image as additional fluorescence image.\n");
  }

  //Cheking for third image type, note that the same string has been read in
  //parametes.txt and steped on with GOptions
  if((strstr(pnt_third_img_label,"nuclear")!=NULL)){
      third_image_type=nuclear_label;
      printf("    Third image will be used to label nucleus.\n");
    }else if(strstr(pnt_third_img_label,"vacuole")!=NULL){
      third_image_type=vacuole_label;
      printf("    Third image will be used to correct for vacuole.\n");
    }else{
      if(third_list_file!=NULL){
      printf("    Unknown third image type, ignoring third images.");
    }
      third_image_type=no_third_image;
    }

  //printing values
  if(align_individual_cells_boundary==1){
    align_individual_cells=1;
      printf("    Will wiggle each cell around to re-align with BF using boundary.\n");
    }else if(align_individual_cells==1){
    printf("    Will wiggle each cell around to re-align with BF.\n");
  }


  if (overall_id_offset>=0){

      strcpy(append,"a");
      do_append=1;
      printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      printf("!!!!              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      printf("!!!!    WARNING   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      printf("!!!!              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      printf("!!!!                                       !!!!!!\n");
      printf("!!!! We're going to OVERWRITE output data. !!!!!!\n");
      printf("!!!!                                       !!!!!!\n");
      printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      printf("All cells will be given ID offset %i.\n",
        overall_id_offset);

  }else if (overall_id_offset==-1){
    overall_id_offset=0;
  }else{
      printf("    Append_output was requested with ID offset of %i.\n",overall_id_offset);
      printf("    Minimum offset is 0.\n");
      exit(0);
    }


  //Read in the files
  n_third=0;
  if (third_list_file!=NULL){//if there is a third image
    if( (fp_in=fopen(third_list_file,"r"))==NULL ){
        printf("Couldnt open file %s.\n",third_list_file);
        return 7;  // Exit code for problems opening a file.
    }
    //Read in names of the third-files to be read in
    i=0;
    while(fscanf(fp_in,"%s",line)==1){
        third_files[i]=(char *)malloc(max_strlen*sizeof(char));
        strcpy(third_files[i],line);
        i++;
    }
    fclose(fp_in);
    n_third=i;
  }

  if ((n_recomb_cuts>0)&&(do_recombination==0)){
    printf("***** Recombination cuts were set but NOT doing the");
    printf("***** recombination check.");
  }

  printf("Using parameters:\n");
  printf("    max_dist_over_minor_axis=%le\n",max_split_d_over_minor);
  printf("    max_dist_over_waist=%le\n",max_d_over_s_cut);
  printf("    background_reject_factor=%le\n",background_reject_factor);
  printf("    max_pixels_per_cell=%i\n",max_pixels_per_cell);
  printf("    min_pixels_per_cell=%i\n",min_pixels_per_cell);
  printf("    I_over_U_for_match (tracking_comparison parameter)=%le\n",I_over_U_for_match);
  max_d_over_s_cut_save=max_d_over_s_cut;
  max_split_d_over_minor_save=max_split_d_over_minor;
  if (max_d_over_s_cut_t0>-10.0){
    printf("   For first image: max_dist_over_waist=%e\n",
       max_d_over_s_cut_t0);
  }
  if (max_split_d_over_minor_t0>-10.0){
    printf("   For first image: max_dist_over_minor_axis=%e\n",
       max_split_d_over_minor_t0);
  }

  printf("\nEnd checking of arguments.\n");
  //End checking of arguments

  //Read in the names of the phase and fluorescence files.
  printf("\nRead in the names of the phase and fluorescence files.\n");

  //loading brightfield file names into array phase_files
  if( (fp_in=fopen(bright_list_file,"r"))==NULL ){
    printf("Couldnt open file %s.\n",bright_list_file);
    return 7;  // Exit code for problems opening a file.
  }
  i=0;
  while(fscanf(fp_in,"%s",line)==1){ //the 1 means it filled the "%s" part
    phase_files[i]=(char *)malloc(max_strlen*sizeof(char));
    strcpy(phase_files[i],line);
    i++;
  }
  fclose(fp_in);
  n_phase=i;

  //loading fluorescence file names into array fluor_files
  if( (fp_in=fopen(fluor_list_file,"r"))==NULL ){
    printf("Couldnt open file %s.\n",fluor_list_file);
    return 7;  // Exit code for problems opening a file.
  }
  i=0;
  while(fscanf(fp_in,"%s",line)==1){ //the 1 means it filled the "%s" part
    fluor_files[i]=(char *)malloc(max_strlen*sizeof(char));
    strcpy(fluor_files[i],line);
    i++;
  }
  fclose(fp_in);
  n_fluor=i;

  if (treat_brightfield_as_fluorescence==1){
    int phase_exists; //mask_mod

    //if mapping list we have to duplicate the phase file list
    if(bf_fl_mapping==bf_fl_mapping_list){
      /*
      mask_mod: commented this section, as it incorrectly copied phase files.
      See comment below

      //adding first bf as fl
      phase_files[n_phase]=(char *)malloc(max_strlen*sizeof(char));
      strcpy(phase_files[n_phase],phase_files[0]);
      fluor_files[n_fluor]=(char *)malloc(max_strlen*sizeof(char));
      strcpy(fluor_files[n_fluor],phase_files[0]);

      n_bf_as_fl=1;

      //adding new bf as fl
      for(i=1;i<n_phase;i++){
          for(j=0;j<i;j++){
              if(strstr(phase_files[i],phase_files[j])==NULL){//new bf
                phase_files[n_phase+n_bf_as_fl]=(char *)malloc(max_strlen*sizeof(char));
                strcpy(phase_files[n_phase+n_bf_as_fl],phase_files[i]);
                fluor_files[n_fluor+n_bf_as_fl]=(char *)malloc(max_strlen*sizeof(char));
                strcpy(fluor_files[n_fluor+n_bf_as_fl],phase_files[i]);
                n_bf_as_fl++;
              }
          }
      }
      n_phase+=n_bf_as_fl;
      n_fluor+=n_bf_as_fl;
      */

      // mask_mod: this corrects the addition of 'phase' images to the 'fluor'
      // list. It first gets a list of unique phase (BF) images, and then
      // inserts them as the last FL channel in fluor_files, and at the
      // corresponding position in phase_files.

      // Get unique phase image names from phase_files
      char *phase_unique[n_phase];
      n_bf_as_fl=0;
      // Loop over all phase images in phase_files
      for(i=0;i<n_phase;i++){
        phase_exists=0;
        for(j=0;j<n_bf_as_fl;j++){
          // Check if phase image already exists in phase_unique
          if(strstr(phase_files[i],phase_unique[j])!=NULL) phase_exists=1;
        }
        // If it does not, add it to phase_unique
        if(phase_exists==0){
          phase_unique[n_bf_as_fl]=(char *)malloc(max_strlen*sizeof(char));
          strcpy(phase_unique[n_bf_as_fl],phase_files[i]);
          n_bf_as_fl++;
        }
      }

      // Add unique phase images into phase_files and fluor_files as last fluor
      // channel
      int insert_pos;
      // Loop over unique phase images in phase_unique
      for(i=0;i<n_bf_as_fl;i++){
        // Find last occurence of phase image in phase_files
        insert_pos=-1;
        for(j=0;j<n_phase;j++) if(strstr(phase_unique[i],phase_files[j])!=NULL) insert_pos=j+1;

        // This should always be true, but just in case...
        if(insert_pos>-1){
          // Allocate space at end of phase_files and fluor_files
          phase_files[n_phase]=(char *)malloc(max_strlen*sizeof(char));
          fluor_files[n_phase]=(char *)malloc(max_strlen*sizeof(char));
          // Move all files behind insert_pos one index forward
          for(j=n_phase;j>insert_pos;j--){
            phase_files[j]=phase_files[j-1];
            fluor_files[j]=fluor_files[j-1];
          }
          // Insert phase image at insert_pos in phase_files and fluor_files
          phase_files[insert_pos]=phase_unique[i];
          fluor_files[insert_pos]=phase_unique[i];
          flag_bf[insert_pos]=1; //mark that is a brightfield
          n_phase++;
          n_fluor++;
        }
      }

    }else{
      /*
      mask_mod: replaced by next for loop. See comment below.
      //Copy brightfield image names into fluor list
      for(i=0;i<n_phase;i++){
        fluor_files[n_fluor]=(char *)malloc(max_strlen*sizeof(char));
        strcpy(fluor_files[n_fluor],phase_files[i]);
        flag_bf[n_fluor]=1; //mark that is a brightfield
        n_fluor++;
      }
      */

      //mask_mod: this is an attempt to fix BF_as_FL when bf_fl_mapping==bf_fl_mapping_time.
      // However, the image sorting by acquisition time doesn't seem to work
      // correctly, so although this fixes the addition of 'phase' images
      // to the 'fluor' list, there are more errors downstream.
      //Copy brightfield image names into fluor list
      for(i=0;i<n_phase;i++){
        phase_exists=0;
        for(j=0;j<n_fluor;j++){
          if(strstr(phase_files[i],fluor_files[j])!=NULL) phase_exists=1;
        }
        if(phase_exists==0){
          fluor_files[n_fluor]=(char *)malloc(max_strlen*sizeof(char));
          strcpy(fluor_files[n_fluor],phase_files[i]);
          flag_bf[n_fluor]=1; //mark that is a brightfield
          n_fluor++;
        }
      }
    }
  }

  //dark files
  if( (fp_in=fopen(dark_list_file,"r"))==NULL ){
    printf("%s not found. No dark-field corrections.\n", dark_list_file);
    n_dark=0;
  }else{
    i=0;
    while(fscanf(fp_in,"%s",line)==1){ //the 1 means it filled the "%s" part
      dark_files[i]=(char *)malloc(max_strlen*sizeof(char));
      strcpy(dark_files[i],line);
      //Get exposure time of each. (Note that we're going to end up
      //reading these in every time, instead of saving them in memory. We
      //only do the correction once and we don't want to have them sitting
      //around in memory probably.)
      exposure[i]=get_exposure(dark_files[i]);
      printf("Found %s with exposure time: %f\n", dark_files[i], exposure[i]);
      i++;
    }
    fclose(fp_in);
    n_dark=i;
  }

  //Files to flatten image
  if( (fp_in=fopen(flat_list_file,"r"))==NULL ){
    printf("%s not found. No flattening-image corrections.\n",flat_list_file);
    n_flat=0;
  }else{
    i=0;
    while(fscanf(fp_in,"%s",line)==1){ //the 1 means it filled the "%s" part
      flat_files[i]=(char *)malloc(max_strlen*sizeof(char));
      strcpy(flat_files[i],line);
      i++;
    }
    fclose(fp_in);
    n_flat=i;
  }

  //We're going to keep track of the ordinal number of
  //each flag. E.g., the first time flag=1 appears we'll
  //increment time_flag[1], etc.
  for(i=0;i<max_files;i++){
    time_flag[i]=0;
    time_index[i]=0;
  }

  //V1.2b I
  if (bf_fl_mapping==bf_fl_mapping_time){ //Mapping by time
    //Loop over all the files and calculate dates and times of each
    //(have to pull this out of metamorph's special tag in the tiff file)
    //V1.2a TODO modify, beware of ph_t, ph_dmin, ph_d arrays for future use
    ph_dmin=(int)(0x1ffffff);
    for(i=0;i<n_phase;i++){
      if(get_date_and_time(phase_files[i],&dtmp,&t,&xstage,&ystage)==0){
        printf("Couldn't get date and time for %s.\n",phase_files[i]);
        dtmp=0.0;
        t=0.0;
        xstage=-99999.0;
        ystage=-99999.0;
      }
      t=t/1000; //Convert to seconds
      ph_t[i]=t;
      if (dtmp<ph_dmin){
        ph_dmin=dtmp;
      }
      ph_d[i]=dtmp;
      printf("Stage position: %e  %e\n",xstage,ystage);fflush(stdout);
    }

   //V1.2a TODO modify, beware of the fl_d, fl_dmin, fl_t arrays
    fl_dmin=(int)(0x1ffffff);
    for(i=0;i<n_fluor;i++){
      if(get_date_and_time(fluor_files[i],&dtmp,&t,&xstage,&ystage)==0){
        printf("Couldn't get date and time for %s.\n",fluor_files[i]);
        dtmp=0.0;
        t=0.0;
        xstage=-99999.0;
        ystage=-99999.0;
      }
      t=t/1000; //Convert to seconds
      if (dtmp<fl_dmin){
        fl_dmin=dtmp;
      }
      fl_d[i]=dtmp;
      fl_t[i]=t;
    }

    //V1.2a TODO modify, beware of the third_t, dmin, d arrays
    third_dmin=(int)(0x1ffffff);
    for(i=0;i<n_third;i++){
      if(get_date_and_time(third_files[i],&dtmp,&t,&xstage,&ystage)==0){
        printf("Couldn't get date and time for %s.\n",third_files[i]);
        dtmp=0.0;
        t=0.0;
        xstage=-99999.0;
        ystage=-99999.0;
      }
      t=t/1000;//Convert to seconds
      third_t[i]=t;
      if (dtmp<third_dmin){
        third_dmin=dtmp;
      }
      third_d[i]=dtmp;
      printf("Thirdtime %i: %i  %i\n",i,third_t[i],third_d[i]);
    }

    flat_dmin=(int)(0x1ffffff);
    for(i=0;i<n_flat;i++){
      if(get_date_and_time(flat_files[i],&dtmp,&t,&xstage,&ystage)==0){
        printf("Couldn't get date and time for %s.\n",flat_files[i]);
        return 5;  // Exit code for problems getting time info from image.
      }else{
        t=t/1000; //Convert to seconds
        flat_t[i]=t;
        if (dtmp<flat_dmin){
            flat_dmin=dtmp;
        }
        flat_d[i]=dtmp;
      }
    }

    //We're going to read in the data in the order of the
    //fluorescence files, and we match other images based on their
    //closeness in time to the fluorescence images. Here, we sort
    //the FL images by time (so even if the user sends them in a different
    //order, they'll always be run time-sorted)--(ie, if you don't want
    //this, comment out the sort here).
    for (i=0;i<(n_fluor-1);i++){
      for(j=i;j<n_fluor;j++){
        if ((((fl_d[i]-fl_dmin)*seconds_per_day)+(fl_t[i]))>
             (((fl_d[j]-fl_dmin)*seconds_per_day)+(fl_t[j]))){
          //Swap them
            d_cur=fl_d[i];
            t_cur=fl_t[i];
            strcpy(line,fluor_files[i]);
            fl_d[i]=fl_d[j];
            fl_t[i]=fl_t[j];
            strcpy(fluor_files[i],fluor_files[j]);
            fl_d[j]=d_cur;
            fl_t[j]=t_cur;
            strcpy(fluor_files[j],line);
            k=flag_bf[i];
            flag_bf[i]=flag_bf[j];
            flag_bf[j]=k;
        }else if ((((fl_d[i]-fl_dmin)*seconds_per_day)+(fl_t[i]))==
                    (((fl_d[j]-fl_dmin)*seconds_per_day)+(fl_t[j]))){
            //This is the case that the times are _equal_. This
            //could happen if the user has, say, symbolic links to
            //the same file to create holder-spaces at time points where
            //certain colors weren't obtained--to simplify the offline
            //data structure, etc. For this case, consider the first
            //letter of the name to sort by
            if ((*fluor_files[i])>(*fluor_files[j])){ //swap
              d_cur=fl_d[i];
              t_cur=fl_t[i];
              strcpy(line,fluor_files[i]);
              fl_d[i]=fl_d[j];
              fl_t[i]=fl_t[j];
              strcpy(fluor_files[i],fluor_files[j]);
              fl_d[j]=d_cur;
              fl_t[j]=t_cur;
              strcpy(fluor_files[j],line);
              k=flag_bf[i];
               flag_bf[i]=flag_bf[j];
              flag_bf[j]=k;
            }//end swap if
        }//end equal time if
      }//end j for
    }//end i for
  }else if (bf_fl_mapping==bf_fl_mapping_list){ //V1.2a Mapping by list order

    if (n_phase!=n_fluor){
      printf("The number of elements in %s (%i) and %s (%i) \n",bright_list_file,n_phase,fluor_list_file,n_fluor);
      printf("must be equal. Make sure not lo leave any blank spaces at the end of the files");
      return 3;  // Exit code for bad parameters.
    }

    //Assining null time values
    ph_dmin=0.0;
    fl_dmin=0.0;
    for(i=0;i<n_phase;i++){
      ph_d[i]=0.0;
      ph_t[i]=0.0;
      fl_d[i]=0.0;
      fl_t[i]=0.0;
    }

    third_dmin=0.0;
    for(i=0;i<n_third;i++){
      third_t[i]=0.0;
      third_d[i]=0.0;
    }

    flat_dmin=0.0;
    for(i=0;i<n_flat;i++){
      flat_t[i]=0.0;
      flat_d[i]=0.0;
    }
  }

  void free(void *ptr); // rcell2: declaration for replacing g_free
  // https://developer.gnome.org/glib/stable/glib-Memory-Allocation.html#g-free
  // https://github.com/GNOME/glib/blob/master/glib/gmem.c
  // https://developer.gnome.org/glib/stable/glib-Basic-Types.html#gpointer
  // https://stackoverflow.com/a/20297598/11524079

  //Do a comparison of names to see if we should give the different
  //fluorescence files different flags.
  //Assume a three-letter prefix on the names identifies what the flag
  //value should be.
  flag[0]=0;
  for(i=1;i<n_fluor;i++){
    //V1.4a file names can have paths
    file_basename= basename(fluor_files[i]);
    c0=(file_basename)[0];
    c1=(file_basename)[1];
    c2=(file_basename)[2];
    // free(file_basename);  // rcell2: g_free replacement, not necessary
                             // https://stackoverflow.com/a/20297598/11524079
    flag[i]=flag[i-1]+1;     //Default to new flag

    //Look for a match among previous files
    // rcell2: we cant disable this for multi-z images, since time courses would be affected.
    // rcell2: the only choice is to increase the comparison up to the first N characters.
    for(j=0;j<i;j++){
            file_basename= basename(fluor_files[j]);
        if (
         (((file_basename)[0])==c0)&&
         (((file_basename)[1])==c1)&&
         (((file_basename)[2])==c2)
         ){
           flag[i]=flag[j]; //Found a match
           break;
      }
        // free(file_basename);  // rcell2: g_free replacement, not necessary
                                 // https://stackoverflow.com/a/20297598/11524079
        }
  }
  //Print out message if we have different flags set.
  j=0;
  for(i=0;i<n_fluor;i++){
    if (flag[i]!=0){
      j=1;
      break;
    }
  }
  if (j!=0){
    for(i=0;i<n_fluor;i++){
      printf("File %s is given flag number %i\n",
         fluor_files[i],flag[i]);
    }
  }

  //Write out the absolute time of the first file.  This is so later
  //we can correct for the differences in t0 from position to position.
  //To enable this function, pass "-z" flag argument to the CellID command-line call.
  if(time_of_t0_flag == 1){
    if( (fp=fopen("time_of_t0.txt","a")) == NULL ){
      printf("Couldnt open file time_of_t0.txt.\n");
    }else{
      dt=((fl_d[0]-fl_dmin)*seconds_per_day)+(fl_t[0]);
      fprintf(fp,"%i\n",dt);
      fclose(fp);
    }
  }

  //Loop over each of the fluorescence files to calculate the
  //fluorescence.  We'll search for cells using the phase file that
  //is closest in time to the fluorescence file
  for(i=0;i<n_fluor;i++){

    printf("----------------------------------------------------\n");
    printf("New Fluorescence image: %s.\n",fluor_files[i]);
    printf("----------------------------------------------------\n");
    fflush(stdout);  // https://stackoverflow.com/a/9469827
    
    //Load dark image to subtract for this exposure time
    free(fl);
    free(dark);
    dark=NULL;
    // Use dark correction on this FL image 
    // only if there is a dark image 
    // with the same exposure time
    // rcell2 TO-DO: this should simply obey the dark file list.
    t_exposure=get_exposure(fluor_files[i]);
    printf("Current fluor file exposure time: %f\n", t_exposure);
    for(j=0;j<n_dark;j++){
      //if (fabs(t_exposure-exposure[j])<0.001){
      printf("Checking dark file %s with exposure time: %f\n", basename(dark_files[j]), exposure[j]);
      if (t_exposure==exposure[j]){
        printf("Correcting with dark image %s\n",dark_files[j]); fflush(stdout);
        dark=get_data_from_tif_file(dark_files[j],0,NULL,&xmax_new,&ymax_new);
        break;
      } else {
        printf("Skipping dark correction with file %s. Exposure times were not equal.\n",dark_files[j]);
      }
    }

    if((fl=get_data_from_tif_file(fluor_files[i],0,dark,
                              &xmax_new,&ymax_new))==NULL){
      printf("Couldn't open tif file %s.\n",fluor_files[i]);
      return 4;  // Exit code for problems opening an image file.
    }
    if (((xmax>0)&&(xmax!=xmax_new))||((ymax>0)&&(ymax!=ymax_new))){
      printf("New file has different dimensions that others\n");
      printf("that were already loaded. (%i,%i) is not (%i,%i)\n",
             xmax,ymax,xmax_new,ymax_new);
      free(fl);
      return 6;  // Exit code for problems with an image's dimensions.
    }
    if ((xmax_new<=0)||(ymax_new<=0)){
      printf("Couldn't get data from tif file %s\n",fluor_files[i]);
      free(fl);
      return 4;  // Exit code for problems opening an image file.
    }
    xmax=xmax_new; //xmax<0 means haven't done yet, so just redefine
    ymax=ymax_new;//here even though usually the same

    if (bf_fl_mapping==bf_fl_mapping_time){
      JulianToYMD(fl_d[i],&year,&month,&day);
      Int_To_Hours_Minutes_Seconds(fl_t[i],&hours,&minutes,&seconds);
      printf("The time stamp of this file is: ");
      date_stamp(year,month,day);
      printf(" at ");
      time_stamp(hours,minutes,seconds);
      printf(".\n");
    }

    d_cur=fl_d[i]; //Date and time of the current fluorescence file
    t_cur=fl_t[i];

    //Load a "should-be-flat" image to use to flatten out this image.
    //Choose the image which has the same three-letter prefix as
    //the current image and is also closest in time.
    j_min=-1;
    dt_min=-1;
    
    //rcell2: the following was bugged.
    // c0=(fluor_files[i])[0];
    // c1=(fluor_files[i])[1];
    // c2=(fluor_files[i])[2];
    
    //rcell2: replaced the bug by this:
    s3=fluor_files[i];
    s4=s3;
    while ((*s3)!='\0'){
      if((*s3)=='/'){
        s4=(s3+1);
      }
      s3++;
    }
    //rcell2: s4 now points to first location after last "/" in name (copied from below).
    //rcell2: save the first three charachters of the basename.
    c0=s4[0];
    c1=s4[1];
    c2=s4[2];
    
    for(j=0;j<n_flat;j++){
      //Only consider images with the current three-character prefix.
      //But don't assume that these images are necessarily the same
      //directory, so go to the end of the directory information in the
      //name.
      // printf("Checking if flat file '%s' can be used for the current image\n", flat_files[j]);
      
      s=flat_files[j];
      s2=s;
      while ((*s)!='\0'){
          if((*s)=='/'){
            s2=(s+1);
          }
          s++;
      }
      //s2 now points to first location after last "/" in name.
      
      //Make sure name is more than three letters.
      if ((*s2)=='\0') continue;
      if ((*(s2+1))=='\0') continue;
      if ((*(s2+2))=='\0') continue;
      if ((c0==s2[0])&&(c1==s2[1])&&(c2==s2[2])){
          // printf("Found a flattening image with matching channel prefix.\n");
          dt=abs((d_cur-flat_d[j])*seconds_per_day+(t_cur-flat_t[j]));
          //Note we converted time to seconds from milliseconds above
          if ((dt_min<0)||(dt<dt_min)){
            dt_min=dt;
            j_min=j;
          }
      }
    }
    if (j_min<0){
      printf("No suitable flattening correction files found, not doing any flattening correction.\n");
    }else{
      printf("Correcting with should-be-flat image %s\n", flat_files[j_min]);
      //Calculate new corrections
      if (flat_cur!=j_min){ //We haven't done this correction yet (note: flat_cur is initialized as "-1").
          flat_cur=j_min;
          
          //Remove directory information
          s=flat_files[flat_cur]; 
          s2=s;
          while ((*s)!='\0'){
            if((*s)=='/'){
              s2=(s+1);
            }
            s++;
          }
          //s2 now points to first location after last "/" in name.
          
          //See if correction was already written to disk
          strcpy(line,s2);
          strcat(line,".corrections.tif");
          free(flat_cors);
          if((flat_cors=get_data_from_tif_file(line,0,NULL,
                                    &xmax_new,&ymax_new))==NULL){
            //File isn't there, so calculate correction
            ftmp=get_data_from_tif_file(flat_files[flat_cur],0,dark,
                      &xmax_new,&ymax_new);
            if (ftmp==NULL){
            printf("Couldnt open corrections file: %s\n",flat_files[flat_cur]);
              goto skip_flat;
            }
            if ((xmax_new!=xmax)||(ymax_new!=ymax)){
              printf("Cant do corrections with this size image: (%i x %i)\n",
                   xmax_new,ymax_new);
              free(ftmp);
              goto skip_flat;
            }
          flat_cors=flatten_image(ftmp,xmax,ymax,return_correction_array,
                                                                x_and_y, nonlinear);
            free(ftmp);
            //Write out correction file in current-directory.
            //use "line" as defined just above
            printf("Writing corrections to output file %s.\n",line);
            //Scale corrections by 20000.0 for writing out
            for(j=0;j<(xmax*ymax);j++)(flat_cors[j])*=(20000.0);
            if(output_data_to_tif_file(line,
                                       flat_cors,
                                       xmax,
                                       ymax,
                                       NULL,  // do not write cell labels for "flat_cors"
                                       0,
                                       16,
                                       0,
                                       mask_output)==0){
              printf("Couldn't output data to tif file %s.\n",line);
            }
          }
          //If we just wrote out corrections, then we just scaled them
          //by 20000 and similarly if we just read them in then they were
          //already scaled by 20000.0.  Unscale them here.
          for(j=0;j<(xmax*ymax);j++)(flat_cors[j])/=20000.0;
      }
      //Apply correction
      for(j=0;j<(xmax*ymax);j++){
          fl[j]*=(flat_cors[j]);
      }

    }
    skip_flat:

    load_global_arrays(1,fl,NULL,xmax,ymax); //Add new fl image to
    //global arrays (won't use FL until after find_cells() unless we're
    //doing a FRET image since we use the fl array to do the alignment)
    //so add it here also (also do it below in case we re-align the
    //FL image) below.
    //--->Actually, load_global_arrays() just copies the pointer, so
    //it doesn't really matter if you do it again or not....


    //Find third-image from third-image list which is closest in time to
    //the fluorescence image. (If there's a third image (ie, n_third!=0)
    if (n_third>0){
      j_min=-1;
      dt_min=(int)(0x1ffffff);
      for(j=0;j<n_third;j++){
          dt=abs((d_cur-third_d[j])*seconds_per_day+(t_cur-third_t[j]));
          //Note we converted time to seconds from milliseconds above
          if(dt<dt_min){
            dt_min=dt;
            j_min=j;
          }
      }
      //If we haven't loaded this file already, then load it
      //-->11/29/04-->Hack to just take next third image

      if (third_cur<0){
          j_min=0;
      }else{
          j_min=third_cur+1;
      }


      if (j_min!=third_cur){
          third_cur=j_min;
          free(third_image); //get_data...() routine will allocate space
          //Load dark image to subtract for this exposure time
          free(dark);
          dark=NULL;
          t_exposure=get_exposure(third_files[third_cur]);
          for(j=0;j<n_dark;j++){
            if (t_exposure==exposure[j]){
              dark=get_data_from_tif_file(dark_files[j],0,NULL,&xmax_new,&ymax_new);
              printf("Correcting third image with dark image %s\n",dark_files[j]);
              break;
            }
          }
          if((third_image=get_data_from_tif_file(third_files[third_cur],0,dark,
                                                         &xmax_new,&ymax_new))==NULL){
            printf("Couldn't open tif file %s.\n",third_files[third_cur]);
            return 4;  // Exit code for problems opening an image file.
          }
          if ((xmax!=xmax_new)||(ymax!=ymax_new)){
            printf("Third file has different dimensions than others\n");
            printf("that were already loaded. (%i,%i) is not (%i,%i)\n",
                xmax,ymax,xmax_new,ymax_new);
            free(third_image);
            return 6;  // Exit code for problems with an image's dimensions.
          }
          //Subtract fluorescence image as a test--for vacuole type
          if(third_image_type==vacuole_label){
            nloop=xmax*ymax;
            tmp1=0.0;
            tmp2=0.0;
            for(iloop=0;iloop<nloop;iloop++){
              tmp1+=(third_image[iloop]);
              tmp2+=(fl[iloop]);
            }
            tmp1/=((float)nloop);
            tmp2/=((float)nloop);
            for(iloop=0;iloop<nloop;iloop++){
              tmp=fl[iloop]-tmp2;
              if (tmp>100.0){
                tmp=10.0;
                third_image[iloop]=(third_image[iloop]-tmp1)/tmp;
                if (third_image[iloop]<=0.0) third_image[iloop]=0.0;
              }else{
                third_image[iloop]=0.0;
              }
            }
          }

          recalculate_internal=1; //Must re-calculate centers
          output_third_image=1; //Haven't written it out yet
          //Make sure to add new array to global variables in segment.
          load_global_arrays(2,third_image,NULL,xmax,ymax);
      }
    }

    if(bf_fl_mapping==bf_fl_mapping_time){
      //Find brightfield file which is closest in time to the current
      //fluorescence image.
      j_min=-1;
      dt_min=(int)(0x1ffffff);
      if (image_type==metamorph_deconvolution){
        //asg--Assume a 1-1 correspondence between the two file sets for
        //the metamorph-deconvolution case (since the decon was probably done
        //off line).
        j_min=i;
      }else{
        for(j=0;j<n_phase;j++){
             dt=abs((d_cur-ph_d[j])*seconds_per_day+(t_cur-ph_t[j]));
            //Note we converted time to seconds from milliseconds above
            if(dt<dt_min){
              dt_min=dt;
              j_min=j;
            }
        }
      }
    } else {
      //Mapping bright field and fluorescence by order.
      j_min=i;
    }

    bf_index[i]=j_min; //keeping track of the bright field image use to detect
                       //the cells of the current fluorescent file.

    new_phase=0;//Whether have a new file or are using the previous
    if((first==0)||(strcmp(phase_files[j_min],phase_files[j_cur])!=0)){ //new file
      //Set the split-cell cuts differently for first image if so
      //requested
      if (first==0){
        if (max_d_over_s_cut_t0>-10.0){
          max_d_over_s_cut=max_d_over_s_cut_t0;
        }
        if (max_split_d_over_minor_t0>-10.0){
          max_split_d_over_minor=max_split_d_over_minor_t0;
        }
      }else{
        max_d_over_s_cut=max_d_over_s_cut_save;
        max_split_d_over_minor=max_split_d_over_minor_save;
      }

      //If "do_recombination" is enabled (default 0, disabled)
      //try doing recombination of previously found cells
      //based on the fluorescence images that we read in since the last
      //call to find_cells(). (If this is the first call it won't do
      //anything since it loops over n_known, which starts at 0).
      if (do_recombination==1){
        printf("Doing recombination check.\n");
        recombination_check(i_last_find_cell_call,
                            n_recomb_cuts,
                            recombination_cuts_type,
                            recombination_cuts_flag,
                            recombination_cuts);
        if (first!=0){
          //Make output with re-combined data
          load_global_arrays(3,NULL,bf_fl_labels,xmax,ymax);  // "just in case"
                                                              // con type==3, bf_fl_labels se _asocia_ con el array "d"
                                                              // bf_fl_labels is declared above as "int *bf_fl_labels=NULL;" and not used in other files

          memset(bf_fl_labels,0,(xmax*ymax*sizeof(int)));     // llenar bf_fl_labels con ceros, esto quizas haga que "d" tenga ceros también

          if(label_cells==1){                  // mask_mod
            add_cell_number_to_the_data(i-1);  // its argument is "int i_t", maybe its the "time"
          }
          add_boundary_points_to_data(NULL);   // only executed for enabled do_recombination (default disabled)

          //Write out the files
          strcpy(line,"COMBINE_");
          strcat(line,phase_files[j_cur]);
          strcat(line,".out.tif");
          printf("Writing found cells and data to output file %s.\n",line);
          if(output_data_to_tif_file(line,
                                     bf,
                                     xmax,
                                     ymax,
                                     bf_fl_labels,
                                     0,
                                     8,
                                     0,
                                     mask_output)==0){
            printf("Couldn't output data to tif file %s.\n",line);
          }
        }
      }

      new_phase=1;            // flag a "new brightfield image"
      recalculate_internal=1; //New cells so must re-do internal calculations
      j_cur=j_min;

      first=1;
      i_last_find_cell_call=i;//How far back to alter cell list with
      //recombination_check_cuts.  i is the loop over the fluorescence
      //images, and each i shows up in update_list_of_found_cells() below

      printf("Doing new cell search.\n");
      printf("----------------------------------------------------\n");
      printf("New brightfield image: %s.\n",phase_files[j_cur]);
      printf("----------------------------------------------------\n");
      fflush(stdout);
      
      if(bf_fl_mapping==bf_fl_mapping_time){
        JulianToYMD(ph_d[j_cur],&year,&month,&day);
        Int_To_Hours_Minutes_Seconds(ph_t[j_cur],&hours,&minutes,&seconds);
        printf("The time stamp of this file is: ");fflush(stdout);
        date_stamp(year,month,day);
        printf(" at ");
        time_stamp(hours,minutes,seconds);
        printf(".\n"); fflush(stdout);
      }

      //Read in brightfield image data for ith file
      free(bf);
      //Do no dark field correction for bright field image....
      if((bf=get_data_from_tif_file(phase_files[j_cur],0,NULL,&xmax_new,&ymax_new))==NULL){
        printf("Couldn't open tif file %s.\n",phase_files[j_cur]);
        free(bf);
        return 4;  // Exit code for problems opening an image file.
      }
      if ((xmax!=xmax_new)||(ymax!=ymax_new)){
              printf("BF file has different dimensions than others\n");
              printf("that were already loaded. (%i,%i) is not (%i,%i)\n",
                  xmax,ymax,xmax_new,ymax_new);
              free(bf);
              return 6;  // Exit code for problems with an image's dimensions.
      }

      //Make sure to add new array to global variables in segment.
      load_global_arrays(0,bf,NULL,xmax,ymax);

      //In case of aligning to bf, load the current bf image
      if (align_fl==2){
        align_image(bf,xmax,ymax,0); //Load bf to compare fl to
      }

      //Search over this file to get cells and fill boundary[] and
      //interior[] in segment.c.

      //Make the labels[] array (only used as work array in find_cells).
      //printf("labels\n");fflush(stdout);
      if (bf_fl_labels==NULL){
        bf_fl_labels=(int *)malloc(xmax*ymax*sizeof(int));  // declarar bf_fl_labels como un array para los pixeles
      }                                                     // ver: www.geeksforgeeks.org/dynamic-memory-allocation-in-c-using-malloc-calloc-free-and-realloc
      memset(bf_fl_labels,0,(xmax*ymax*sizeof(int)));       // luego lo llena de ceros (porque se inicializa con "junk")
      load_global_arrays(3,NULL,bf_fl_labels,xmax,ymax);    // con type==3, bf_fl_labels se _asocia_ con el array "d"

      //Now if we have a fret image, find the split between the
      //upper and lower using the current fluorescence image.
      if (fret_image==1){
         free(fret_label_array);
         fret_label_array=find_split_regions(fl,xmax,ymax); //using FL[] array
         //fret_labels() calculates top and bottom regions using fl[] array.
         load_global_arrays(4,NULL,fret_label_array,xmax,ymax);
      }else{
         load_global_arrays(4,NULL,NULL,xmax,ymax);//just in case to
         //force a crash if there's something weird going on.
      }

      //Make sure call load_global_arrays() first:
      n_found_cells=find_cells(&found_boundary_points,&found_interior_points);
      //We now have a list of interior and boundary points as well as the
      //mean x and mean y of the interior points. (These are saved in
      //global variables in find_cells().)

      //Free pixels, keeping only the latest that were recently put
      //in the cell lists (the cells just found in find_cells() haven't
      //been added to the list yet).

      if(out_mask!=1){                           // mask_mod: this must be messing up the TSV output
      	free_pixels_from_earlier_time_points();  // mask_mod: disabled on mask output option
      }
    }

    //Check if we want to align to BF image (or also to previous
    //fl image) to correct for bad repositioning.
    if ((align_fl==1)||(align_fl==2)){ //Re-align fl to first fl image
      align_image(fl,xmax,ymax,1); //First call will automatically
      //save fl for later alignments if hasn't been called yet
      if ((i==0)&&(align_fl==2)){ //Only compare first image to BF
          //align_image(fl,xmax,ymax,0);
          //printf("Only comparing first FL image to BF.\n");
      }

    }
    load_global_arrays(1,fl,NULL,xmax,ymax); // Make sure to add new fl[] array to global variables in segment.c
                                             // se pasa fl a v1, y despues fl=v1; o sea que fl = fl ?? no entiendo

    //Each cell may have wiggled a bit. If the correction option was
    //passed in parameters.txt then re-align individual cells with
    //the current fluorescence image.
    // mask_mod: This bit aligns/offsets points in many arrays such as boundary[i], for each "n_found"
    if (align_individual_cells==1){
      if (align_individual_cells_boundary==1){
        align_found_cells_to_fl(2); //2 to use boundary
      }else{
        align_found_cells_to_fl(0);
      }
    }

    memset(bf_fl_labels,0,(xmax*ymax*sizeof(int)));  // mask_mod: reset bf_fl_labels to 0

    add_boundary_points_to_data(NULL); // mask_mod: this is the key call to add boundaries to FL.out images

    //Check for nucleus or vacuole, etc, using third image.
    if ((third_image_type!=no_third_image)||(fret_image==1)){

      //Conditions to find internal structure, etc
      if (third_labels==NULL){
          third_labels=(int *)malloc(xmax*ymax*sizeof(int));
          memset(third_labels,0,(xmax*ymax*sizeof(int)));
      }

      //The third_labels[] may be used in calculate_fluorescence...
      //to determine which pixels to ignore.
      load_global_arrays(3,NULL,third_labels,xmax,ymax);  // con type==3, bf_fl_labels se _asocia_ con el array "d"
      memset(third_labels,0,(xmax*ymax*sizeof(int)));

      //If there is no third image, then we use the top part of the fret
      //fl image. Whether to use third_image or not is decided by existence
      //of third image in internal_structure().
      //In other words, if you take out the third image for a FRET image
      //it defaults to using the upper part of the image for the nucleus.
      if (third_image_type==no_third_image){
          recalculate_internal=1; //Always re-do if not using a third image
      }
      if(nuclear_fret_lower==0){
          internal_structure(0,1); //Use lower part for fret image
      }else{
          internal_structure(1,1); //Use upper part for fret image
      }
    }else{
      //do internal_structure() always using FL image

      recalculate_internal=1; //Re-do search for centers
      if (force_nucleus_in_center==1){
          internal_structure(0,0);
      }else{
          internal_structure(0,1);
      }

    }

    //Calculate how much fluorescence is within each cell.  (This will
    //fill an array in segment.c.)
    //Do this after internal_structure() since, for the vacuole-third-image
    //types, we use the third_labels[] array to remove pixels associated
    //with the vacuole.
    //Note that the third_labels[] arrays gets sent to calculate_flu...
    //through the load_global_arrays() subroutine.

    calculate_fluorescence_with_r_info();

    background_level(i);  //Calculate mode of fl[] array of pixels not
                          //included in any found cell

    next_prev_fl_comparison();  //Collect information by comparing
                                //pixel-by-pixel the found cells and the same
                                //pixel locations in the previous image.

    //Add the cells from the latest file to the list of known cells using time
    //index i, and also the number of seconds since the first image
    dt=(d_cur-fl_d[0])*seconds_per_day+(t_cur-fl_t[0]);
    update_list_of_found_cells(i,dt,flag[i]);

    //time_index[i] is the time-integer for this flag value.
    //E.g. if this is the 2nd time we've seen this flag value,
    //then time_index[i] will be set to 1 (it starts at 0).
    time_index[i]=time_flag[flag[i]];
    time_flag[flag[i]]++;

    //Only output fluor_files[] if it isn't actually the brightfield image.
    if ((treat_brightfield_as_fluorescence==0)||(flag_bf[i]==0)){
      //Write out the files

      if (output_individual_cells==1){
          strcpy(line,"cells/");
          strcat(line,fluor_files[i]);
          if(output_individual_cells_to_file(i,line,fl,xmax,ymax,0,8,0,mask_output)==0){
            printf("Couldn't output individual to tif file %s.\n",line);
          }
      }

      strcpy(line,fluor_files[i]);
      strcat(line,".out.tif");
      printf("Writing found cells and data to output file %s.\n",line);

      if(output_data_to_tif_file(line,
                                 fl,
                                 xmax,
                                 ymax,
                                 bf_fl_labels,
                                 1,
                                 8,
                                 0,
                                 mask_output)==0){
          printf("Couldn't output data to tif file %s.\n",line);
        }
  }

  if (output_third_image==1){
    if (output_individual_cells==1){
        strcpy(line,"cells/");
        strcat(line,third_files[third_cur]);
        if(output_individual_cells_to_file(i,
                                           line,
                                           third_image,
                                           xmax,
                                           ymax,
                                           0,
                                           8,
                                           0,
                                           mask_output)==0){

        printf("Couldn't output individual to tif file %s.\n",line);
        }
    }
    output_third_image=0;
    strcpy(line,third_files[third_cur]);
    strcat(line,".out.tif");
    printf("Writing found cells and data to output file %s.\n",line);
    if(output_data_to_tif_file(line,
                               third_image,
                               xmax,ymax,
                               third_labels,
                               2,
                               8,
                               0,
                               mask_output)==0){

      printf("Couldn't output data to tif file %s.\n",line);
    }
  }

  if (align_individual_cells==1){
    align_found_cells_to_fl(1); //Re-set the BF boundaries, etc
  }

  if(new_phase==1){ //Haven't written out yet
    //We do it down here so we can write it out after we've added
    //the numbers labelling the cells.  We use the updated list of
    //known cells for the numbers, so we need to have run over the
    //fluorescence stuff already.
    //Label each cell in the tiff file with a number

    load_global_arrays(3,NULL,bf_fl_labels,xmax,ymax);  //just in case
                                                        // con type==3, bf_fl_labels se _asocia_ con el array "d"

    memset(bf_fl_labels,0,(xmax*ymax*sizeof(int)));


    // mask_mod:
    if(mask_output>0){
      add_cell_mask_data(NULL, i, fill_interior, label_cells, interior_offset);
    }
    else{
      add_cell_number_to_the_data(i);  // its argument is "int i_t"
      add_boundary_points_to_data(NULL);
    }

    if (output_individual_cells==1){
      //Write out the files
      strcpy(line,"cells/");
      strcat(line,phase_files[j_cur]);
      if(output_individual_cells_to_file(i,line,bf,xmax,ymax,0,8,0,mask_output)==0){
        printf("Couldn't output individual to tif file %s.\n",line);
      }
    }

    //Write out the files
    strcpy(line,phase_files[j_cur]);
    strcat(line,".out.tif");
    //strcat(line,".out_cfp.tif");
    printf("Writing found cells and data to output file %s.\n",line);
    if(output_data_to_tif_file(line,
                               bf,
                               xmax,
                               ymax,
                               bf_fl_labels,
                               0,
                               8,
                               0,
                               mask_output)==0){

      printf("Couldn't output data to tif file %s.\n",line);
    }

  }

  }//end loop over the number of fluorescence files

  //We're done the loop over the fluorescence image. We've been calling
  //the recombinatino_check_cuts() before each find_cell() call to account
  //for the previous call+fluorescence images. So we have to do it one
  //more time for these last images.
  if (do_recombination==1){

    recombination_check(i_last_find_cell_call,
                        n_recomb_cuts,
                        recombination_cuts_type,
                        recombination_cuts_flag,
                        recombination_cuts);

    //Make output with re-combined data
    load_global_arrays(3,NULL,bf_fl_labels,xmax,ymax);  //just in case
                                                        // con type==3, bf_fl_labels se _asocia_ con el array "d"
    memset(bf_fl_labels,0,(xmax*ymax*sizeof(int)));

    if(label_cells==1){                  // mask_mod: optional labeling through "-l" command line option
      add_cell_number_to_the_data(i-1);  // mask_mod: its argument is "int i_t"
    }
    add_boundary_points_to_data(NULL);   // mask_mod: this chunk is only executed for enabled do_recombination (default disabled)

    //Write out the files
    strcpy(line,"COMBINE_");
    strcat(line,phase_files[j_cur]);
    strcat(line,".out.tif");
    //strcat(line,".out_cfp.tif");
    printf("Writing found cells and data to output file %s.\n",line);
    if(output_data_to_tif_file(line,
                               bf,
                               xmax,
                               ymax,
                               bf_fl_labels,
                               0,
                               8,
                               0,
                               mask_output)==0){
      printf("Couldn't output data to tif file %s.\n",line);
    }
  }

  //Now write out the list of cells
  if(output_basename==NULL){
    //The basename of the output file was passed in as an
    //argument.  It should contain the directory path of where to put
    //the output files.
    output_basename="Output";
  }

  printf("Writing output to files to directory %s.\n",output_basename);

  //if(output_cells(line2,append,time_index)==0){
  if(paw_output==1){
    printf("Output file in PAW format.\n");
    if(output_cells(output_basename,append,time_index)==0){
      printf("Couldn't open X output files: %s.\n",line2);
    }
  } else {
    printf("Output file in R format.\n");
    if(output_cells_single_file(output_basename,append,time_index, out_mask)==0){
      printf("Couldn't open X output files: %s.\n",line2);
    }
  }

  //V1.2a Creating brightfield - fluorescent mapping file
  strcpy(bf_fl_file_name, output_basename);
  strcat(bf_fl_file_name, "_bf_fl_mapping");

  if((bf_fl_file=fopen(bf_fl_file_name,"w"))==NULL){
    printf("Couldn't open file %s\n",bf_fl_file_name);
    fflush(stdout);
    return 7;  // Exit code for problems opening a file.
  }

  fprintf(bf_fl_file,"fluor\tflag\tt.frame\tbright\tbf.as.fl\n");
  for(i=0;i<n_fluor;i++){
    fprintf(bf_fl_file,"%s\t%i\t%i\t%s\t%i\n",fluor_files[i],flag[i],
                          time_index[i],phase_files[bf_index[i]],flag_bf[i]);
  }

  if (bf_fl_file!=NULL)fclose(bf_fl_file);

  printf("\n\n******** CellID is done! enjoy your single cell datasets :) ********\n");
  fflush(stdout);
  return 0;
}
