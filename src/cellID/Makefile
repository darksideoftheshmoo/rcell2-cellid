# gcc -Wall enables all compiler's warning messages.
CFLAGS = -Wall -Wno-unused-variable -Wno-unused-but-set-variable

# Guia para incluir cosas:
## -I <include directory>
## -L <library directory>
## -l <library name>
## The include directory will be location of the header files. 
## The library directory will be the location of the library (.a or .so)
## and the library name will be the name of the library file, without the leading 'lib' prefix or its extension (i.e. -lsal rather than -l libsal.a ).

# Librerías de C para CellID
## Para cellid hace falta: "-lm" y ahora que comprime los TSV necesita "-lz"
CLIBScellid = -lm -lz

# Conditional external -ltiff test. Check result from Makevars exported environment:
$(info "Makefile got tiflibs < $(tiflibs) > with tiflinks < $(tiflinks) > from environment.")

# Check test result and define appropriate links to the libtiff library.
# Note: tiflibs es una variable exportada del Makevars.
ifeq (tiflibs, ${tiflibs})
$(info "Using bundled libtiff")
# Ahora hay que hacer los links necesarios, de tipo "dinamic" o "static".

# Para dynamic linking, los requerimients se consultan con: $  pkg-config --libs libtiff-4
## En Arch solo hace falta "-ltiff".
## No uso este approach porque requeriría dejar el folder con libtiff.

# Para incluir libtiff "estatico" (el libtiff.a en vez del libtiff.so) tengo que saber que
## otras flags "-l" necesito
## Eso se puede consultar con: $ pkg-config --static --libs libtiff-4
## Luego para linkear libtiff "estatico" hay que poner dos cosas:
## 1) Linkear el path donde está el archivo "libtiff.a" así: "-L../../inst/tiff4/r_build/libtiff"
## 2) y luego la librería usando: "-ltiff"
## Además, según pkg-config, para linking necesito también: "-llzma -ljbig -ljpeg -lz -lm",
## y parece que también necesitaba: "-lzstd" y "-lwebp",
## pero en Windows alcanzó con "-ltiff -ljpeg -lz".

# Link to the static "libtiff.a".
# Note: tiflinks es una variable exportada del Makevars.
# It may contain just "-ltiff" (to link the system tiff library),
# or "-ltiff -lz -lm"  (to link the bundled tiff library),
CLIBStiff = -L../../inst/tiff4/lib ${tiflinks}
# We also need to include the ".h" header files
CLIBStiffh = -I../../inst/tiff4/include

else
$(info "Using system's libtiff")
# Simple link to system's libtiff:
CLIBStiff = ${tiflinks}
# No header files are needed.
# However "tiffio.h" was not found in MacOS.
# Added "-I/usr/local/include" as possible fix:
CLIBStiffh = -I/usr/local/include
endif

# Una variable con los targets-requisitos de cell,
# que también se pasa a "gcc" al compilar ese target ("cell"),
# y aparece en el target "clean" para limpiar el output.
objects = cell.o segment.o tif.o nums.o date_and_time.o fit.o fft.o fft_stats.o split_and_overlap.o contiguous.o fl_dist.o align_image.o flatten.o

# El orden de los argumentos en gcc es importante.
## Las libs deben ir despues de los objetos .c/.o/... o se descartan
## Leer:  https://stackoverflow.com/questions/2624238/c-undefined-references-with-static-library
cell : $(objects)
	echo Building Cell-ID with LD_LIBRARY_PATH: $(LD_LIBRARY_PATH)
	gcc -o cell $(CFLAGS) $(objects) $(CLIBStiff) $(CLIBScellid)

# Target original
#cell : $(objects)
#	gcc -o cell $(CFLAGS) -lm -l$(tiflib) $(objects)

cell.o: segment.h tif_routines.h date_and_time.h nums.h point.h image_type.h align_image.h split_and_overlap.h parameters.h

tif.o: tif_routines.h
	gcc $(CFLAGS) -c -o tif.o tif.c $(CLIBStiffh) $(CLIBStiff)  # Not sure if this is important, and the "-c" in theory ignores the "-l".

segment.o: segment.h nums.h fit.h tif_routines.h fft_stats.h point.h image_type.h split_and_overlap.h contiguous.h fl_dist.h parameters.h flatten.h

contiguous.o: contiguous.h

nums.o: nums.h tif_routines.h

date_and_time.o: date_and_time.h

fit.o: fit.h

fft.o: fft.h complex.h

fft_stats.o: segment.h fft.h fft_stats.h point.h complex.h

split_and_overlap.o: split_and_overlap.h contiguous.h

fl_dist.o: fl_dist.h contiguous.h

align_image.o: align_image.h

flatten.o: flatten.h

clean:
	rm -f $(objects)

