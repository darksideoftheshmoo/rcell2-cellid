#' A function to download the latest workflow template in Rmarkdown
#' 
#' Will download the .Rmd file to the current working directory.
#' 
#' @param file_name File name for the downloaded workflow template.
#' @param open.template Try using file.edit to open the file in RStudio after getting it.
#' 
#' @export
get_workflow_template_cellid <- function(
    file_name = "rcell2.cellid_workflow_template.Rmd",
    open.template = T){
  
  if(file.exists(file_name)) stop("get_workflow_template error: file", file_name, "exists.")
  
  workflow.file <- system.file(
    "rmarkdown/templates/rmd_template.cellid/skeleton/skeleton.Rmd",
    package = "rcell2.cellid"
  )

  if(file.exists(workflow.file)){
    file.copy(from = workflow.file, to = file_name)
  } else {
    download.file(url = paste0("https://raw.githubusercontent.com/darksideoftheshmoo/rcell2-cellid/main/",
                               "inst/rmarkdown/templates/rmd_template.cellid/skeleton/skeleton.Rmd"), 
                  destfile = file_name)
  }
  
  if(open.template){
    file.edit(file_name)
  }
}
