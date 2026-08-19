#' General a summary report.
#' @description An HTML file will be created in the working directory.
#' @param assoc_file path to the association file
#' @param output_html_file path for the output HTML file
#' @param force_overwrite should existing Rmd/HTML output file be overwritten?
#' @return a character with the report filename
#' @name make_report
#' @export
make_report <- function(assoc_file, output_html_file='stoat.summary.report.html',
                        force_overwrite=FALSE) {
  ## prepare paths for new files
  output_rmd_file = gsub('\\.html$', '.Rmd', output_html_file)

  if (file.exists(output_rmd_file) & !force_overwrite) {
    stop(paste(output_rmd_file, 'exists. Run with force_overwrite=TRUE to overwrite it.'))
  }
  if (file.exists(output_html_file) & !force_overwrite) {
    stop(paste(output_html_file, 'exists. Run with force_overwrite=TRUE to overwrite it.'))
  }
  
  work.dir = getwd()

  ## read and edit report content
  report.src = file.path(system.file("reports", package = "StoatPlot"), 'stoat.summary.report.Rmd')
  rmd.content = scan(report.src, '', quiet=TRUE, sep='\n', blank.lines.skip=FALSE)
  rmd.content = gsub('\\{\\{assoc_file\\}\\}', paste0('"', assoc_file, '"'), rmd.content)

  ## write Rmd report
  write(rmd.content, file=output_rmd_file, sep='\n')

  ## compile report
  rmarkdown::render(output_rmd_file, knit_root_dir=work.dir,
                    output_file=basename(output_html_file))

  return(output_html_file)
}


