#' Manhattan plot with all tested snarls
#' @description Generate a Manhattan plot showing the p-values across the pangenome.
#'
#' @importFrom ggplot2 ggplot aes geom_hline scale_x_continuous scale_size_continuous labs theme_bw theme element_text ggsave facet_grid guides scale_shape_manual element_blank
#' @importFrom utils read.table head write.table
#' @importFrom rlang .data
#'
#' @param assoc Either the data.frame imported by *import_assoc*, or the path to STOAT's output (*assoc.pvalues.tsv.gz)
#' @param output_file If not NULL, the name of the output image where to save the plot (image type guessed from the file name).
#' @param p_threshold P-value threshold for the horizontal significance line (default: 1e-5).
#' @param output_file If not NULL, the name of the output image where to save the plot (image type guessed from the file name).
#' @param chr_order list of chromosome names in the desired order for the plot. If NULL, will try to guess.
#' @param show_all_points should we plot all the points? Default is FALSE (cluster close-by points instead). TRUE is slower (and heavier in PDF) for large genomes.
#' @param show_top_points how many top associations should be shown as single points? Default is 0. 
#' @param wrap_by_rows how many rows to wrap the chromosome panels. Default: 1 (no wrapping)
#' @param no_plotting if TRUE, no graph is plotted and the ggplot/gtable is just returned. Only useful in very specific situations, with multi-row plots and where you want to handle the plotting yourself.
#'
#' @return A ggplot object. Saves a file if output_file is provided.
#' @name manhattan_plot
#' @export

manhattan_plot <- function(assoc, p_threshold=1e-5, output_file=NULL, chr_order=NULL,
                           show_all_points=FALSE, show_top_points=0, wrap_by_rows=1,
                           no_plotting=FALSE) {

  ## potentially load the association results
  if(is.character(assoc) && length(assoc) == 1) {
    assoc = import_assoc(assoc)
  }

  ## which column to use for the pvalue?
  pv_col = ifelse(any(colnames(assoc) == 'P'), 'P', 'P_CHI2')

  ## handle null pvalues
  null_pvs = which(assoc[, pv_col] == 0)
  if (length(null_pvs) > 0) {
    min.pv = min(assoc[-null_pvs, pv_col])
    min.pv = 10^(floor(log10(min.pv) - 5))
    warning(length(null_pvs), ' p-values changed from 0 to ', min.pv)
    assoc[null_pvs, pv_col] = min.pv
  }

  ## prepare -log10 pvalues
  assoc$logp = -log10(assoc[,pv_col, TRUE])
  
  ## set the chromosome order
  if (is.null(chr_order)) {
    ## try to guess
    exhaustive_suffix = c(1:100, 'X', 'Z', 'Y', 'W', 'O', '0', 'U', 'V', 'M', 'MT')
    exhaustive_prefix = c('', 'chr', 'Chr', 'CHR')
    exhaustive_order = unlist(lapply(exhaustive_prefix,
                                     function(pref) paste0(pref, exhaustive_suffix)))
    ## keep the ones in common
    all_chrs = unique(assoc$CHR)
    chr_order = intersect(exhaustive_order, all_chrs)
    ## eventually add other ones in alphabetical order
    if(!all(all_chrs %in% chr_order)){
      chr_order = c(chr_order, sort(setdiff(all_chrs, chr_order)))
    }
  }
  assoc$CHR = factor(assoc$CHR, levels=chr_order)

  ## plot
  ggp <- ggplot()

  if (!show_all_points & nrow(assoc) > 1e6) {
    ## prepare an approximate location/pvalue to group points
    max.logp = max(assoc$logp)
    logp.bks = seq(-1, max.logp+1, by=min(max.logp, 10) / 100)
    logq.bks.vals = (logp.bks[-1] + logp.bks[-length(logp.bks)] ) / 2
    assoc$logp.a = logq.bks.vals[as.numeric(cut(assoc$logp, logp.bks))]

    pos.bks = seq(-1, max(assoc$START_OFFSET)+1, length.out=50)
    pos.bks.vals = (pos.bks[-1] + pos.bks[-length(pos.bks)] ) / 2
    assoc$pos.a = pos.bks.vals[as.numeric(cut(assoc$START_OFFSET, pos.bks))]

    if (show_top_points > 0) {
      ## force top associations as single points
      assoc$logp.a = ifelse(rank(assoc$P, ties.method='min') > show_top_points,
                            assoc$logp.a, assoc$logp)
    }
    
    assoc = dplyr::summarize(dplyr::group_by(assoc, .data$CHR, .data$pos.a, .data$logp.a),
                             logp=mean(.data$logp),
                               nvariants=dplyr::n(), 
                             START_OFFSET=mean(.data$START_OFFSET),
                             .groups='drop')    
  }

  ## function to make a manhattan plot for a set of chromosomes
  make_plot <- function(chrs.names=NULL){
    if (!is.null(chrs.names)) {
      assoc = assoc[which(assoc$CHR %in% chrs.names),]
    }
    ggp = ggplot()
    if (any(colnames(assoc) == 'nvariants')){
      ## points were clustered
      ggp = ggplot(assoc, aes(x=.data$START_OFFSET/1e6, y=.data$logp,
                              size=.data$nvariants, shape=.data$nvariants>1)) +
        geom_point(alpha=.3) +
        scale_size_continuous(range=c(1,5)) + guides(size='none') +
        scale_shape_manual(values=c(17, 16), name='variant', labels=c('single', 'cluster'))
    } else {
      ## show all points
      ggp = ggplot(assoc, aes(x=.data$START_OFFSET/1e6, y=.data$logp)) +
        geom_point(alpha=.3)  
    }
    ## add panels, horizontal line, legends, etc
    ggp + facet_grid(.~.data$CHR, scales='free', space='free') +
      scale_x_continuous(n.breaks=3) + 
      geom_hline(yintercept=-log10(p_threshold),
                 color="red", linetype="dashed") +
      labs(x='position (Mbp)', y=expression(-log[10](P))) +
      theme_bw() + 
      theme(strip.text.x=element_text(angle=90),
            panel.spacing.x=grid::unit(0, 'cm'),
            legend.position='bottom')
  }

  ## either make multiple rows or just one
  if (wrap_by_rows > 1){

    ## compute the cumulative chromosome length and cut in rows of
    ## about equal sizes
    chrs.len = dplyr::summarize(dplyr::group_by(assoc, .data$CHR),
                                size=max(.data$START_OFFSET)/1e6)
    chrs.len$cum_size = cumsum(chrs.len$size)
    chrs.len$row = cut(chrs.len$cum_size, breaks=seq(0, max(chrs.len$cum_size), length.out=wrap_by_rows+1))

    ggp.l = lapply(1:nlevels(chrs.len$row), function(row_lvl_idx) {
      cur_chrs = chrs.len$CHR[which(as.numeric(chrs.len$row) == row_lvl_idx)]      
      ggp = make_plot(cur_chrs)
      ## remove the title of the x-axis and legend, except for the last row
      if (row_lvl_idx != nlevels(chrs.len$row)) {
        ggp = ggp + theme(axis.title.x=element_blank())
        if (any(colnames(assoc) == 'nvariants')) {
          ggp = ggp + theme(legend.position='none')
        }
      }
      return(ggp)
    })

    ## combine in rows
    row.h = rep(1, length(ggp.l))
    row.h[length(row.h)] = 1 + wrap_by_rows / 10
    ggp = gridExtra::arrangeGrob(grobs=ggp.l, heights=row.h)

    if (!no_plotting) {
      ## plot the multi-panel plot
      gridExtra::grid.arrange(ggp)
    }
    
  } else {
    ## no wrapping, just one row
    ggp = make_plot()
  }

  if(!is.null(output_file)) {
    ggsave(output_file, plot=ggp, width=12, height=4)
  }

  return(ggp)
}
