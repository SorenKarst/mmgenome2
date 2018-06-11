#' @title Visualise metagenomes contigs in circos plot
#' 
#' @import circlize
#' @import dplyr
#' @importFrom Biostrings names
#' @importFrom Biostrings width
#' @importFrom Biostrings length
#' @importFrom Biostrings letterFrequencyInSlidingView
#' @importFrom RColorBrewer brewer.pal
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
mmcircos <- function(mm,
                     assembly,
                     coverage_data,
                     link_data,
                     window_plotsize = 1000, 
                     label_size = 1,
                     coverage_height = 0.1,
                     assembly_height = 0.05
                     ) {
  # Checks and error messages before anything else
  if(nrow(mm) > 500)
    stop("Number of scaffolds in data subset is too,
         high for meaningful plotting in circos. Reduce
         to 500 or fewer scaffolds", call. = FALSE)
  
  # Custom function
  bins <- function(position, step){
    rgn <- unique(c(seq(0, max(position), step), max(position)))
    lab <- rgn[-length(rgn)] + (rgn[-1] - rgn[-length(rgn)])/2
    cut(position, rgn, lab, include.lowest = T) %>%
      as.character() %>%
      as.numeric()
  }
  
  # Import data
  a <- assembly[as.character(mm$scaffold)]
  c <- coverage_data %>%
    filter(scaffold %in% mm$scaffold)
  l <- link_data %>%
    filter(scaffold %in% mm$scaffold)
  w <- window_plotsize
  
  # Format assembly data for plotting
  ad <- data_frame(
    name = Biostrings::names(a),
    start = 0,
    end = Biostrings::width(a)
  )

  ad_col <- rep(RColorBrewer::brewer.pal(12, "Paired"),
                ceiling(nrow(ad)/12))
  
  # Calculate GC data for plotting
  gc <- lapply(a, function(x){
    if (Biostrings::length(x) >= w) {window = w} else {window = Biostrings::length(x)}
    x %>%
      Biostrings::letterFrequencyInSlidingView(., window, c("G", "C")) %>%
      as.data.frame() %>%
      transmute(pos = row_number(), gc = rowSums(.)/window) %>%
      filter(row_number() %% window == 1)}) %>%
    bind_rows(.id = "name") %>%
    transmute(name, start = pos, end = pos, gc)
  
  # Prepare coverage data
  cd <- c %>%
    group_by(scaffold) %>%
    mutate(bin = bins(position, w)) %>%
    group_by(scaffold, bin) %>%
    summarise(coverage = mean(coverage)) %>%
    ungroup() %>%
    transmute(name = scaffold, start = bin, end = bin, coverage = coverage)
  
  # Prepare links data
  ld <-  l %>%
    transmute(rname_1 = as.character(scaffold1),
              rname_2 = as.character(scaffold2),
              pos_1 = pos1,
              pos_2 = pos2,
              connections) %>%
    as.data.frame()
  
  # Build plot
  circos.par(gap.after = 0.5,
             points.overflow.warning = F,
             start.degree = 90)
  circos.genomicInitialize(ad, plotType = NULL)
  
  # Assembly labels
  circos.track(ylim = c(0,05),
               panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter,
                             CELL_META$ycenter,
                             CELL_META$sector.index, 
                             facing = "clockwise",
                             niceFacing = T,
                             cex = label_size)},
               track.height = 0.05,
               bg.border = NA)
  # Coverage
  circos.genomicTrack(cd, ylim = c(min(cd$coverage), max(cd$coverage)),
                      panel.fun = function(region, value, ...) {
                        xlim = CELL_META$xlim
                        circos.genomicLines(region,
                                            value,
                                            col = "lightgrey",
                                            lwd = 0.01,
                                            area = T)},
                      track.height = coverage_height,
                      bg.border = NA)
  
  # GC
  circos.genomicTrack(gc, ylim = c(min(gc$gc), max(gc$gc)), track.index = 2,
                      panel.fun = function(region, value, ...) {
                        xlim = CELL_META$xlim
                        circos.genomicLines(region,
                                            value,
                                            col = rgb(0.46, 0.77, 0.47),
                                            lwd = 0.01)},
                      track.height = coverage_height,
                      bg.border = NA)
  # Assembly
  circos.track(ylim = c(0, 1),
               cell.padding = c(0.02, 1.00, 0.005, 1.00),
               panel.fun = function(x, y) {
                 index = CELL_META$sector.numeric.index
                 xlim = CELL_META$xlim
                 circos.rect(xlim[1],0, xlim[2], 1,
                             col = acol[index],
                             lwd = 0.01)
                 circos.axis(h = "top",
                             major.at = seq(0, xlim[2], 25000),
                             major.tick.length = 0.2,
                             minor.ticks = 0,
                             labels = NULL,
                             lwd = 0.01)
               },
               track.height = assembly_height,
               bg.border = NA)
  
  # Links
  circos.genomicLink(ld[,c("rname_1", "pos_1", "pos_1")],
                     ld[,c("rname_2", "pos_2", "pos_2")],
                     col = rgb(0.42, 0.68, 0.84, ld$n/max(ld$n)),
                     lwd = ld$n/max(ld$n) * 5,
                     border = NA)
  
  # Finish
  circos.clear()
  ```
}