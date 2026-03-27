
# region colors

geo.colours <- c(
  "South America"   = "#E41A1C",
  "Central America" = "#FFAA10", 
  "Caribbean"       = "#FFD92F",
  "Micronesia"      = "#A99060",
  "Hawaii"          = "#66C2A5",
  "Indonesia"       = "#A2D2E2",
  "Australia"       = "#002FA7",  
  "Taiwan"          = "#710193",
  "Africa"          = "#FA9591",
  "unknown"         = "grey"
)


lineage_colors <- c(
  LAC = "#FFCA10",
  Tw1 = "#161308",
  Tw2 = "#FF66CC",
  Indo1 = "#3390FF",
  Tw3 = "#D1007A",
  Indo2 = "#4EA3C8",
  Tw4 = "#AA3399",
  Mic1 = "#7A541F",
  Mic2 = "#A99020",
  HC = "#D9EF8B",
  Hw1 = "#1B4D1B",
  Tw5 = "#6600CC",
  Au = "#002FA7",
  Tw6 = "#A366FF",
  Hw2 = "#4C8C4C",
  Af = "#F4B6B6",
  Indo3 = "#C8E6F0",
  Tw7 = "#D9B3FF",
  Hw3 = "#00FF00"
)


species_color <- c("C. elegans"   = "#DB6333",
                   "C. tropicalis" = "#0719BC",
                   "C. briggsae"   = "#53886C")


genome_domain_colors <- c("Tip" = "#5E3C99", 
                          "Center" = "#FDB863", 
                          "Arm" = "#4393C3")



library(pals)
# Cluster color palettes 
group_col <- kelly()[c(3, 4, 5, 6, 7, 8, 10, 9)]
names(group_col) <- c("1", "2", "3", "4", "5", "6", "7", "outlier")
group_col8 <- kelly()[c(3, 4, 5, 6, 7, 8, 10, 11, 9)]
names(group_col8) <- c("1", "2", "3", "4", "5", "6", "7", "8", "outlier")
group_col9 <- kelly()[c(3, 4, 5, 6, 7, 8, 10, 11, 12, 9)]
names(group_col9) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "outlier")


# Hawaii island color palettes
isl_cols <- alphabet()[c(1, 2, 3, 4, 6, 5)]
names(isl_cols) <- c("Kauai", "Oahu", "Molokai", "Maui", "Big Island", "unknown")







# Function to grab ggrepel boxes
findboxes <- function(
    df, xcol, ycol,
    box_padding_x, box_padding_y,
    point_padding_x, point_padding_y,
    xlim, ylim,
    force = 1e-7, maxiter = 20000
) {
  
  # x and y posiitons as a dataframe
  posdf <- df[c(xcol, ycol)]
  
  # returnd a df where columns are points
  boxdf <- apply(posdf, 1, function(row) {
    xval <- row[xcol]
    yval <- row[ycol]
    return(c(
      xval - box_padding_x / 2,
      yval - box_padding_y / 2,
      xval + box_padding_x / 2,
      yval + box_padding_y / 2
    ))
  })
  # columns are x1,y1,x2,y2
  boxmatrix <- as.matrix(t(boxdf))
  
  moved <- ggrepel:::repel_boxes(
    data_points = as.matrix(posdf),
    point_padding_x = point_padding_x,
    point_padding_y = point_padding_y,
    boxes = boxmatrix,
    xlim = xlim,
    ylim = ylim,
    hjust = 0.5,
    vjust = 0.5,
    force = force,
    maxiter = maxiter
  )
  
  finaldf <- cbind(df$isotype, posdf, moved)
  names(finaldf) <- c("isotype", "x1", "y1", "x2", "y2")
  return(finaldf)
}



# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}









