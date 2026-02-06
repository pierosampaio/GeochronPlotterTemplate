# set wd

# Import required libraries 
library(IsoplotR)
library(dplyr)
library(glue)



# define a function that generates a palette for a given number of letter codes
# based on predefined colours

map_colours <- function(letters_vec) {
  # Validation
  if (!is.character(letters_vec)) {
    stop("Error: input must be a character vector.")
  }
  if (length(letters_vec) == 0) {
    stop("Error: input vector is empty.")
  }
  
  # mapping
  colour_map <- c(
    I = "#f6eb1480", # I
    M = "#6ebe4480", # M
    Y = "#21aae280", # Y
    S = "#ee2a2980", # S
    X = "#231f2580", # X
    P = "#9467bd80", # P
    D = "#FFFFFF00"  # D
  )
  
  colours_vec <- colour_map[letters_vec]
  
  return(unique(unname(colours_vec)))
}

# Read df
file = "Siliciclastic_test.csv"
df <- read.csv(file)
df_sorted <- df[order(df$Group), ] # Sort by letter code
palette <- map_colours(df_sorted$Group)



# Extract and process group column
group <- df_sorted[, grep("(?i)^group$",names(df_sorted))]
group_clean <- toupper(trimws(as.character(group)))
value_map <- setNames(seq_along(unique(group_clean)), unique(group_clean))
GroupNum <- unname(value_map[group_clean])

# Remove group column before passing it to isoplotR

df$Group <- NULL

cols <- c("U238Pb206","errU238Pb206","Pb207Pb206","errPb207Pb206")
mat <- unname(as.matrix(df[, cols, drop = FALSE]))
storage.mode(mat) <- "double"

SampleName <- strsplit(file,"_")[[1]][1]
n_analyses <- length(GroupNum)
n_grains <- 65

# Pass to isoplotR

data <- IsoplotR::as.UPb(mat, format = 2, ierr = 1)

## Prepare export pdf
#cairo_pdf(glue("{SampleName}_TW.pdf"), width = 14, height = 9.5,
#          family = "Arial")

par(tck = 0.02)

xmin <- min(mat[,1])*0.9
xmax <- max(mat[,1])*1.1
ymin <- min(mat[,3])*0.9
ymax <- max(mat[,3])*1.1


## Suppress default axes
concordia(
  data,
  xlim = c(xmin, xmax),
  ylim = c(ymin, ymax),
  type = 2,
  levels = GroupNum,
  ellipse.fill = palette,
  ellipse.stroke = "black",
  oerr = 1,
  exterr = TRUE,
  ticks = 15,
  concordia.col = "salmon",
  show.numbers = FALSE
)

# Title and subtitle
usr <- par("usr")
x_pos <- usr[2] - 0.25 * diff(usr[1:2])
y_pos <- usr[4] - 0.25 * diff(usr[3:4])

text(x_pos, y_pos, glue("{SampleName}"), font = 2, cex = 1.5, adj = 0.5)
text(x_pos,
     y_pos - 0.03 * diff(usr[3:4]),
     glue("{n_analyses} analyses of {n_grains} zircons"),
     font = 2,
     cex = 1.0,
     adj = 0.5)

mtext("data-point error ellipses are 1s",
      side = 3,
      adj = 1,
      line = -1.0,
      cex = 0.8)
dev.off()

## export png
png_name = glue("{SampleName}_TW.png")
png(png_name,width=14,height=9.5,units="in",res=600)

par(tck = 0.02)

concordia(
  data,
  xlim = c(xmin, xmax),
  ylim = c(ymin, ymax),
  type = 2,
  levels = GroupNum,
  ellipse.fill = palette,
  ellipse.stroke = "black",
  oerr = 1,
  exterr = TRUE,
  ticks = 15,
  concordia.col = "#cccccc",
  show.numbers = FALSE
)


# Title and subtitle
usr <- par("usr")
x_pos <- usr[2] - 0.25 * diff(usr[1:2])
y_pos <- usr[4] - 0.25 * diff(usr[3:4])

text(x_pos, y_pos, glue("{SampleName}"), font = 2, cex = 1.5, adj = 0.5)
text(x_pos,
     y_pos - 0.03 * diff(usr[3:4]),
     glue("{n_analyses} analyses of {n_grains} zircons"),
     font = 2,
     cex = 1.0,
     adj = 0.5)

mtext("data-point error ellipses are 1s",
      side = 3,
      adj = 1,
      line = -1.0,
      cex = 0.8)

dev.off()


KDE <- kde(data,
    bty="o",
    nmodes="all",
    xlab = "Age (Ma)",
    normalise=TRUE)
KDE
