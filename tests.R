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


min(age(data, method="U238-Pb206",
        discordance = discfilter("r",before=TRUE,c(-5,5)))[,"t.68"])
