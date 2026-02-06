library(shiny)
library(IsoplotR)
library(glue)

# Function to map colours to groups based on the number of groups present in
# the sample
map_colours <- function(letters_vec) {
  if (!is.character(letters_vec)) {
    stop("Error: input must be a character vector.")
  }
  if (length(letters_vec) == 0) {
    stop("Error: input vector is empty.")
  }
  
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

# Input options
ui <- fluidPage(
  titlePanel("GSWA Plotter"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Upload UPb CSV File", accept = ".csv"),
      hr(),
      radioButtons("plot_type", "Select Plot Type:",
                   choices = c("Concordia" = "concordia", 
                               "KDE" = "kde"),
                   selected = "concordia"),
      numericInput("ngrains","Number of grains",value=0),
      numericInput("Threshold", "Discordance threshold (%)\nCurrent value plots all grains", value = 1000),
      numericInput("n_modes", "Number of labelled modes for KDE", value = 5),
      numericInput("pdf_w", "PDF Width", value = 14),
      numericInput("pdf_h", "PDF Height", value = 9.5),
      downloadButton("download_pdf", "Download PDF")
    ),
    mainPanel(
      plotOutput("concordiaPlot", height = "700px")
    )
  )
)


server <- function(input, output) {
  
  # Reactive logic
  processed_bundle <- reactive({
    req(input$file1)
    
    # Read the file
    df <- read.csv(input$file1$datapath)
    
    # Sort by Group
    df_sorted <- df[order(df$Group), ] 
    
    # Generate palette (Returns unique colors)
    palette_vec <- map_colours(df_sorted$Group)
    
    # Process group column
    group_idx <- grep("(?i)^group$", names(df_sorted))
    group_clean <- toupper(trimws(as.character(df_sorted[, group_idx])))
    
    # Value map and GroupNum
    # unique() on sorted data preserves the alphabetical encounter order
    u_groups <- unique(group_clean)
    value_map <- setNames(seq_along(u_groups), u_groups)
    GroupNum_vec <- unname(value_map[group_clean])
    
    # Final matrix prep
    df_for_mat <- df_sorted
    df_for_mat$Group <- NULL
    cols <- c("U238Pb206","errU238Pb206","Pb207Pb206","errPb207Pb206")
    mat <- unname(as.matrix(df_for_mat[, cols, drop = FALSE]))
    storage.mode(mat) <- "double"
    data <- IsoplotR::as.UPb(mat, format = 2, ierr = 1)

    # Return as a LIST
    list(
      data = data,
      palette = palette_vec,
      GroupNum = GroupNum_vec,
      unique_groups = u_groups, # Added for debugging
      xmin = min(mat[,1]) * 0.9,
      xmax = max(mat[,1]) * 1.1,
      ymin = min(mat[,3]) * 0.9,
      ymax = max(mat[,3]) * 1.1,
      SampleName = strsplit(input$file1$name, "_")[[1]][1],
      n_analyses = length(GroupNum_vec)
    )
  })
  
  # 2. Debug: Print the GroupNum Vector and Palette Strings
#  output$debugPalette <- renderPrint({
    # Note the () after processed_bundle!
#    p <- processed_bundle() 
    
#    cat("--- UNIQUE GROUPS FOUND ---\n")
#    print(p$unique_groups)
    
#    cat("\n--- PALETTE (HEX CODES) ---\n")
#    print(p$palette)
    
#    cat("\n--- FULL GROUPNUM VECTOR ---\n")
#    print(p$GroupNum)
#  })


  # The Plotting Logic
  # Function to plot concordia
  draw_concordia <- function() {
    
    # Processed data
    p <- processed_bundle()
    
    # Call concordia
    par(tck = 0.02)
    concordia(
      p$data,
      xlim = c(p$xmin, p$xmax),
      ylim = c(p$ymin, p$ymax),
      type = 2,
      levels = p$GroupNum,
      ellipse.fill = p$palette,
      ellipse.stroke = "black",
      oerr = 1,
      exterr = TRUE,
      ticks = 10,
      concordia.col = "#cccccc",
      show.numbers = FALSE,
      clabel = NA,
      bty = "n",
      ann = FALSE
    )
    
    # Labels for common-Pb corrected ratios
    title(xlab = expression(paste(""^238, "U/"^206, "Pb*")), line = 2.5, cex.lab = 1.5)
    title(ylab = expression(paste(""^207, "Pb*/"^206, "Pb*")), line = 2.5, cex.lab = 1.5)
    
    usr <- par("usr")
    x_pos <- usr[2] - 0.25 * diff(usr[1:2])
    y_pos <- usr[4] - 0.25 * diff(usr[3:4])
    
    # Sample name and number of analyses
    text(x_pos, y_pos, glue("{p$SampleName}"), font = 2, cex = 1.5, adj = 0.5)
    text(x_pos,
         y_pos - 0.035 * diff(usr[3:4]),
         glue("{p$n_analyses} analyses of {input$ngrains} grains"),
         font = 2, cex = 1.0, adj = 0.5)
    
    #mtext("data-point error ellipses are 1s", 
    #      side = 3, adj = 1, line = -1.0, cex = 0.8)
  }
  
  # Function to plot KDE
  draw_kde <- function() {
    p <- processed_bundle()
    
    
    modes <- round(kde(
      data,
      cutoff.disc = discfilter(
        option="r",before=TRUE,cutoff=c(-input$Threshold,input$Threshold)),
      nmodes = input$n_modes,
      plot=FALSE
    )$modes[,1],0)
     
    par(tck = 0.02)
    kde(
      p$data,
      bty="o",
      kde.col= NA,
      hist.col = "#60cdf740",
      cutoff.disc = discfilter(
        option="r",before=TRUE,cutoff=c(-input$Threshold,input$Threshold)),
      xlab = "Age (Ma)",
      ylab = "n Analyses",
      binwidth = 50,
      ann=FALSE,
      sigdig=3
    )
    
    modes <- sort(modes)
    
    # Get plot boundaries: [x_min, x_max, y_min, y_max]
    coords <- par("usr")
    x_left <- coords[2] - 0.10 * diff(coords[1:2]) 
    y_top  <- coords[4] - 0.05 * diff(coords[3:4])
    line_height <- 0.04 * diff(coords[3:4])
    col_width   <- 0.10 * diff(coords[1:2]) # Horizontal distance between columns
    
    # 2. Determine if we have multiple columns
    num_cols <- ceiling(length(modes) / 8)
    # Calculate the center point based on how many columns exist
    # If 1 col, center is x_left. If 2 cols, center is x_left + (col_width/2)
    header_x_center <- x_left + ((num_cols - 1) * col_width) / 2
    
    # 3. Draw Centered Headers
    # adj = 0.5 centers the text on the calculated X coordinate
    text(header_x_center, y_top, labels = p$SampleName, font = 2, adj = 0.5, cex = 1.2)
    text(header_x_center, y_top - line_height, labels = "maxima", font = 1, adj = 0.5, cex = 1.0)
    
    # 4. Draw the Peak Ages (Left-aligned under the header area)
    for (i in seq_along(modes)) {
      col_idx <- (i - 1) %/% 8  # 0 for first 8, 1 for next 8, etc.
      row_idx <- (i - 1) %% 8   # 0 through 7
      
      x_pos <- x_left + (col_idx * col_width)
      y_pos <- y_top - (line_height * (row_idx + 2)) # +2 to start below "maxima"
      
      text(x_pos, y_pos, labels = modes[i], adj = 0, cex = 0.9)
    }
    
    
  }
  
  # Draw based on selected plot type
  output$concordiaPlot <- renderPlot({
    if (input$plot_type == "concordia"){
    draw_concordia()}
    else {
      draw_kde()
    }
  })
  
  
  
  output$download_pdf <- downloadHandler(
    filename = function() { 
      # Dynamically name the file based on the plot type
      paste0(processed_bundle()$SampleName, "_", input$plot_type, ".pdf") 
    },
    content = function(file) {
      # 1. Start the PDF device
      pdf(file, width = input$pdf_w, height = input$pdf_h)
      
      # 2. Check the toggle inside the content function
      if (input$plot_type == "concordia") {
        draw_concordia()
      } else {
        draw_kde()
      }
      
      # 3. Close the device
      dev.off()
    }
  )
}

shinyApp(ui, server)

