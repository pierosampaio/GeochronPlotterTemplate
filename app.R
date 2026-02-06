library(shiny)
library(IsoplotR)
library(glue)
library(DT)


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
      radioButtons("Pbcorr","Is the data common-Pb corrected?",
                   choices = c("Yes"="Yes","No"="No"),selected="Yes"),
      radioButtons("CalcAge", "Calculate age?",
                   choices = c("No"=0,
                               "Concordia age"=1,
                               "Discordia"=2),
                   selected = 0),
      numericInput("ngrains","Number of grains",value=0),
      numericInput("Threshold", "Discordance threshold (%)\n(Current value plots all grains)", value = 1000),
      numericInput("n_modes", "Number of labelled modes for KDE", value = 5),
      numericInput("w_bins", "Width of bins for histogram on KDE", value = 50),
      numericInput("pdf_w", "PDF Width", value = 14),
      numericInput("pdf_h", "PDF Height", value = 9.5),
      downloadButton("download_csv", "Download Edited CSV"),
      downloadButton("download_pdf", "Download PDF")
    ),
    mainPanel(
      plotOutput("concordiaPlot", height = "700px"),
      hr(),
      h3("Edit groups below:"),
      DT::DTOutput("editable_table")
    )
  )
)


server <- function(input, output) {
  
  # 1. Initialize reactive storage for the raw data matrix
  v <- reactiveValues(data_matrix = NULL, df_metadata = NULL)
  
  # 2. Upload Logic: Store the raw dataframe and matrix
  observeEvent(input$file1, {
    df <- read.csv(input$file1$datapath)
    
    # Ensure a Group column exists
    if (!"Group" %in% names(df)) df$Group <- "I" # If there is no group assign I
    
    v$df_metadata <- df[order(df$Group), ]
    v$data_matrix <- as.matrix(v$df_metadata[, c("U238Pb206",
                                                 "errU238Pb206",
                                                 "Pb207Pb206",
                                                 "errPb207Pb206")])
  })
  
  # 3. Render the Editable Table
  output$editable_table <- DT::renderDT({
    req(v$df_metadata)
    DT::datatable(
      v$df_metadata,
      editable = 'cell',
      selection = "multiple",
      options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # 4. CAPTURE EDITS: This is the trigger for reactivity
  observeEvent(input$editable_table_cell_edit, {
    info <- input$editable_table_cell_edit
    # Update the local dataframe
    v$df_metadata[info$row, info$col] <- info$value
    # If the "Group" column was edited, the plot will naturally redraw 
    # because 'processed_bundle' depends on 'v$df_metadata'
  })
  

  
  # Reactive logic
  processed_bundle <- reactive({
    req(input$file1)
    req(v$df_metadata)
    df_sorted <- v$df_metadata
    # Read the file
    #df <- read.csv(input$file1$datapath)
    
    # Sort by Group
    #df_sorted <- df[order(df$Group), ] 
    
    palette_vec <- map_colours(df_sorted$Group)
    
    # Process group column to numeric levels for IsoplotR
    group_clean <- toupper(trimws(as.character(df_sorted$Group)))
    u_groups <- unique(group_clean)
    value_map <- setNames(seq_along(u_groups), u_groups)
    GroupNum_vec <- unname(value_map[group_clean])
    
    mat <- unname(as.matrix(df_sorted[, c("U238Pb206","errU238Pb206","Pb207Pb206","errPb207Pb206")]))
    storage.mode(mat) <- "double"
    data <- IsoplotR::as.UPb(mat, format = 2, ierr = 1)
    
    list(
      data = data,
      palette = palette_vec,
      GroupNum = GroupNum_vec,
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
    selected_rows <- input$editable_table_rows_selected
    plot_levels <- p$GroupNum
    
    highlight_col <- "#FF00FFFF"
    current_palette <- p$palette
    
    if (length(selected_rows) > 0) {
      highlight_idx <- max(p$GroupNum) + 1
      plot_levels[selected_rows] <- highlight_idx
      current_palette <- c(current_palette, highlight_col)
    }
    
    # Call concordia
    #par(tck = 0.02)
    concordia(
      p$data,
      xlim = c(p$xmin, p$xmax),
      ylim = c(p$ymin, p$ymax),
      cutoff.disc = discfilter(
        option="r",before=TRUE,cutoff = c(-input$Threshold,input$Threshold)
      ),
      type = 2,
      levels = plot_levels,
      ellipse.fill = current_palette,
      ellipse.stroke = "black",
      oerr = 1,
      exterr = TRUE,
      ticks = 10,
      concordia.col = "#cccccc40",
      show.numbers = FALSE,
      clabel = NA,
      show.age = as.numeric(input$CalcAge),
      bty = "n",
      ann = FALSE
    )
    
    # Labels for common-Pb corrected ratios
    if (input$Pbcorr == "Yes") {
    title(xlab = expression(paste(""^238, "U/"^206, "Pb*")), line = 2.5, cex.lab = 1.5,
          family = "sans")
    title(ylab = expression(paste(""^207, "Pb*/"^206, "Pb*")), line = 2.5, cex.lab = 1.5,
          family = "sans")
    } 
    else {
    title(xlab = expression(paste(""^238, "U/"^206, "Pb")), line = 2.5, cex.lab = 1.5,
          family = "sans")
    title(ylab = expression(paste(""^207, "Pb/"^206, "Pb")), line = 2.5, cex.lab = 1.5,
          family = "sans")  
    }
    
    usr <- par("usr")
    x_pos <- usr[2] - 0.15 * diff(usr[1:2])
    y_pos <- usr[4] - 0.15 * diff(usr[3:4])
    
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
      p$data,
      cutoff.disc = discfilter(
        option="r",before=TRUE,cutoff=c(-input$Threshold,input$Threshold)),
      nmodes = input$n_modes,
      plot=FALSE,
      show.n = FALSE
    )$modes[,1],0)
     
    par(tck = 0.02)
    old_col <- par("col.main")
    par(col.main = "#FFFFFF00")
    kde(
      p$data,
      bty="o",
      kde.col= NA,
      hist.col = "#60cdf740",
      cutoff.disc = discfilter(
        option="r",before=TRUE,cutoff=c(-input$Threshold,input$Threshold)),
      xlab = "Age (Ma)",
      ylab = "Number of analyses",
      binwidth = input$w_bins,
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
      par(family = "sans")
      draw_concordia()}
    else {
      par(family="sans")
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
  
  output$download_csv <- downloadHandler(
    filename = function() {
      # Dynamically name the file based on the original sample name
      paste0(processed_bundle()$SampleName, "_edited.csv")
    },
    content = function(file) {
      # Ensure there is data to download
      req(v$df_metadata)
      
      # Write the reactive dataframe to the download file path
      # row.names = FALSE prevents an extra index column from being added
      write.csv(v$df_metadata, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)

