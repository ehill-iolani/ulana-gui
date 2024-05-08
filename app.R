library(shiny)
library(ggplot2)
library(dplyr)
library(stringr)
library(plotly)
library(shinydashboard)
library(shinyalert)
library(DT)
library(htmlwidgets)

# UI for the ULANA whole genome assembler
ui <- dashboardPage(skin = "red",
  dashboardHeader(title = "ULANA"),
  dashboardSidebar(sidebarMenu(
      menuItem("Upload Data", tabName = "upload", icon = icon("file")),
      menuItem("QC and Filtering", tabName = "qc", icon = icon("filter")),
      menuItem("Flye Assembly", tabName = "flye", icon = icon("dna")),
      menuItem("Medaka Polished Assembly", tabName = "medaka", icon = icon("dna")),
      menuItem("CheckM Completion", tabName = "checkm", icon = icon("check")),
      menuItem("Prokka Annotation", tabName = "prokka", icon = icon("database")),
      menuItem("Identifier Genes", tabName = "idg", icon = icon("dna"))
  )),

  # Formatting to make it iolani colors
  dashboardBody(tags$style(HTML("
  .box.box-solid.box-primary>.box-header {
    color:#fff;
    background:#000000
  }

  .box.box-solid.box-primary {
    border-bottom-color:#000000;
    border-left-color:#000000;
    border-right-color:#000000;
    border-top-color:#000000;
  }

  .box.box-primary>.box-header {
    color:#000000;
    background:#fff
  }

  .box.box-primary {
    border-bottom-color:#000000;
    border-left-color:#000000;
    border-right-color:#000000;
    border-top-color:#000000;
  }

  .skin-red .main-sidebar {
    background-color: #000000;
  }")),
    tabItems(
      tabItem(tabName = "upload",
            fluidRow(
            box(
              title = "Fastq Upload", status = "primary", solidHeader = TRUE, width = 12,
              fileInput("fastqs", "Input fastq files", accept = c("fastq", "fastq.gz"), multiple = TRUE),
              actionButton("fastqhelp", "Help!")),
            box(
              title = "Ulana Parameters", status = "primary", solidHeader = TRUE, width = 12,
              numericInput("lowerlength", "Minimum read length", value = 1000),
              numericInput("quality", "Minimum Quality Score", value = 9),
              numericInput("threads", "Threads", value = 4),
              checkboxInput("medaka", "Polish with Medaka", value = FALSE),
                conditionalPanel(
                  condition = "input.medaka == true",
                  selectInput("medakaModel", "Select Medaka Model", choices = c("r1041_e82_400bps_fast_g632", "r1041_e82_400bps_hac_v4.1.0", "r1041_e82_400bps_sup_v4.1.0"))
                ),
              checkboxInput("prokka", "Annotate with Prokka", value = FALSE),
              checkboxInput("checkm", "Check completion with CheckM", value = FALSE),
              checkboxInput("idg", "Extract Identifier genes", value = FALSE)),
            box(
              title = "Run Workflow", status = "primary", solidHeader = TRUE, width = 12,
              actionButton("run", "Run ULANA!"),
              actionButton("ulanahelp", "Help!")))),
      tabItem(tabName = "qc",
        fluidRow(
          box(
            title = "QC and Filtering",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotlyOutput("qc_plot"),
            downloadButton("qc_plot_download", "Download QC Report"),
            downloadButton("filtered_fastq", "Download Filtered Fastq")),
          box(
            title = "QC Table",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            DTOutput("qc"),
            downloadButton("qc_download", "Download QC Table")))),
      tabItem(tabName = "flye",
        fluidRow(
          box(
            title = "Assembly Stats",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DTOutput("ass_info"),
            downloadButton("ass_info_download", "Download Assembly Stats")),
          box(
            title = "Assembly Visualization",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            imageOutput("gfa"),
            downloadButton("ass_plot_download", "Download Assembly Visualization"),
            downloadButton("ass_download", "Download Flye Assembly")))),
      tabItem(tabName = "medaka",
        fluidRow(
          box(
            title = "Medaka Polishing",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            downloadButton("ass_polished_download", "Download Medaka Polished Assembly")))),
      tabItem(tabName = "checkm",
        fluidRow(
          box(
            title = "Genome Completeness",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            div(style = 'overflow-x: auto;', DTOutput("ass_complete")),
            downloadButton("ass_complete_download","Download CheckM Completion Analysis")))),
      tabItem(tabName = "prokka",
        fluidRow(
          box(
            title = "Genome Annotation",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DTOutput("prokka_annotation"),
            downloadButton("prokka_annotation_download", "Download Prokka Annotation")))),
      tabItem(tabName = "idg",
        fluidRow(
          box(
            title = "Identifier Genes",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            downloadButton("idg_download", "Download Identifier Genes"))))
)))

server <- shinyServer(function(input, output, session) {
  # Create a reactive value to track the session status
  session_status <- reactiveValues(new_session = TRUE)

  # Create an observer that checks if it's a new session at the start of the server function
  observe({
    if (isolate(session_status$new_session)) {
      system("rm -rf /home/processing/*")
      isolate(session_status$new_session <- FALSE)
    }
  }) 
  
  # Set large upload size limit (server side)
  options(shiny.maxRequestSize = 250 * 1024^2)

  observeEvent(input$run, {
    showModal(modalDialog("Running the ULANA pipeline, preprocessing the FASTQs", footer = NULL))

    # Make a directory to store input files
    dir.create(file.path("/home/processing/fastq"))

    # Write to fastq directory
    files <- input$fastqs$name
    for (i in 1:length(files)) {
      file.copy(input$fastqs$datapath[i], file.path("/home/processing/fastq", files[i]))
    }

    # Check if the files are fastq or fastq.gz; unzip if necessary
    setwd("/home/processing")
    if (any(str_detect(files, ".gz"))) {
      system("gunzip fastq/*.gz")
    }

    # Cat the files together and rezip
    system("cat fastq/*.fastq > wgs.fastq")
    system("gzip wgs.fastq")
    system("rm -r fastq wgs.fastq")
    removeModal()

    ######################################
    ### QC and filtering using chopper ###
    ######################################
    showModal(modalDialog("Running the ULANA pipeline, filtering FASTQs with chopper", footer = NULL))

    # Run chopper with user input
    chopper <- paste("gunzip -c wgs.fastq.gz | chopper -q", input$quality, "--minlength", input$lowerlength, "-t", input$threads, "| gzip > wgs_qc.fastq.gz")
    system(chopper, intern = TRUE)

    # Summarize QC
    raw <- readLines("wgs.fastq.gz")
    raw <- data.frame(matrix(unlist(strsplit(raw, "\t")), ncol = 1, byrow = TRUE))
    seq <- raw[seq(2, nrow(raw), 4), ]
    seq_len <- nchar(seq)
    seq_len_dt <- length(seq_len)

    filt <- readLines("wgs_qc.fastq.gz")
    filt <- data.frame(matrix(unlist(strsplit(filt, "\t")), ncol = 1, byrow = TRUE))
    seq2 <- filt[seq(2, nrow(filt), 4), ]
    seq_len2 <- nchar(seq2)
    seq_len2_dt <- length(seq_len2)

    # Create a data frame for the QC
    qc <- data.frame(seq_len_dt, seq_len2_dt)
    colnames(qc) <- c("Raw", "Filtered")
    rownames(qc) <- c("Read Number")

    # Display the QC table
    output$qc <- renderDataTable(qc)

    # Prepare QC table for download
    output$qc_download <- downloadHandler(
      filename = "qc_table.csv",
      content = function(file) {
        write.table(qc, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )

    # Plot the results of the QC
    qc_plot <- ggplot() +
      geom_histogram(aes(seq_len, fill = "Raw"), binwidth = 1000) +
      geom_histogram(aes(seq_len2, fill = "Filtered"), binwidth = 1000) +
      ggtitle("Raw vs Filtered Sequence Lengths") +
      xlab("Sequence Length") +
      ylab("Number of Reads") +
      scale_fill_manual(values = c("Raw" = "grey", "Filtered" = "green")) +
      labs(fill = "Read Catagory") +
      theme_minimal()
    output$qc_plot <- renderPlotly(ggplotly(qc_plot))

    # Prepare QC plot for download
    output$qc_plot_download <- downloadHandler(
      filename = "qc_plot.html",
      content = function(file) {
        saveWidget(ggplotly(qc_plot), file)
      }
    )

    # Prepare filtered fastq for download
    output$filtered_fastq <- downloadHandler(
      filename = "wgs_qc.fastq.gz",
      content = function(file) {
        file.copy("wgs_qc.fastq.gz", file)
      },
      contentType = "application/gzip"
    )

    removeModal()

    ###########################
    ### Assembly with  Flye ###
    ###########################
    showModal(modalDialog("Running the ULANA pipeline, assembling genome with Flye", footer = NULL))

    # Run Flye
    flye <- paste("flye --nano-hq wgs_qc.fastq.gz --out-dir flye_assembly --threads", input$threads)
    system(flye)

    # Display the assembly stats
    ass_info <- read.csv("./flye_assembly/assembly_info.txt", sep = "\t")
    ass_info <- as.data.frame(ass_info)
    output$ass_info <- renderDataTable(ass_info)

    # Prepare assembly states for download
    output$ass_info_download <- downloadHandler(
      filename = "assembly_stats.csv",
      content = function(file) {
        write.table(ass_info, file, sep = "\t", row.names = FALSE)
      }
    )

    # Display the assembly gfa
    system("Bandage image ./flye_assembly/assembly_graph.gfa ./flye_assembly/assembly_graph.png --lengths --depth --fontsize 6")

    output$gfa <- renderImage({
      list(src = "./flye_assembly/assembly_graph.png",
        contentType = "image/png",
        height = "100%",
        alt = "Run the pipeline for visualization")
    }, deleteFile = FALSE)

    # Prepare asssembly gfa for download
    output$ass_plot_download <- downloadHandler(
      filename = "assembly_graph.png",
      content = function(file) {
        file.copy("./flye_assembly/assembly_graph.png", file)
      }
    )

    removeModal()

    ######################################
    ### Polish with Medaka if selected ###
    ######################################
    observeEvent(input$medaka, {
      if (input$medaka) {
        showModal(modalDialog("Running the ULANA pipeline, polishing the assembly with Medaka", footer = NULL))

        medaka_consensus <- paste("medaka_consensus -i wgs_qc.fastq.gz -d ./flye_assembly/assembly.fasta -o medaka_polished_assembly -t 2 -m", input$medakaModel)
        system(medaka_consensus)
        medaka_stitch <- paste("medaka stitch ./medaka_polished_assembly/*.hdf ./flye_assembly/assembly.fasta ./medaka_polished_assembly/polished_consensus.fasta")
        system(medaka_stitch)

        # Prepare polished assembly for download
        output$ass_polished_download <- downloadHandler(
          filename = "polished_assembly.fasta",
          content = function(file) {
            file.copy("./medaka_polished_assembly/polished_consensus.fasta", file)
          }
        )

        removeModal()
      }
    })

    ########################################
    ### Annotate with prokka if selected ###
    ########################################
    observeEvent(input$prokka, {
      if (input$prokka) {
        showModal(modalDialog("Running the ULANA pipeline, annotating the assembly with PROKKA", footer = NULL))
        
        # Run prokka based on if medaka was selected
        if (input$medaka) {
          prokka_command <- paste("/tools/prokka/bin/prokka --outdir prokka_out ./medaka_polished_assembly/polished_consensus.fasta --prefix prokka_annotation --cpus", input$threads)
        } else {
          prokka_command <- paste("/tools/prokka/bin/prokka --outdir prokka_out ./flye_assembly/assembly.fasta --prefix prokka_annotation --cpus", input$threads)
        }
        system(prokka_command)

        # Display the prokka annotation
        prokka_annotation <- read.csv("./prokka_out/prokka_annotation.tsv", sep = "\t")
        prokka_annotation <- as.data.frame(prokka_annotation)
        output$prokka_annotation <- renderDataTable(prokka_annotation)

        # Prepare prokka annotation for download
        output$prokka_annotation_download <- downloadHandler(
          filename = "prokka_annotation.tsv",
          content = function(file) {
            write.table(prokka_annotation, file, sep = "\t", row.names = FALSE, quote = FALSE)
          }
        )

        removeModal()
      }
    })

    ####################################
    ### Check completion with CheckM ###
    ####################################
    observeEvent(input$checkm, {
      if (input$checkm) {
        showModal(modalDialog("Running the ULANA pipeline, checking completion with CheckM", footer = NULL))

        if (input$medaka) {
          checkm_command <- paste("checkm lineage_wf --reduced_tree -t", input$threads, "-x fasta ./medaka_polished_assembly checkm")
        } else {
          checkm_command <- paste("checkm lineage_wf --reduced_tree -t", input$threads, "-x fasta ./flye_assembly checkm")
        }
        system(checkm_command)
        system("mkdir checkm/summary")
        checkm_summary <- paste("checkm qa -o 2 --tab_table -t", input$threads, "-f checkm/summary/summary.tsv checkm/lineage.ms checkm")
        system(checkm_summary)

        # Display the checkm summary
        checkm_summary <- read.csv("checkm/summary/summary.tsv", sep = "\t")
        checkm_summary <- as.data.frame(checkm_summary)
        output$ass_complete <- renderDataTable(checkm_summary, options = list(scrollX = TRUE))

        # Prepare checkm summary for download
        output$ass_complete_download <- downloadHandler(
          filename = "checkm_summary.tsv",
          content = function(file) {
            write.table(checkm_summary, file, sep = "\t", row.names = FALSE)
          }
        )

        removeModal()
      }
    })

    ################################
    ### Extract Identifier genes ###
    ################################
    observeEvent(input$idg, {
      if (input$idg & input$prokka) {
        showModal(modalDialog("Running the ULANA pipeline, extracting identifier genes", footer = NULL))

        # Prepare to extract identifier genes
        system("mkdir identifier_genes")
        rrna <- paste("echo '>16S ribosomal RNA' > identifier_genes/16S_rRNA.fasta && awk '/16S ribosomal RNA/{flag=1;next}/^>/{flag=0}flag' prokka_out/prokka_annotation.ffn >> identifier_genes/16S_rRNA.fasta")
        dnaa <- paste("echo '>Chromosomal replication initiator protein DnaA' > identifier_genes/DnaA.fasta && awk '/Chromosomal replication initiator protein DnaA/{flag=1;next}/^>/{flag=0}flag' prokka_out/prokka_annotation.ffn >> identifier_genes/DnaA.fasta")
        rpob <- paste("echo '>DNA-directed RNA polymerase subunit beta' > identifier_genes/RpoB.fasta && awk '/DNA-directed RNA polymerase subunit beta/{flag=1;next}/^>/{flag=0}flag' prokka_out/prokka_annotation.ffn >> identifier_genes/RpoB.fasta")

        # Extract identifier genes
        system(rrna)
        system(dnaa)
        system(rpob)

        # Prepare identifier genes for download
        output$idg_download <- downloadHandler(
          filename = "identifier_genes.zip",
          content = function(file) {
            system("zip -r identifier_genes.zip identifier_genes")
            file.copy("identifier_genes.zip", file)
          },
          contentType = "application/zip"
        )

        removeModal()
      }
    })

  })

  # Set new_session to TRUE when the session ends
  session$onSessionEnded(function() {
    session_status$new_session <- TRUE
  })
})

shinyApp(ui = ui, server = server)
