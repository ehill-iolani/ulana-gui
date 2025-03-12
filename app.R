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
      menuItem("Identifier Genes", tabName = "idg", icon = icon("dna")),
      menuItem("Antimicrobial Resistance Genes", tabName = "amr", icon = icon("dna")),
      menuItem("Download Results Package", tabName = "dl_all", icon = icon("download"))
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
  }
  ")),

    tabItems(
      tabItem(tabName = "upload",
        fluidRow(
          box(
            title = "Fastq Upload", status = "primary", solidHeader = TRUE, width = 12,
            fileInput("fastqs", "Input fastq files", accept = c("fastq", "fastq.gz"), multiple = TRUE),
            actionButton("fastqhelp", "Help!"))),
        fluidRow(
          box(
            title = "Ulana Assembly Parameters", status = "primary", solidHeader = TRUE, width = 4,
            numericInput("lowerlength", "Minimum read length", value = 1000),
            numericInput("quality", "Minimum Quality Score", value = 9),
            numericInput("threads", "Threads", value = 4),
            actionButton("ulanahelp", "Help!")),
          box(
            title = "Post Assembly Analyses (optional)", status = "primary", solidHeader = TRUE, width = 4,
            checkboxInput("medaka", "Polish Flye assembly with Medaka", value = FALSE),
              conditionalPanel(
                condition = "input.medaka == true",
                selectInput("medakaModel", "Select Medaka Model", choices = c("r1041_e82_400bps_fast_g632", "r1041_e82_400bps_hac_v4.3.0", "r1041_e82_400bps_sup_v4.3.0"))),
            checkboxInput("prokka", "Annotate with Prokka", value = FALSE),
            checkboxInput("checkm", "Check completion with CheckM", value = FALSE),
            checkboxInput("idg", "Extract Identifier genes", value = FALSE),
            checkboxInput("amr", "Extract antimicrobial resistance genes with AMRFinder", value = FALSE),
            actionButton("postasshelp", "Help!")),
          box(
            title = "Run Workflow", status = "primary", solidHeader = TRUE, width = 4,
            textInput("sample_name", "Sample Name", value = "sample1"),
            actionButton("run", "Run ULANA!")))),
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
            downloadButton("idg_download", "Download Identifier Genes")))),
      tabItem(tabName = "amr",
        fluidRow(
          box(
            title = "AMR genes",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DTOutput("amr_table"),
            downloadButton("amr_download", "Download AMR protein sequences")))),
      tabItem(tabName = "dl_all",
        fluidRow(
          box(
            title = "Download Results Package",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            downloadButton("dl_all", "Download All Results"))))
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
    system("rm -r fastq")
    removeModal()

    ######################################
    ### QC and filtering using chopper ###
    ######################################
    showModal(modalDialog("Running the ULANA pipeline, filtering FASTQs with chopper", footer = NULL))

    # Run chopper with user input
    chopper <- paste("chopper -q", input$quality, "--minlength", input$lowerlength, "--threads", input$threads, "-i wgs.fastq | pigz --fast -p", input$threads, "> wgs_qc.fastq.gz")
    system(chopper, intern = TRUE)

    # Summarize QC
    raw <- readLines("wgs.fastq")
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

    # Prepare assembly for download
    output$ass_download <- downloadHandler(
      filename = "flye_assembly.fasta",
      content = function(file) {
        file.copy("./flye_assembly/assembly.fasta", file)
      }
    )

    removeModal()

    ######################################
    ### Polish with Medaka if selected ###
    ######################################
    observeEvent(input$medaka, {
      if (input$medaka) {
        showModal(modalDialog("Running the ULANA pipeline, polishing the assembly with Medaka", footer = NULL))

        medaka_consensus <- paste("medaka_consensus -i wgs_qc.fastq.gz -d ./flye_assembly/assembly.fasta -o medaka_polished_assembly -t", input$threads, "-m", input$medakaModel)
        system(medaka_consensus)
        medaka_stitch <- paste("medaka sequence ./medaka_polished_assembly/*.hdf ./flye_assembly/assembly.fasta ./medaka_polished_assembly/polished_consensus.fasta")
        system(medaka_stitch)

        # Prepare only the polished assembly for download
        system("mkdir ./medaka_polished_assembly_nobam")
        system("cp ./medaka_polished_assembly/polished_consensus.fasta ./medaka_polished_assembly_nobam/polished_consensus.fasta")

        # Prepare whole polished assembly results for download
        output$ass_polished_download <- downloadHandler(
          filename = "medaka_polished_assembly.zip",
          content = function(file) {
            system("zip -r medaka_polished_assembly.zip medaka_polished_assembly")
            file.copy("medaka_polished_assembly.zip", file)
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
          filename = "prokka_annotation.zip",
          content = function(file) {
            system("zip -r prokka_annotation.zip prokka_out")
            file.copy("prokka_annotation.zip", file)
          },
          contentType = "application/zip"
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

    ##################################
    ### Find and extract AMR genes ###
    ##################################
    observeEvent(input$amr, {
      if (input$amr & input$prokka) {
        showModal(modalDialog("Running the ULANA pipeline, extracting antimicrobial resistance genes using amrfinder", footer = NULL))

        # Prepare to extract AMR genes
        system("mkdir amr_genes")
        amr <- paste("amrfinder --threads", input$threads, "-a prokka -p ./prokka_out/prokka_annotation.faa -g ./prokka_out/prokka_annotation.gff --print_node --plus --protein_output ./amr_genes/amrfinder_pro.fasta >./amr_genes/amrfinder_pro_results.tsv")

        # Extract AMR genes
        system(amr)

        # Display the AMR genes
        amr_genes <- read.csv("amr_genes/amrfinder_pro_results.tsv", sep = "\t")
        amr_genes <- as.data.frame(amr_genes)
        output$amr_table <- renderDataTable(amr_genes, options = list(scrollX = TRUE))

        # Prepare AMR genes for download
        output$amr_download <- downloadHandler(
          filename = "amr_genes.zip",
          content = function(file) {
            system("zip -r amr_genes.zip amr_genes")
            file.copy("amr_genes.zip", file)
          },
          contentType = "application/zip"
        )

        removeModal()
      }
    })

    ############################
    ### Download all results ###
    ############################
    output$dl_all <- downloadHandler(
      filename = paste(input$sample_name, "_results.zip"),
      content = function(file) {
        # Show a modal dialog indicating that the files are being compressed
        shinyalert::shinyalert("Compressing files, please wait...", type = "info", showConfirmButton = FALSE)

        # Create a temporary directory to store the results
        temp_dir <- tempdir()
        temp_zip <- file.path(temp_dir, "results.zip")

        # Zip all the results if they exist
        system(paste("zip -r", temp_zip, "flye_assembly medaka_polished_assembly_nobam prokka_out checkm identifier_genes amr_genes"))

        # Copy the zip file to the final destination
        file.copy(temp_zip, file)

        # Close the modal dialog
        shinyalert::closeAlert()

        # Alert user that the download is ready
        shinyalert::shinyalert("Download ready!", "Your download is ready. Click the button below to download the results.", type = "success")
      },
      contentType = "application/zip"
    )
  })

  ####################
  ### Help buttons ###
  ####################
  observeEvent(input$fastqhelp, {
    shinyalert(title = "FASTQ Upload", text = "Upload your FASTQ files here. You can upload multiple files at once. Please only process one bacterial genome at a time.", type = "info")
  })

  observeEvent(input$ulanahelp, {
    shinyalert(title = "ULANA Assembly Parameters", text = "Set the minimum read length, quality score, and number of threads for the ULANA pipeline. Increasing the read length typically results in a less fragmented assembly, at the cost of lower the read count. Increaseing the quality score will typically result in a more accurate assembly at the cost of read count. Increasing the number of threads will decrease run time (this should match the number of cores either on your computer or on the number of CPUs requested when starting the app.)", type = "info")
  })

  observeEvent(input$postasshelp, {
    shinyalert(title = "Post Assembly Analyses", text = "Select the post-assembly analyses you would like to run. Medaka - polishes the Flye assembly for a more accurate consensus genome (will increase run time). Prokka - annotates the genome assembly, if the assembly was polished using Medaka the polished assembly will be used otherwise the Flye assembly will be used. CheckM - checks the completeness of the assembly, if the assembly was polished using Medaka the polished assembly will be used otherwise the Flye assembly will be used. Identifier genes - will find and extract 16S, rpoB and dnaA. AMRFinder - will find and extract confirmed antimicrobial genes if they exist in the genome. ", type = "info")
  })

  # Set new_session to TRUE when the session ends
  session$onSessionEnded(function() {
    session_status$new_session <- TRUE
  })
})

shinyApp(ui = ui, server = server)
