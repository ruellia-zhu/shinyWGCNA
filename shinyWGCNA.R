###############################
#	prj: shiny app
#	Assignment: WGCNA by click shiny app
#	Author: Shawn Wang
#	Date: Jan 12, 2021
# V3.0 Updater : Yuntao Zhu
# Update date: 2026/05/08
# Update V3.0: configure WGCNA server threads, support FPKM/TPM labels, restore cleaned genes, and add TOMplot
###############################
select = dplyr::select
### You must set this before app is loaded, or hub gene will not work - by yuntao
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require('devtools')) install.packages('devtools');
if (!require('DESeq2')) BiocManager::install('DESeq2',update = FALSE);
if (!require('shinyjs')) install.packages('shinyjs');
if (!require('dashboardthemes')) install.packages('dashboardthemes');
if (!require('shinydashboard')) install.packages('shinydashboard');
if (!require("DT")) install.packages('DT');
if (!require('shiny')) install.packages('shiny');
if (!require('ggpmisc')) install.packages('ggpmisc');
if (!require('dplyr')) install.packages('dplyr');
if (!require('GO.db')) BiocManager::install('GO.db',update = FALSE);
if (!require('WGCNA')) BiocManager::install('WGCNA',update = FALSE);
if (!require('ComplexHeatmap')) BiocManager::install('ComplexHeatmap',update = FALSE);
if (!require('circlize')) BiocManager::install('circlize',update = FALSE);
if (!require('stringr')) install.packages('stringr');
if (!require('ape')) install.packages('ape');
if (!require('reshape2')) install.packages('reshape2');
if (!require('edgeR')) BiocManager::install('edgeR',update = FALSE);
if (!require('shinythemes')) install.packages('shinythemes');
if (!require('ggplotify')) install.packages('ggplotify');
if (!require('ggprism')) install.packages('ggprism');
if (!require('ggpubr')) install.packages('ggpubr');
if (!require('patchwork')) install.packages('patchwork');
if (!require('tidyverse')) install.packages('tidyverse');
if (!require('shinyjqui')) install.packages('shinyjqui');
if (!require('colourpicker')) install.packages('colourpicker');
suppressMessages(library(devtools))
if (!require('ShinyWGCNA')) devtools::install_github("ShawnWx2019/WGCNAShinyFun",ref = "master");
suppressMessages(library(ShinyWGCNA))
suppressMessages(library(shinyjs))
suppressMessages(library(dashboardthemes))
suppressMessages(library(shinydashboard))
suppressMessages(library(DT))
suppressMessages(library(shiny))
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(library(WGCNA))
suppressMessages(library(stringr))
suppressMessages(library(ape))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(reshape2))
suppressMessages(library(edgeR))
suppressMessages(library(shinythemes))
suppressMessages(library(ggplotify))
suppressMessages(library(ggprism))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(shinyjqui))
suppressMessages(library(ggpubr))
suppressMessages(library(dplyr))
options(shiny.maxRequestSize = 300*1024^2)
options(scipen = 6)
select =  dplyr::select
# type = "unsigned"
# corType = "pearson"
# maxPOutliers = ifelse(corType=="pearson",1,0.05)
# robustY = ifelse(corType=="pearson",T,F)
enableWGCNAThreads(nThreads = 20)
# functions =========================

# 01. UI =========================
## logo
customLogo <- shinyDashboardLogoDIY(

  boldText = "Interactive"
  ,mainText = "WGCNA with GUI"
  ,textSize = 14
  ,badgeText = "v3.0"
  ,badgeTextColor = "white"
  ,badgeTextSize = 2
  ,badgeBackColor = "#40E0D0"
  ,badgeBorderRadius = 3

)

ui <- shinyUI(
  navbarPage(theme = shinytheme("journal"),
             customLogo,
             tabPanel(
               title = "Before Use",
               icon = icon("tree"),
               sidebarLayout(
                 div(id = "Sidebar0",
                     sidebarPanel(
                       width = 3,
                       p("WGCNA Thread Settings", style = "color: #000000;font-size: 24px; font-style:bold"),
                       sliderInput(
                         inputId = "wgcna.threads",
                         label = "Server threads used by WGCNA",
                         min = 4,
                         max = 64,
                         value = 20,
                         step = 1
                       ),
                       helpText("Set how many server CPU threads WGCNA can use for network analysis. Choose a value according to server load and available cores.",
                         style = "color: #7a8788;font-size: 12px; font-style:Italic"),
                       helpText("Default is 20 threads. Lower values reduce server pressure; higher values may speed up large analyses but can affect other users.",
                                style = "color: #7a8788;font-size: 12px; font-style:Italic"),
                       
                     )# sidebarPanel
                 ),# div
                 mainPanel(
                   p("WGCNA threads selected", style = "color: #000000;font-size: 24px; font-style:bold"),
                   textOutput("show.wgcna.threads"),
                   actionButton("set.wgcna.threads",
                                "Apply thread setting"),
                   p("Click the apply button before running WGCNA steps so the selected thread count is used by WGCNA.", style = "color: #000000;font-size: 16px; font-style:bold"),
                   p("WGCNA threads actually configured", style = "color: #000000;font-size: 24px; font-style:bold"),
                   textOutput("current.wgcna.threads"),
                   p("NOTE", style = "color: #000000;font-size: 24px; font-style:bold"),
                   p("1. The old Windows-only R memory-limit adjustment has been removed because it is obsolete in modern R environments.", style = "color: #000000;font-size: 16px; font-style:bold"),
                   p("2. This server-oriented setting controls the number of CPU threads WGCNA may use; it does not change the server memory limit.", style = "color: #000000;font-size: 16px; font-style:bold"),
                   p("3. This app's primary content and original version is developed by Shawn Wang. You can attach Wang's original version in the URL below.", style = "color: #000000;font-size: 16px; font-style:bold"),
                   p("4. This app has never been sold. ", style = "color: #000000;font-size: 16px; font-style:bold"),
                   p("Shawn Wang's original version: https://github.com/ShawnWx2019/WGCNA-shinyApp")
                 )
               ) # sidebarLayout
             ),##tabPanel
             tabPanel(
               useShinyjs(),
               title = "Data import and cleaning",
               icon = icon("file-upload"),
               sidebarLayout(
                 div(id = "Sidebar",
                     sidebarPanel(
                       width = 2,
                       fileInput(
                         inputId = "ExpMat",
                         label = "Upload expression matrix",
                         accept = c(".txt", ".tsv", ".tab", ".csv", ".xls")
                       ),
                       p("accept Tab-delimited .txt/.tsv/.tab/.xls files and comma-separated .csv files",
                         style = "color: #7a8788;font-size: 12px; font-style:Italic"),
                       radioButtons(
                         inputId = "format",
                         label = "Format",
                         choices = c("count", "FPKM/TPM" = "FPKM"),
                         selected = "FPKM"
                       ),
                       p("Normalized expression input includes FPKM, TPM, RPKM, CPM and other already-normalized expression values",
                         style = "color: #7a8788;font-size: 12px; font-style:Italic"),
                       selectInput(
                         inputId = "method1",
                         label = "Normalized method",
                         choices = c(VST = "varianceStabilizingTransformation",
                                     lgCPM = "lgcpm",
                                     FPKM_TPM = "rawFPKM",
                                     logFPKM_TPM = "lgFPKM"),
                         selected = "rawFPKM"
                       ),
                       HTML('<font color = #FF6347  size = 3.2><b>First Time filter</b></font>'),
                       textInput(
                         inputId = "SamPer",
                         label = "Sample percentage",
                         value = "0.9"
                       ),
                       textInput(
                         inputId = "RCcut",
                         label = "Expression Cutoff",
                         value = "1"
                       ),
                       p("Noise removal: for count data, for example, remove features with count less than 10 in more than 90% of samples; for FPKM/TPM, use expression cutoff 1 by default.",
                         style = "color: #7a8788;font-size: 12px; font-style:Italic"),
                       br(),
                       HTML('<font color = #FF6347 size = 3.2><b>Second Time filter</b></font>'),
                       radioButtons(
                         inputId = "CutMethod",
                         label = "Filter Method",
                         choices = c("MAD","Var"),
                         selected = "MAD"
                       ),
                       textInput(
                         inputId = "remain",
                         label = "Reserved genes Num.",
                         value = "8000"
                       ),
                       p("Probesets or genes may be filtered by mean expression or variance (or their robust analogs such as median and median absolute deviation, MAD) since low-expressed or non-varying genes usually represent noise. Whether it is better to filter by mean expression or variance is a matter of debate; both have advantages and disadvantages, but more importantly, they tend to filter out similar sets of genes since mean and variance are usually related.",
                         style = "color: #7a8788;font-size: 12px; font-style:Italic"),
                     )# sidebarPanel
                 ),# div
                 mainPanel(
                   fluidPage(
                     actionButton("toggleSidebar",
                                  "Toggle sidebar"),
                     actionButton("action1", "Update information!"),
                     tabsetPanel(
                       tabPanel(title = "Input file check",height = "500px",width = "100%",
                                icon = icon("check-circle"),
                                htmlOutput("Inputcheck"),
                                htmlOutput("filter1"),
                                htmlOutput("filter2"),

                       ),
                       tabPanel(title = "Preview of Input",height = "500px",width = "100%",
                                icon = icon("table"),
                                DT::dataTableOutput("Inputbl"),
                       ),
                       tabPanel(title = "Cleaned genes",height = "500px",width = "100%",
                                icon = icon("recycle"),
                                p("Search and select genes removed by cleaning, then add them back to the Input gene table.",
                                  style = "color: #7a8788;font-size: 12px; font-style:Italic"),
                                DT::dataTableOutput("washedGenes"),
                                actionButton("addWashedGenes", "Add selected cleaned genes to Input"),
                                htmlOutput("addWashedGenesInfo")
                       ),
                       tabPanel(title = "SampleCluster",height = "500px",width = "100%",
                                icon = icon("tree"),
                                jqui_resizable(
                                  plotOutput("clustPlot")
                                ),
                                downloadButton("downfig1","Download")

                       )# tabPanel
                     )

                   )# fluidPage
                 )#mainPanel
               ) # sidebarLayout
             ),##tabPanel
             tabPanel(
               useShinyjs(),
               title = "SFT and Power Select",
               icon = icon("play-circle"),
               sidebarLayout(
                 div(id = "Sidebar2",
                     sidebarPanel(
                       width = 2,
                       sliderInput(
                         inputId = "CutoffR",
                         label = HTML('R<sup>2</sup> cutoff'),
                         min = 0,
                         max = 1,
                         value = 0.8
                       ),
                       br(),
                       HTML('<font size = 2.5 color = #7a8788><i>WGCNA will generate a recommended power value. If it does not match, a power will be given according to the experience list in the WGCNA FAQ. I don’t like this experience power very much. <font color = blue>If you find that the R <sup>2</sup> value corresponding to experience power given by the software lower than your setting Threshold </font>,<font color = purple><b> please select a customized power based on the SFT plot.</b></i></font></font>'),
                       radioButtons(
                         inputId = "PowerTorF",
                         label = "Power type",
                         choices = c("Recommended","Customized"),
                         selected = "Recommended"
                       ),

                       sliderInput(
                         inputId = "PowerSelect",
                         label = "Final Power Selection",
                         min = 1,
                         max = 33,
                         value = 6
                       )
                     )
                 ),
                 mainPanel(
                   fluidPage(
                     #### output field
                     actionButton("toggleSidebar2",
                                  "Toggle sidebar"),
                     tabsetPanel(
                       tabPanel(title = "Select Power",height = "500px",width = "100%",
                                icon = icon("th"),
                                actionButton("Startsft","Start analysis"),
                                htmlOutput("powerout"),
                                jqui_resizable(
                                  plotOutput("sftplot")
                                ),
                                textInput(inputId = "width2",
                                          label = "width",
                                          value = 10),
                                textInput(inputId = "height2",
                                          label = "height",
                                          value = 10),
                                actionButton("adjust2","Set fig size"),
                                downloadButton("downfig2","Download")


                       ),
                       tabPanel(title = "Information of sft table",height = "500px",width = "100%",
                                icon = icon("table"),
                                DT::dataTableOutput("sfttbl")
                       ),
                       tabPanel(title = "scale free estimate",height = "500px",width = "100%",
                                icon = icon("chart-bar"),
                                actionButton("Startcheck","Check Scale-free"),
                                jqui_resizable(
                                  plotOutput("sfttest")
                                ),
                                textInput(inputId = "width3",
                                          label = "width",
                                          value = 10),
                                textInput(inputId = "height3",
                                          label = "height",
                                          value = 10),
                                actionButton("adjust3","Set fig size"),
                                downloadButton("downfig3","Download")

                       )## tabPanel
                     )## tabsetPanel
                   )## fluidPage
                 )
               )
             ),##tabPanel
             tabPanel(
               useShinyjs(),
               title = "Module-net",
               icon = icon("play-circle"),
               sidebarLayout(
                 div(id = "Sidebar3",
                   sidebarPanel(
                     width = 2,
                     sliderInput(
                       inputId = "minMsize",
                       label = "min Module Size",min = 0,max = 200,value = 30
                     ),
                     sliderInput(
                       inputId = "mch",
                       label = "module cuttree height",
                       min = 0, max = 1, value = 0.25
                     ),
                     textInput(inputId = "blocksize",
                               label = "select max blocksize",
                               value = 5000),
                     p("MaxBlockSize, The default was 5000, 4GB memory could handle 8000-10000 genes, for 16GB memory you can select at most of 24000 genes in one block, 32GB should be enough for 30000-40000. Try to keep all selected genes in one block",
                       style = "color: #7a8788;font-size: 12px; font-style:Italic")
                   )
                 ),
                 mainPanel(
                   fluidPage(
                     actionButton("toggleSidebar3",
                                  "Toggle sidebar"),
                     tabsetPanel(
                       tabPanel(
                         title = "Cluster",height = "500px",width = "100%",
                         icon = icon("table"),
                         actionButton("Startnet","Start"),
                         jqui_resizable(
                           plotOutput("cluster")
                         ),
                         textInput(inputId = "width4",
                                   label = "width",
                                   value = 10),
                         textInput(inputId = "height4",
                                   label = "height",
                                   value = 10),
                         actionButton("adjust4","Set fig size"),
                         downloadButton("downfig4","Download"),

                         br(),
                         br(),
                         tableOutput("m2num")
                       ),
                       tabPanel(
                         title = "Eigengene adjacency heatmap",height = "500px",width = "100%",
                         icon = icon("th"),
                         jqui_resizable(
                           plotOutput("eah")
                         ),
                         textInput(inputId = "width5",
                                   label = "width",
                                   value = 10),
                         textInput(inputId = "height5",
                                   label = "height",
                                   value = 10),
                         actionButton("adjust5","Set fig size"),
                         downloadButton("downfig5","Download")
                       ),
                       tabPanel(
                         title = "TOMplot",height = "500px",width = "100%",
                         icon = icon("th-large"),
                         HTML('<font size = 2.5 color = #7a8788><i>TOMplot can be very resource-intensive for large gene sets. Run it only when needed.</i></font>'),
                         actionButton("StartTOMplot", "Plot TOMplot"),
                         jqui_resizable(
                           plotOutput("tomplot")
                         ),
                         textInput(inputId = "width9",
                                   label = "width",
                                   value = 20),
                         textInput(inputId = "height9",
                                   label = "height",
                                   value = 15),
                         actionButton("adjust9","Set fig size"),
                         downloadButton("downfig9","Download")
                       ),
                       tabPanel(
                         title = "Gene to module",height = "500px",width = "100%",
                         icon = icon("table"),
                         DT::dataTableOutput("g2m"),
                         downloadButton("downtbl2","download")
                       )
                     )

                   )
                 )
               )
             ),##tabPanel
             tabPanel(
               useShinyjs(),
               title = "Module-trait",
               icon = icon("star-of-david"),
               sidebarLayout(
                 div(id = "Sidebar4",
                   sidebarPanel(
                     width = 2,
                     fileInput(
                       inputId = "traitData",
                       label = "Upload expression matrix",
                       accept = c(".txt",".csv",".xls")
                     ),
                     colourpicker::colourInput(inputId = "colormin",
                                               label = "Minimum",
                                               value = "blue"),
                     colourpicker::colourInput(inputId = "colormid",
                                               label = "Middle",
                                               value = "white"),
                     colourpicker::colourInput(inputId = "colormax",
                                               label = "Maxmum",
                                               value = "red"),
                     textInput(
                       inputId = "xangle",
                       label = "x axis label angle",
                       value = 0
                     ),
                     actionButton("starttrait","Start analysis")
                   )
                 ),
                 mainPanel(
                   actionButton("toggleSidebar4",
                                "Toggle sidebar"),
                   fluidPage(
                     tabsetPanel(
                       tabPanel(
                         title = "Module to trait",height = "500px",width = "100%",
                         icon = icon("ht"),
                         jqui_resizable(
                           plotOutput("mtplot")
                         ),
                         textInput(inputId = "width6",
                                   label = "width",
                                   value = 10),
                         textInput(inputId = "height6",
                                   label = "height",
                                   value = 10),
                         actionButton("adjust6","Set fig size"),
                         downloadButton("downfig6","Download")
                       ),
                       tabPanel(
                         title = "Module-trait matrix",height = "500px",width = "100%",
                         icon = icon("table"),
                         DT::dataTableOutput("traitmat"),
                         DT::dataTableOutput("traitp")
                       ),
                       tabPanel(
                         title = "eigengene-based connectivities,KME",height = "500px",width = "100%",
                         icon = icon("table"),
                         DT::dataTableOutput("KME"),
                         downloadButton("downtbl3","download")
                       )
                     )
                   )
                 )
               )
             ),##tabPanel
             tabPanel(
               useShinyjs(),
               title = "Interested module",
               icon = icon("broom"),
               sidebarLayout(
                 div(id = "sidebar5",
                   sidebarPanel(
                     width = 2,
                     selectInput(
                       inputId = "strait",
                       label = "Select traits",
                       choices = c("select a trait","trait2","..."),
                       selected = "select a trait",
                       multiple = F
                     ),
                     selectInput(
                       inputId = "smodule",
                       label = "Select module",
                       choices = c("red","black","..."),
                       selected = "red",
                       multiple = F
                     ),
                     actionButton("InterMode","Start Analysis")
                   )
                 ),
                 mainPanel(
                   actionButton("toggleSidebar5",
                                "Toggle sidebar"),
                   fluidPage(
                     tabsetPanel(
                       tabPanel(
                         title = "GS-Connectivity",height = "500px",width = "100%",
                         icon = icon("chart-line"),
                         jqui_resizable(
                           plotOutput("GSCon")
                           ),
                         textInput(inputId = "width7",
                                   label = "width",
                                   value = 10),
                         textInput(inputId = "height7",
                                   label = "height",
                                   value = 10),
                         actionButton("adjust7","Set fig size"),
                         downloadButton("downfig7","Download")
                         ),
                       tabPanel(
                         title = "Heatmap",height = "500px",width = "100%",
                         icon = icon("buromobelexperte"),
                         jqui_resizable(
                           plotOutput("heatmap")
                           ),
                         textInput(inputId = "width8",
                                   label = "width",
                                   value = 10),
                         textInput(inputId = "height8",
                                   label = "height",
                                   value = 10),
                         actionButton("adjust8","Set fig size"),
                         downloadButton("downfig8","Download")
                       ),
                       tabPanel(
                         title = "MM vs GS all",height = "500px",width = "100%",
                         icon = icon("chart-line"),
                         jqui_resizable(
                           plotOutput("GSMM.all")
                         ),
                         textInput(inputId = "width10",
                                   label = "width",
                                   value = 10),
                         textInput(inputId = "height10",
                                   label = "height",
                                   value = 10),
                         actionButton("adjust10","Set fig size"),
                         downloadButton("downfig10","Download")
                       )
                     )
                   )
                 )
               ),
             ),##tabPanel
             tabPanel(
               useShinyjs(),
               title = "hub gene",
               icon = icon("star"),
               sidebarLayout(
                 div(id = "sidebar6",
                     sidebarPanel(
                       width = 2,
                       selectInput(
                         inputId = "hubtrait",
                         label = "Select trait",choices = c("Select a trait","..."),multiple = F
                       ),
                       selectInput(
                         inputId = "hubmodule",
                         label = "Select module",choices = c("Select a module","..."),multiple = F
                       ),
                       actionButton("starthub","Start Analysis")
                     )
                 ),
                 mainPanel(
                   actionButton("toggleSidebar6",
                                "Toggle sidebar"),
                   fluidPage(
                     tabsetPanel(
                       tabPanel(
                         title = " choose Top Hub In Each Module (Not recommended)",
                         icon = icon("sad-cry"),
                         DT::dataTableOutput("cthub")
                       ),
                       tabPanel(
                         title = "By kME and GS (Yes!)",
                         icon = icon("smile"),
                         sliderInput(
                           inputId = "kMEcut",
                           label = "cutoff of  absolute value of kME",
                           min = 0,max = 1,step = 0.01,
                           value = 0.5
                         ),
                         sliderInput(
                           inputId = "GScut",
                           label = "cutoff of  absolute value of GS",
                           min = 0,max = 1,step = 0.01,
                           value = 0.5
                         ),
                         DT::dataTableOutput("kMEhub"),
                         downloadButton("downtbl4","download")
                       ),
                       tabPanel(
                         title = "Cytoscape output",
                         icon = icon("dna"),
                         textInput(
                           inputId = "threshold",
                           label = "weight threshold",
                           value = 0.02
                         ),
                         actionButton("threadd","choose the threshold"),
                         DT::dataTableOutput("edgeFile"),
                         DT::dataTableOutput("nodeFile"),
                         downloadButton("downtbl5","download edgefile"),
                         downloadButton("downtbl6","download nodefile")
                       )
                     )
                   )
                 )
               )
             )##tabPanel

  )## navbarPage
)## UI

readExpressionMatrix <- function(upload) {
  ext <- tolower(tools::file_ext(upload$name))
  first_lines <- readLines(upload$datapath, n = 10, warn = FALSE, encoding = "UTF-8")
  first_lines <- sub("^\\ufeff", "", first_lines)
  first_line <- first_lines[nzchar(first_lines)][1]
  if(is.na(first_line)) {
    first_line <- ""
  }
  countDelimiter <- function(delimiter) {
    matches <- gregexpr(delimiter, first_line, fixed = TRUE)[[1]]
    if(identical(matches, -1L)) {
      return(0)
    }
    length(matches)
  }
  sep <- if(ext == "csv") {
    ","
  } else if(ext %in% c("txt", "tsv", "tab", "xls")) {
    "\t"
  } else if(countDelimiter("\t") >= countDelimiter(",")) {
    "\t"
  } else {
    ","
  }

  encodings <- c("UTF-8-BOM", "UTF-8", "GB18030", "")
  last_error <- NULL
  for(file_encoding in encodings) {
    read_args <- list(file = upload$datapath,
                      sep = sep,
                      header = TRUE,
                      stringsAsFactors = FALSE,
                      check.names = FALSE,
                      quote = "\"",
                      comment.char = "",
                      na.strings = c("NA", "", "NaN"))
    if(nzchar(file_encoding)) {
      read_args$fileEncoding <- file_encoding
    }
    table <- tryCatch(do.call(read.table, read_args),
                      error = function(e) {
                        last_error <<- e
                        NULL
                      })
    if(!is.null(table)) {
      stripUtf8BomAndWhitespace <- function(x) {
        utf8_bom <- rawToChar(as.raw(c(0xEF, 0xBB, 0xBF)))
        trimws(gsub(utf8_bom, "", x, fixed = TRUE))
      }
      names(table) <- stripUtf8BomAndWhitespace(names(table))
      if(ncol(table) > 0 && (is.character(table[[1]]) || is.factor(table[[1]]))) {
        table[[1]] <- stripUtf8BomAndWhitespace(as.character(table[[1]]))
      }

      first_column_is_gene_id <- ncol(table) > 1 &&
        (is.character(table[[1]]) || is.factor(table[[1]]))

      if(first_column_is_gene_id) {
        gene_ids <- table[[1]]
        if(any(is.na(gene_ids) | !nzchar(gene_ids))) {
          stop(paste(
            "Gene ID column contains empty or missing gene IDs.",
            expression_format_message
          ), call. = FALSE)
        }
        duplicated_gene_ids <- unique(gene_ids[duplicated(gene_ids)])
        if(length(duplicated_gene_ids) > 0) {
          stop(paste(
            "Gene ID column contains duplicate gene IDs:",
            paste(utils::head(duplicated_gene_ids, 5), collapse = ", "),
            "Please make gene IDs unique before uploading.",
            expression_format_message
          ), call. = FALSE)
        }
      }

      expression_column_names <- names(table)[-1]
      for(column_name in expression_column_names) {
        original_column <- table[[column_name]]
        converted_column <- suppressWarnings(as.numeric(as.character(original_column)))
        if(any(is.na(converted_column) & !is.na(original_column))) {
          stop(paste(
            "Expression column", shQuote(column_name),
            "contains non-numeric expression values. Please check the non-numeric expression values and upload again.",
            expression_format_message
          ), call. = FALSE)
        }
      }
      rownames(table) <- NULL
      return(table)
    }
  }
  stop(last_error)
}

server <- function(input, output, session){
  observeEvent(input$toggleSidebar, {
    shinyjs::toggle(id = "Sidebar")
  })
  observeEvent(input$toggleSidebar2, {
    shinyjs::toggle(id = "Sidebar2")
  })
  observeEvent(input$toggleSidebar3, {
    shinyjs::toggle(id = "Sidebar3")
  })
  observeEvent(input$toggleSidebar4, {
    shinyjs::toggle(id = "Sidebar4")
  })
  observeEvent(input$toggleSidebar5, {
    shinyjs::toggle(id = "Sidebar5")
  })
  observeEvent(input$toggleSidebar6, {
    shinyjs::toggle(id = "Sidebar6")
  })
  data <- reactive({
    file1 <- input$ExpMat
    if(is.null(file1)){return()}
    readExpressionMatrix(file1)
  })

  inputMatrixValidation <- function(input_data) {
    if(is.null(input_data)) {
      return(list(valid = FALSE, message = "No expression matrix was uploaded."))
    }
    if(ncol(input_data) < 2) {
      return(list(
        valid = FALSE,
        message = "Your expression matrix has fewer than 2 columns after reading. Expected format: first column is gene IDs; all remaining columns are numeric sample expression values. Please check whether the file delimiter is correct and upload a gene-by-sample table again."
      ))
    }

    gene_ids <- as.character(input_data[[1]])
    if(any(is.na(gene_ids) | !nzchar(trimws(gene_ids)))) {
      return(list(
        valid = FALSE,
        message = "The first column must contain non-empty gene IDs. Please remove blank gene IDs and upload the file again."
      ))
    }
    duplicated_gene_ids <- unique(gene_ids[duplicated(gene_ids)])
    if(length(duplicated_gene_ids) > 0) {
      return(list(
        valid = FALSE,
        message = paste0(
          "The first column must contain unique gene IDs. Duplicate IDs include: ",
          paste(utils::head(duplicated_gene_ids, 5), collapse = ", "),
          ". Please make gene IDs unique before uploading."
        )
      ))
    }

    expression_data <- input_data[-1]
    numeric_columns <- vapply(expression_data, is.numeric, logical(1))
    if(!all(numeric_columns)) {
      return(list(
        valid = FALSE,
        message = paste0(
          "Every column after gene_id must contain numeric sample expression values. Non-numeric columns include: ",
          paste(names(expression_data)[!numeric_columns], collapse = ", "),
          "."
        )
      ))
    }
    if(any(is.na(as.matrix(expression_data)))) {
      return(list(
        valid = FALSE,
        message = "Columns after gene_id contain blank, non-numeric, or otherwise invalid NA expression values. Please clean the expression matrix and upload the file again."
      ))
    }

    list(valid = TRUE, message = "There is no problem with your expression matrix format. Expected format confirmed: first column contains unique non-empty gene IDs and every remaining column contains numeric sample expression values without invalid NA values. Please proceed to the next step.")
  }

  expressionMatrixForProcessing <- function(input_data) {
    validation <- inputMatrixValidation(input_data)
    if(!validation$valid) {
      return(NULL)
    }
    expression_data <- input_data[-1]
    rownames(expression_data) <- as.character(input_data[[1]])
    expression_data
  }

  output$Inputcheck = renderUI({
    if(is.null(data())){return()}
    validation <- inputMatrixValidation(data())
    if(validation$valid) {
      HTML(paste0('<font color = red><b>Congratulations!,</b></font> ', validation$message))
    } else {
      HTML(paste0('<font color = blue><b>Sorry!</b></font> ', validation$message))
    }
  })
  
  ## set WGCNA threads

  select = dplyr::select
  wgcna_thread <- reactiveValues(current = 20)
  output$show.wgcna.threads = renderText({
    paste(as.character(input$wgcna.threads), "threads")
  })
  output$current.wgcna.threads = renderText({
    paste(as.character(wgcna_thread$current), "threads")
  })
  observeEvent(input$set.wgcna.threads, {
    thread_count <- as.integer(input$wgcna.threads)
    enableWGCNAThreads(nThreads = thread_count)
    wgcna_thread$current <- thread_count
    message(paste("WGCNA threads set to", thread_count))
  })
  

  
  ## count number
  fmt = reactive({
    input$format
  })
  observe({
    if(fmt() == "count") {
      updateSelectInput(session, "method1",choices = c(VST = "varianceStabilizingTransformation",
                                                       logCPM = "lgcpm"))
      updateTextInput(session,"RCcut",value = 10)
    } else if(fmt() == "FPKM") {
      updateSelectInput(session, "method1",choices = c(FPKM_TPM = "rawFPKM",
                                                       logFPKM_TPM = "lgFPKM"))
      updateTextInput(session,"RCcut",value = 1)
    }

  })
  mtd = reactive({
    input$method1
  })

  sampP = reactive({
    as.numeric(input$SamPer)
  })
  rccutoff = reactive({
    as.numeric(input$RCcut)
  })
  GNC = reactive({
    as.numeric(input$remain)
  })
  cutmethod = reactive({
    input$CutMethod
  })
  ## set reactiveValues
  exp.ds<-reactiveValues(data=NULL)
  downloads <- reactiveValues(data = NULL)
  observeEvent(
    input$action1,
    {
      if(is.null(data())){return()}
      processing_data <- expressionMatrixForProcessing(data())
      if(is.null(processing_data)) {return()}
      exp.ds$table = data.frame()
      exp.ds$table2 = data.frame()
      exp.ds$param = list()
      exp.ds$layout = as.character(input$treelayout)

      output$filter1 = renderUI({
        input$action1
        p_mass = c("Processing step1, remove very low expressed genes",
                   paste("Processing step2, pick out high variation genes via",cutmethod()))
        withProgress(
          message = "Raw data normlization",
          value = 0,{
            incProgress(1/2,detail = p_mass[1])
            exp.ds$table = getdatExpr(rawdata = processing_data,
                                      RcCutoff = rccutoff(),samplePerc = sampP(),
                                      datatype = fmt(),method = mtd())
            Sys.sleep(0.1)
            if(is.null(exp.ds$table) || nrow(exp.ds$table) == 0 || ncol(exp.ds$table) == 0) {
              exp.ds$table2 = data.frame()
              return(HTML(paste0('<font color = blue><b>No genes remain after the first filter.</b></font><br/>',
                                 'For FPKM/TPM input, please use a smaller Expression Cutoff (default 1) or lower the Sample percentage.')))
            }
            reserved_gene_num <- min(GNC(), nrow(exp.ds$table))
            gene_num_cut <- 1 - reserved_gene_num/nrow(exp.ds$table)
            incProgress(1/2,detail = p_mass[2])
            exp.ds$table2 = getdatExpr2(datExpr = exp.ds$table,
                                        GeneNumCut = gene_num_cut,cutmethod = cutmethod())
            if(is.null(exp.ds$table2) || nrow(exp.ds$table2) == 0 || ncol(exp.ds$table2) == 0) {
              return(HTML('<font color = blue><b>No genes remain after the second filter.</b></font> Please reduce the Reserved genes Num. or adjust the filter settings.'))
            }
            exp.ds$param = getsampleTree(exp.ds$table2,layout = exp.ds$layout)
            Sys.sleep(0.1)
          }
        )
        isolate(HTML(paste0('<font color = red> <b>After filtered by conditions:</b> </font>removing all features with expression/count less than <font color = red><b>',rccutoff(),'</b></font> in more than <font color = red> <b>',100*sampP(),'% </b></font> of the samples','<br/>',
                            '<font color = red> <b>Remaining Gene Numbers: </b> </font>',nrow(exp.ds$table),'<br/>',
                            '<font color = red> <b>After filtered by conditions:</b> </font>Genes with <font color = red><b>',cutmethod(),'</b></font> ranked top <font color = red> <b>',min(GNC(), nrow(exp.ds$table)),' </b></font> of all expressed genes','<br/>',
                            '<font color = red> <b>Remaining Gene Numbers: </b> </font>',ncol(exp.ds$table2))))
      })
    }
  )


  ## summary num
  firstFilterGeneIds <- reactive({
    if(is.null(exp.ds$table)){return(character(0))}
    if(is.null(exp.ds$table2)){return(character(0))}
    table_rows <- rownames(exp.ds$table)
    table_cols <- colnames(exp.ds$table)
    if(!is.null(table_rows) && length(intersect(table_rows, colnames(exp.ds$table2))) > 0) {
      return(table_rows)
    }
    table_cols
  })

  cleanedGeneIds <- reactive({
    if(is.null(exp.ds$table)){return(character(0))}
    if(is.null(exp.ds$table2)){return(character(0))}
    setdiff(firstFilterGeneIds(), colnames(exp.ds$table2))
  })

  cleanedGeneExpression <- function(gene_ids) {
    if(length(gene_ids) == 0){return(data.frame())}
    if(all(gene_ids %in% rownames(exp.ds$table))) {
      exp.ds$table[gene_ids, , drop = FALSE]
    } else {
      as.data.frame(t(exp.ds$table[, gene_ids, drop = FALSE]))
    }
  }

  cleanedGenesForInput <- function(gene_ids) {
    if(length(gene_ids) == 0){return(data.frame())}
    if(all(gene_ids %in% rownames(exp.ds$table))) {
      add_back <- as.data.frame(t(exp.ds$table[gene_ids, , drop = FALSE]))
      colnames(add_back) <- gene_ids
      add_back
    } else {
      exp.ds$table[, gene_ids, drop = FALSE]
    }
  }
  hasNonEmptyDimensions <- function(input_data) {
    !is.null(input_data) && length(dim(input_data)) == 2 && all(dim(input_data) > 0)
  }

  inputPreviewSummary <- function(input_data, preview_label) {
    first_column_name <- if(ncol(input_data) > 0) colnames(input_data)[1] else "None"
    first_column_unique <- if(ncol(input_data) > 0) length(unique(input_data[[1]])) else 0
    na_count <- sum(is.na(input_data))
    paste0(
      preview_label,
      " | Rows read: ", nrow(input_data),
      " | Columns read: ", ncol(input_data),
      " | First column: ", first_column_name,
      " | Unique first-column IDs: ", first_column_unique,
      " | NA values: ", na_count
    )
  }

  output$Inputbl = DT::renderDataTable({
    input_data <- data()
    if(is.null(input_data)){return()}

    raw_preview <- head(input_data, 10)
    preview_message <- inputPreviewSummary(
      input_data,
      "Raw uploaded matrix preview. Confirm the TSV delimiter and first-column gene IDs before filtering."
    )

    validation <- inputMatrixValidation(input_data)
    if(!validation$valid) {
      preview_message <- paste(preview_message, validation$message, sep = " | ")
    }

    if(hasNonEmptyDimensions(exp.ds$table2)) {
      filtered_preview <- head(as.data.frame(exp.ds$table2), 10)
      filtered_message <- inputPreviewSummary(
        filtered_preview,
        paste0(
          "Filtered expression matrix preview after clicking Update information! ",
          "The raw upload preview keeps gene_id as the first column; this filtered preview is transposed for downstream processing, so rows correspond to samples and columns correspond to filtered genes."
        )
      )
      filtered_message <- paste(
        filtered_message,
        paste0("Raw upload contained ", nrow(input_data), " gene rows and ", ncol(input_data), " columns including gene_id."),
        sep = " | "
      )
      return(DT::datatable(
        filtered_preview,
        caption = htmltools::tags$caption(style = "caption-side: top; text-align: left;", filtered_message),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = TRUE
      ))
    }

    preview_message <- paste(
      preview_message,
      "Filtered matrix has not been generated yet; click Update information! to run filtering.",
      sep = " | "
    )
    DT::datatable(
      raw_preview,
      caption = htmltools::tags$caption(style = "caption-side: top; text-align: left;", preview_message),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$washedGenes = DT::renderDataTable({
    if(is.null(data())){return()}
    if(!inputMatrixValidation(data())$valid) {return()}
    washed_ids <- cleanedGeneIds()
    if(length(washed_ids) == 0){
      return(data.frame(gene_id = character(0)))
    }
    washed_data <- as.data.frame(cleanedGeneExpression(washed_ids))
    washed_data <- cbind(gene_id = rownames(washed_data), washed_data)
    rownames(washed_data) <- NULL
    DT::datatable(
      washed_data,
      selection = list(mode = "multiple", target = "row"),
      filter = "top",
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })

  observeEvent(input$addWashedGenes, {
    if(is.null(exp.ds$table)){return()}
    if(is.null(exp.ds$table2)){return()}
    selected_rows <- input$washedGenes_rows_selected
    if(length(selected_rows) == 0){
      output$addWashedGenesInfo <- renderUI({
        HTML('<font color = blue><b>No cleaned genes selected.</b></font>')
      })
      return()
    }
    washed_ids <- cleanedGeneIds()
    selected_genes <- washed_ids[selected_rows]
    selected_genes <- selected_genes[!is.na(selected_genes)]
    selected_genes <- setdiff(selected_genes, colnames(exp.ds$table2))
    if(length(selected_genes) == 0){
      output$addWashedGenesInfo <- renderUI({
        HTML('<font color = blue><b>Selected genes are already in Input.</b></font>')
      })
      return()
    }
    exp.ds$table2 <- cbind(exp.ds$table2, cleanedGenesForInput(selected_genes))
    exp.ds$param <- getsampleTree(exp.ds$table2, layout = exp.ds$layout)
    exp.ds$sft <- NULL
    exp.ds$power <- NULL
    exp.ds$net <- NULL
    exp.ds$tomDiss <- NULL
    exp.ds$tomGeneTree <- NULL
    output$addWashedGenesInfo <- renderUI({
      HTML(paste0('<font color = red><b>', length(selected_genes), '</b></font> cleaned genes were added back to Input. Please rerun downstream WGCNA analyses.'))
    })
  })
  ## sample tree
  output$clustPlot = renderPlot({
    if(is.null(data())){return()}
    if(length(which(is.na(data()))) != 0) {return()}
    validate(
      need(!is.null(exp.ds$table2) &&
             length(dim(exp.ds$table2)) == 2 &&
             nrow(exp.ds$table2) > 0 &&
             ncol(exp.ds$table2) > 0,
           "请先成功完成 Update information / filtering"),
      need(!is.null(exp.ds$param$sampleTree),
           "请先成功完成 Update information / filtering")
    )
    plot(exp.ds$param$sampleTree,main = "Sample clustering to detect outlier", sub = "", xlab = "")
  })
  ## download sample tree

  rscut = reactive({
    as.numeric(input$CutoffR)
  })

  observeEvent(
    input$Startsft,
    {
      if(is.null(exp.ds$table2)){return()}
      exp.ds$sft = list()
      output$powerout = renderUI({
        sft_mess = c("pick soft threshold in processing ...",
                     "Finish.")
        withProgress(message = 'SFT selection', value = 0,
                     expr = {
                       for (i in 1:2) {
                         incProgress(1/2, detail = sft_mess[i] )
                         if (i == 1) {
                           exp.ds$sft = getpower(datExpr = exp.ds$table2,rscut = rscut())
                         } else {
                           return()
                         }

                       }
                     })
        isolate(HTML(paste0('<font color = red> <b>The power recommended by WGCNA is:</b> </font><font color = bule><b>',exp.ds$sft$power,'</b></font> ','<br/>',
                            '<font color = pink> <i>If all power values lower than the R square threshold which you set, it means that the power value is an empirical value. At this time, you need to infer a power value based on the results on this plot and check whether it can form a scale-free network. </i> </font>')))
      })
    }
  )


  ## outsft
  output$sftplot = renderPlot({
    if(is.null(exp.ds$table2)){return()}
    input$Startsft
    if(length(exp.ds$sft) == 0){return()}
    exp.ds$sft$plot
  })
  ## outtbl
  output$sfttbl = DT::renderDataTable({
    if(is.null(exp.ds$table2)){return()}
    input$Startsft
    if(length(exp.ds$sft) == 0){return()}
    as.data.frame(exp.ds$sft$sft)
  })
  ## test sft
  pcus = reactive({
    as.numeric(input$PowerSelect)
  })
  PowerTorF = reactive({
    input$PowerTorF
  })
  observeEvent(
    input$Startcheck,
    {
      if(is.null(exp.ds$table2)){return()}
      if(is.null(exp.ds$sft)){return()}
      sftcheck_mess = c("Checking scale free network ...",
                        "Finish.")
      exp.ds$power = exp.ds$sft$power
      exp.ds$cksft = list()
      withProgress(message = 'SFT selection', value = 0,
                   expr = {
                     for (i in 1:2) {
                       incProgress(1/2, detail = sftcheck_mess[i] )
                       if (i == 1) {
                         if(PowerTorF() == "Recommended"){
                           exp.ds$power = exp.ds$sft$power
                           exp.ds$cksft = powertest(power.test = exp.ds$sft$power,datExpr = exp.ds$table2,nGenes = exp.ds$param$nGenes)
                         } else if (PowerTorF() == "Customized"){
                           exp.ds$power = pcus()
                           exp.ds$cksft = powertest(power.test = pcus(),datExpr = exp.ds$table2,nGenes = exp.ds$param$nGenes)
                         }
                       } else {
                         return()
                       }

                     }
                   })
    }
  )

  output$sfttest = renderPlot({
    if(is.null(exp.ds$sft)){return()}
    input$Startcheck
    if(length(exp.ds$cksft ) == 0){return()}
    exp.ds$cksft
  })
  mms = reactive({
    as.numeric(input$minMsize)
  })
  mch = reactive({
    as.numeric(input$mch)
  })
  blocksize = reactive({
    as.numeric(input$blocksize)
  })
  observeEvent(
    input$Startnet,
    {
      if(is.null(exp.ds$table2)){return()}
      if(is.null(exp.ds$power)){return()}
      exp.ds$netout = getnetwork(datExpr = exp.ds$table2,power = exp.ds$power,
                                 minModuleSize = mms(),mergeCutHeight = mch(),maxBlocksize = blocksize())
      exp.ds$nSamples = nrow(exp.ds$table2)
      exp.ds$net = exp.ds$netout$net
      exp.ds$moduleLabels = exp.ds$netout$moduleLabels
      exp.ds$moduleColors = exp.ds$netout$moduleColors
      exp.ds$MEs_col = exp.ds$netout$MEs_col
      exp.ds$MEs = exp.ds$netout$MEs
      exp.ds$Gene2module = exp.ds$netout$Gene2module
    }
  )
  output$cluster = renderPlot({
    input$Startnet
    if(is.null(exp.ds$net)){return()}
    plotDendroAndColors(exp.ds$net$dendrograms[[1]], exp.ds$moduleColors[exp.ds$net$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
  })
  output$m2num = renderTable({
    input$Startnet
    if(is.null(exp.ds$net)){return()}
    table(exp.ds$moduleColors)
  })
  output$eah = renderPlot({
      input$Startnet
      if(is.null(exp.ds$net)){return()}
      plotEigengeneNetworks(exp.ds$MEs_col, "Eigengene adjacency heatmap",
                            marDendro = c(3,3,2,4),
                            marHeatmap = c(3,4,2,2), plotDendrograms = T,
                            xLabelsAngle = 90)
    })

  observeEvent(input$StartTOMplot, {
    if(is.null(exp.ds$net)){return()}
    if(is.null(exp.ds$table2)){return()}
    if(is.null(exp.ds$power)){return()}
    withProgress(message = 'TOMplot', value = 0, {
      incProgress(1/3, detail = "Calculating adjacency matrix")
      adjacency_matrix <- adjacency(exp.ds$table2, power = exp.ds$power)
      incProgress(1/3, detail = "Calculating TOM similarity")
      TOM <- TOMsimilarity(adjacency_matrix)
      exp.ds$tomDiss <- 1 - TOM
      incProgress(1/3, detail = "Clustering genes for TOMplot")
      exp.ds$tomGeneTree <- hclust(as.dist(exp.ds$tomDiss), method = "average")
    })
  })

  output$tomplot = renderPlot({
    input$StartTOMplot
    if(is.null(exp.ds$tomDiss)){return()}
    if(is.null(exp.ds$tomGeneTree)){return()}
    plotTOM <- exp.ds$tomDiss^7
    TOMplot(plotTOM, exp.ds$tomGeneTree, exp.ds$moduleColors, main = "Network heatmap plot, all genes")
  })

  output$g2m = DT::renderDataTable({
    input$Startnet
    if(is.null(exp.ds$net)){return()}
    exp.ds$Gene2module
  })

  phen <- reactive({
    file2 <- input$traitData
    if(is.null(file2)){return()}
    read.csv(file = file2$datapath,
               sep=",",
               header = T,
               stringsAsFactors = F)
  })


  observeEvent(
    input$starttrait,
    {
      if(is.null(phen())){return()}
      if (ncol(phen()) == 2) {
        x <- phen()
        Tcol = as.character(unique(x[,2]))
        b <- list()
        for (i in 1:length(Tcol)) {
          b[[i]] = data.frame(row.names = x[,1],
                              levels = ifelse(x[,2] == Tcol[i],1,0))
        }
        c <- bind_cols(b)
        c <- data.frame(row.names = x$name,
                        c)
        colnames(c) = Tcol
        rownames(c) = phen()[,1]
        exp.ds$phen<- c
      } else {
        exp.ds$phen = data.frame(row.names = phen()[,1],
                            phen()[,-1])
      }
      exp.ds$phen =  exp.ds$phen[match(rownames(exp.ds$table2),rownames(exp.ds$phen)),]
      exp.ds$traitout = getMt(phenotype = exp.ds$phen,MEs_col = exp.ds$MEs_col,
                              nSamples = exp.ds$nSamples,moduleColors = exp.ds$moduleColors,datExpr = exp.ds$table2)
      exp.ds$xangle = as.numeric(input$xangle)
      exp.ds$c_min = as.character(input$colormin)
      exp.ds$c_mid = as.character(input$colormid)
      exp.ds$c_max = as.character(input$colormax)
      exp.ds$modTraitCor = exp.ds$traitout$modTraitCor
      exp.ds$modTraitP = exp.ds$traitout$modTraitP
      exp.ds$textMatrix = exp.ds$traitout$textMatrix
      exp.ds$KME = getKME(datExpr = exp.ds$table2,moduleColors = exp.ds$moduleColors,MEs_col = exp.ds$MEs_col)
      exp.ds$mod_color = gsub(pattern = "^..",replacement = "",rownames(exp.ds$modTraitCor))
      exp.ds$mod_color_anno = setNames(exp.ds$mod_color,rownames(exp.ds$modTraitCor))
      exp.ds$Left_anno = rowAnnotation(
        Module = rownames(exp.ds$modTraitCor),
        col = list(
          Module = exp.ds$mod_color_anno
        ),
        show_legend = F,
        show_annotation_name = F
      )
    }
  )

    output$mtplot = renderPlot({
      input$starttrait
      if(is.null(phen())){return()}
      if(is.null(exp.ds$phen)){return()}
      Heatmap(
        matrix = exp.ds$modTraitCor,
        cluster_rows = F, cluster_columns = F,
        left_annotation = exp.ds$Left_anno,
        cell_fun = function(j,i,x,y,width,height,fill) {
          grid.text(sprintf(exp.ds$textMatrix[i,j]),x,y,gp = gpar(fontsize = 12))
        },
        row_names_side = "left",
        column_names_rot = exp.ds$xangle,
        heatmap_legend_param = list(
          at = c(-1,-0.5,0,0.5, 1),
          labels = c("-1","-0.5", "0","0.5", "1"),
          title = "",
          legend_height = unit(9, "cm"),
          title_position = "lefttop-rot"
        ),
        rect_gp = gpar(col = "black", lwd = 1.2),
        column_title = "Module-trait relationships",
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        col = colorRamp2(c(-1, 0, 1), c(exp.ds$c_min, exp.ds$c_mid, exp.ds$c_max))
      )

    })



  output$traitmat = DT::renderDataTable({
    input$starttrait
    if(is.null(phen())){return()}
    if(is.null(exp.ds$phen)){return()}
    as.data.frame(exp.ds$modTraitCor)
  })

  output$traitp = DT::renderDataTable({
    input$starttrait
    if(is.null(phen())){return()}
    if(is.null(exp.ds$phen)){return()}
    as.data.frame(exp.ds$modTraitP)
  })

  output$KME = DT::renderDataTable({
    input$starttrait
    if(is.null(phen())){return()}
    if(is.null(exp.ds$phen)){return()}
    as.data.frame(exp.ds$KME)
  })
  s_mod = reactive({
    gsub(pattern = "^..",replacement = "",rownames(exp.ds$modTraitP))
  })
  observe({
    updateSelectInput(session, "smodule",choices = s_mod())
  })

  s_trait = reactive({
    colnames(exp.ds$modTraitP)
  })
  observe({
    updateSelectInput(session, "strait",choices = s_trait())
  })

  observeEvent(
    input$InterMode,
    {
      if(is.null(phen())){return()}
      if(is.null(exp.ds$phen)){return()}
      exp.ds$GSout = getMM(datExpr = exp.ds$table2,MEs_col = exp.ds$MEs_col,nSamples = exp.ds$nSamples,corType = "pearson")
      exp.ds$MM = exp.ds$GSout$MM
      exp.ds$MMP = exp.ds$GSout$MMP
      exp.ds$sml = as.character(input$smodule)
      exp.ds$st = as.character(input$strait)
      exp.ds$Heatmap = moduleheatmap(datExpr = exp.ds$table2,MEs = exp.ds$MEs_col,which.module = exp.ds$sml,
                                     moduleColors = exp.ds$moduleColors)

    }
  )

  output$GSCon = renderPlot({
    input$InterMode
    if(is.null(exp.ds$st)){return()}
    if(is.null(exp.ds$sml)){return()}
    getverboseplot(datExpr = exp.ds$table2,module = exp.ds$sml,pheno = exp.ds$st,MEs = exp.ds$MEs_col,
                   traitData = exp.ds$phen,moduleColors = exp.ds$moduleColors,
                   geneModuleMembership = exp.ds$MM,nSamples = exp.ds$nSamples)
  })

  output$heatmap = renderPlot({
    input$InterMode
    if(is.null(exp.ds$st)){return()}
    if(is.null(exp.ds$sml)){return()}
    exp.ds$Heatmap
  })

  output$GSMM.all = renderPlot({
    input$InterMode
    if(is.null(exp.ds$st)){return()}
    if(is.null(exp.ds$sml)){return()}
    MMvsGSall(which.trait = exp.ds$st,
              traitData = exp.ds$phen,
              datExpr = exp.ds$table2,
              moduleColors = exp.ds$moduleColors,
              geneModuleMembership = exp.ds$MM,
              MEs = exp.ds$MEs_col,
              nSamples = exp.ds$nSamples)
  })

  observe({
    updateSelectInput(session, "hubmodule",choices = s_mod())
  })

  observe({
    updateSelectInput(session, "hubtrait",choices = s_trait())
  })

  observeEvent(
    input$starthub,
    {
      exp.ds$hubml = as.character(input$hubmodule)
      exp.ds$hubt = as.character(input$hubtrait)
      exp.ds$kMEcut = as.numeric(input$kMEcut)
      exp.ds$GScut = as.numeric(input$GScut)
      print(exp.ds$hubml)
      exp.ds$hub.all = hubgenes(datExpr = exp.ds$table2,
                                mdl = exp.ds$hubml,
                                power = exp.ds$power,
                                trt = exp.ds$hubt,
                                KME = exp.ds$KME,
                                GS.cut = exp.ds$GScut,
                                kME.cut =exp.ds$kMEcut,
                                datTrait = exp.ds$phen,
                                g2m = exp.ds$Gene2module
                                )
    }
  )

  observeEvent(
    input$threadd,
    {
      exp.ds$threshold = as.numeric(input$threshold)
      exp.ds$cyt = cytoscapeout(datExpr = exp.ds$table2,
                                power = exp.ds$power,module = exp.ds$hubml,
                                moduleColors = exp.ds$moduleColors,
                                threshold = exp.ds$threshold)
    }
  )
  # checkAdjMat
  output$cthub = DT::renderDataTable({
    input$starthub
    if(is.null(exp.ds$hubml)){return()}
    if(is.null(exp.ds$hubt)){return()}
    exp.ds$hub.all$hub1
  })

  output$kMEhub = DT::renderDataTable({
    input$starthub
    if(is.null(exp.ds$hubml)){return()}
    if(is.null(exp.ds$hubt)){return()}
    if(is.null(exp.ds$kMEcut)){return()}
    if(is.null(exp.ds$GScut)){return()}
    exp.ds$hub.all$hub3
  })

  output$edgeFile = DT::renderDataTable({
    input$threadd
    if(is.null(exp.ds$hubml)){return()}
    if(is.null(exp.ds$threshold)){return()}
    exp.ds$cyt[[1]]
  })

  output$nodeFile = DT::renderDataTable({
    input$threadd
    if(is.null(exp.ds$hubml)){return()}
    if(is.null(exp.ds$threshold)){return()}
    exp.ds$cyt[[2]]
  })

# download ----------------------------------------------------------------


  observeEvent(
    input$adjust1,
    {
      downloads$width1 <- as.numeric(input$width1)
      downloads$height1 <-  as.numeric(input$height1)
    }
  )
  observeEvent(
    input$adjust2,
    {
      downloads$width2 <- as.numeric(input$width2)
      downloads$height2 <-  as.numeric(input$height2)
    }
  )
  observeEvent(
    input$adjust3,
    {
      downloads$width3 <- as.numeric(input$width3)
      downloads$height3 <-  as.numeric(input$height3)
    }
  )
  observeEvent(
    input$adjust4,
    {
      downloads$width4 <- as.numeric(input$width4)
      downloads$height4 <- as.numeric(input$height4)
    }
  )
  observeEvent(
    input$adjust5,
    {
      downloads$width5 <- as.numeric(input$width5)
      downloads$height5 <-  as.numeric(input$height5)
    }
  )
  observeEvent(
    input$adjust6,
    {
      downloads$width6 <- as.numeric(input$width6)
      downloads$height6 <-  as.numeric(input$height6)
    }
  )
  observeEvent(
    input$adjust7,
    {
      downloads$width7 <- as.numeric(input$width7)
      downloads$height7 <-  as.numeric(input$height7)
    }
  )
  observeEvent(
    input$adjust8,
    {
      downloads$width8 <- as.numeric(input$width8)
      downloads$height8 <-  as.numeric(input$height8)
    }
  )
  observeEvent(
    input$adjust9,
    {
      downloads$width9 <- as.numeric(input$width9)
      downloads$height9 <-  as.numeric(input$height9)
    }
  )
  observeEvent(
    input$adjust10,
    {
      downloads$width10 <- as.numeric(input$width10)
      downloads$height10 <-  as.numeric(input$height10)
    }
  )
  library(ape)
  output$downfig1 = downloadHandler(
    filename = function() {
      "01.SampleCluster.nwk"
    },
    content = function(file) {
      validate(
        need(!is.null(exp.ds$param$tree),
             "请先成功完成 Update information / filtering")
      )
      write.tree(phy = exp.ds$param$tree,file = file)
    }
  )
  output$downfig2 = downloadHandler(
    filename = function() {
      "02.SftResult.pdf"
    },
    content = function(file) {
      ggsave(plot = exp.ds$sft$plot,filename = file,width = downloads$width2,height = downloads$height2)
    }
  )
  output$downfig3 = downloadHandler(
    filename = function() {
      "03.CheckSft.pdf"
    },
    content = function(file) {
      ggsave(plot = exp.ds$cksft,filename = file,width = downloads$width3,height = downloads$height3)
    }
  )
  output$downfig4 = downloadHandler(
    filename = function() {
      "04.ClusterDendrogram.pdf"
    },
    content = function(file) {
      pdf(file = file,width = downloads$width4, height = downloads$height4 )
      plotDendroAndColors(exp.ds$net$dendrograms[[1]], exp.ds$moduleColors[exp.ds$net$blockGenes[[1]]],
                          "Module colors",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05)
      dev.off()
    }
  )
  output$downfig5 = downloadHandler(
    filename = function() {
      "05.EigengeneadJacencyHeatmap.pdf"
    },
    content = function(file) {
      pdf(file = file,width = downloads$width5, height = downloads$height5)
      plotEigengeneNetworks(exp.ds$MEs_col, "Eigengene adjacency heatmap",
                            marDendro = c(3,3,2,4),
                            marHeatmap = c(3,4,2,2), plotDendrograms = T,
                            xLabelsAngle = 90)
      dev.off()
    }
  )
  output$downfig9 = downloadHandler(
    filename = function() {
      "06.TOMplot.pdf"
    },
    content = function(file) {
      pdf(file = file,width = downloads$width9, height = downloads$height9)
      plotTOM <- exp.ds$tomDiss^7
      TOMplot(plotTOM, exp.ds$tomGeneTree, exp.ds$moduleColors, main = "Network heatmap plot, all genes")
      dev.off()
    }
  )
  output$downfig6 = downloadHandler(
    filename = function() {
      "06.Module2Trait.pdf"
    },
    content = function(file) {
      pdf(file = file,width = downloads$width6, height = downloads$height6)

      print(Heatmap(
        matrix = exp.ds$modTraitCor,
        cluster_rows = F, cluster_columns = F,
        left_annotation = exp.ds$Left_anno,
        cell_fun = function(j,i,x,y,width,height,fill) {
          grid.text(sprintf(exp.ds$textMatrix[i,j]),x,y,gp = gpar(fontsize = 12))
        },
        row_names_side = "left",
        column_names_rot = exp.ds$xangle,
        heatmap_legend_param = list(
          at = c(-1,-0.5,0,0.5, 1),
          labels = c("-1","-0.5", "0","0.5", "1"),
          title = "",
          legend_height = unit(9, "cm"),
          title_position = "lefttop-rot"
        ),
        rect_gp = gpar(col = "black", lwd = 1.2),
        column_title = "Module-trait relationships",
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        col = colorRamp2(c(-1, 0, 1), c(exp.ds$c_min, exp.ds$c_mid, exp.ds$c_max))
      ))

      dev.off()
    }
  )
  output$downfig7 = downloadHandler(
    filename = function() {
      paste0("07.GS",exp.ds$sml,"-",exp.ds$st,"-Connectivity.pdf")
    },
    content = function(file) {
      pdf(file = file,width = downloads$width7, height = downloads$height7)
     print(getverboseplot(datExpr = exp.ds$table2,module = exp.ds$sml,pheno = exp.ds$st,MEs = exp.ds$MEs_col,
                     traitData = exp.ds$phen,moduleColors = exp.ds$moduleColors,
                     geneModuleMembership = exp.ds$MM,nSamples = exp.ds$nSamples))
      dev.off()
    }
  )
  output$downfig8 = downloadHandler(

    filename = function() {
      paste0("08.",exp.ds$sml,"-",exp.ds$st,"MEandGeneHeatmap.pdf")
    },
    content = function(file) {
      pdf(file = file,width = downloads$width8, height = downloads$height8)
      print(exp.ds$Heatmap)
      dev.off()
    }
  )
  output$downfig10 = downloadHandler(

    filename = function() {
      "09.GSvsMM.all.pdf"
    },
    content = function(file) {
      pdf(file = file,width = downloads$width8, height = downloads$height8)
      print(MMvsGSall(which.trait = exp.ds$st,
                traitData = exp.ds$phen,nSamples = exp.ds$nSamples,
                datExpr = exp.ds$table2,
                moduleColors = exp.ds$moduleColors,
                geneModuleMembership = exp.ds$MM,MEs = exp.ds$MEs_col))
      dev.off()
    }
  )
  output$downtbl2 = downloadHandler(

    filename = function() {
      if(is.null(exp.ds$net)){return()}
      "01.Gene2Module.xls"
    },
    content = function(file) {
      write.table(x = exp.ds$Gene2module,file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtbl3 = downloadHandler(
    filename = function() {
      "02.KMEofAllGenes.xls"
    },
    content = function(file) {
      write.table(x = exp.ds$KME,file = file,sep = "\t",row.names = T,quote = F)
    }
  )
  output$downtbl4 = downloadHandler(
    filename = function() {
      paste0("03.",exp.ds$hubml,"-",exp.ds$hubt,"hubgene_by_GS_MM.xls")
    },
    content = function(file) {
      write.table(x = exp.ds$hub.all$hub3,file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtbl5 = downloadHandler(
    filename = function() {
      paste0("04.",exp.ds$hubml,".edge.xls")
    },
    content = function(file) {
      write.table(x = exp.ds$cyt[[1]],file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  output$downtbl6 = downloadHandler(
    filename = function() {
      paste0("04.cyt",exp.ds$hubml,".node.xls")
    },
    content = function(file) {
      write.table(x = exp.ds$cyt[[2]],file = file,sep = "\t",row.names = F,quote = F)
    }
  )
  shinyServer(function(input, output, session){
    session$onSessionEnded(function() {
      stopApp()
    })
  })

}

shinyServer(function(input, output, session){
  session$onSessionEnded(function() {
    stopApp()
  })
})

shinyApp(ui,server)
