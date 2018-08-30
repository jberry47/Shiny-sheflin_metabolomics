library(reshape2)
library(grid)
library(stringr)
library(scales)
library(gridExtra)
library(shiny)
library(shinyBS)
library(shinydashboard)
library(DT)
library(shinyjs)
library(plyr)
library(ggplot2)
library(lme4)
library(RSQLite)
library(FactoMineR)
library(factoextra)
library(ggridges)

ui <- dashboardPage(skin="black", title="Metabolomics",
                    dashboardHeader(
                      title = tagList(
                        tags$span(
                          class = "logo-lg", "Metabolomics"
                        )
                      ),
                      titleWidth = 450
                    ),
                    dashboardSidebar(
                      tags$script(HTML("$('body').addClass('sidebar-mini');")),
                      width = 150,
                      sidebarMenu(
                        menuItem("Overview", tabName = "overview"),
                        menuItem("Get Started",tabName = "get_started")
                      )
                    ),
                    dashboardBody(
                      fluidRow(
                        tags$head(tags$style("#container * {display: inline;}")),
                        tags$style(HTML("
      .tabbable > .nav > li[class=active]    > a {background-color: #444444; color:white}
      .multicol{
      -webkit-column-count: 4; /* Chrome, Safari, Opera */
      -moz-column-count: 4; /* Firefox */
      column-count: 4;
      }
      .twocol{
      -webkit-column-count: 2; /* Chrome, Safari, Opera */
      -moz-column-count: 2; /* Firefox */
      column-count: 2;
      }
      .warning { 
      color: red;
      }"
                        )),
                        tabItem(tabName = "get_started",
                                box(width=10,title = "Import File",solidHeader = T,status = 'success',collapsible = TRUE,
                                    fileInput("metab_file", "Choose file",
                                              multiple = F,
                                              accept = c(".csv")),
                                    uiOutput("metab_go_ui")
                                ),
                                uiOutput("meta_parse"),
                                uiOutput("pca_output")
                        )
                      )
                    )
)

server <- function(input, output){
  options(shiny.maxRequestSize=2000*1024^2)
  output$metab_go_ui <- renderUI({
    b <- c(input$metab_file$name)
    if(length(b) == 1){
      actionButton("metab_merge","Merge Data")
    }
  })
  
  
  metab <- reactiveValues(data=NULL)
  
  observeEvent(input$metab_merge,{
    id <- showNotification(h3("Importing file..."),duration = NULL)
    metab$data <- read.csv(input$metab_file$datapath,header = T, stringsAsFactors = T)
    removeNotification(id)
  })
  output$meta_parse <- renderUI({
    if(!is.null(metab$data)){
      box(style = "overflow-y:scroll",width=12,title = "Subsetting Data",solidHeader = T,status = 'success',collapsible = TRUE,
          dataTableOutput("meta_config"),
          hr(),
          actionButton("use_config","Use subset")
      )   
    }
  })
  output$meta_config <- renderDataTable({
    datatable(metab$data[,1:11],rownames = F,selection = "none",filter = 'top',options = list(sDom  = '<"top">lrt<"bottom">ip',pageLength = 10, autoWidth = TRUE))
  })

  meta_sub <- reactiveValues(data = NULL)
  observeEvent(input$use_config,{
    rows <- input$meta_config_rows_all
    meta_sub$data <- metab$data[rows,]
    
  })
  
  output$pca_output <- renderUI({
    if(!is.null(meta_sub$data)){
      box(style = "overflow-y:scroll",width = 5,title = "PCA plots",solidHeader = T,status = 'success',collapsible = TRUE,
          selectInput("pca_color_by",width = 300,label = "Color By",choices = colnames(meta_sub$data[1:11]),selected = "genotype"),
          plotOutput("pca_plot_ind")
      )
    }
  })
  
  pca_df <- reactive({
    comps <- meta_sub$data
    pca <- PCA(comps[,12:ncol(comps)],graph = F)
    pca_df <- data.frame("Treatment"=comps[,input$pca_color_by],
                         "PC1"=pca$ind$coord[,1],
                         "PC2"=pca$ind$coord[,2])
    varexp <- signif(c(pca$eig[1,2],pca$eig[2,2]),3)
    list(pca_df,varexp)
  })
  
  output$pca_plot_ind <- renderPlot({
    pca_df <- pca_df()
    ggplot(data=pca_df[[1]], aes(PC1,PC2))+
      geom_point(aes(color=Treatment),alpha=0.6,size=2)+
      stat_ellipse(aes(fill=Treatment,color=Treatment),geom = "polygon",alpha=0.25)+
      xlab(paste("PC1 (",pca_df[[2]][1],"%)",sep = ""))+
      ylab(paste("PC2 (",pca_df[[2]][2],"%)",sep = ""))+
      geom_vline(xintercept = 0,linetype="dashed")+
      geom_hline(yintercept = 0,linetype="dashed")+
      theme_bw()+
      theme(strip.background=element_rect(fill="gray50"),
            strip.text.x=element_text(size=14,color="white"),
            strip.text.y=element_text(size=14,color="white"))+
      theme(axis.text = element_text(size = 14),
            axis.title= element_text(size = 18))+
      theme(axis.ticks.length=unit(0.2,"cm"))+
      theme(panel.border = element_rect(colour = "gray60", fill=NA, size=1,linetype = 1))
  })
  
}


shinyApp(ui, server)





