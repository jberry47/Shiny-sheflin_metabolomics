library(reshape2)
library(grid)
library(scales)
library(gridExtra)
library(shiny)
library(shinyBS)
library(shinydashboard)
library(DT)
library(shinyjs)
library(plyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(mdatools)

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
                        tabItems(
                          tabItem(tabName = "overview",
                            box(width=10,title = "Test Box",solidHeader = T,status = 'success',collapsible = TRUE,
                              p("test")
                            )
                          ),
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
)

server <- function(input, output){
  options(shiny.maxRequestSize=2000*1024^2)
  output$metab_go_ui <- renderUI({
    b <- c(input$metab_file$name)
    if(length(b) == 1){
      actionButton("metab_merge","Import Data")
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
      box(style = "overflow-y:scroll",width = 7,title = "SIMCA plots",solidHeader = T,status = 'success',collapsible = TRUE,
          column(width = 6,
              selectInput("simca_grouping_1",width = 280,label = "Grouping 1",choices = colnames(meta_sub$data[1:11]),selected = "genotype")
          ),
          column(width = 6,
                 selectInput("simca_grouping_2",width = 280,label = "Grouping 2",choices = colnames(meta_sub$data[1:11]),selected = "treatment")
          ),
          column(width = 12,
                 hr(),
                 plotOutput("simca_plot_ind"),
                 hr()
          ),
          textInput("simca_plot_name","Plot Name (ex: E17_leafOnly_cbTreatment.png)",width=400),
          downloadButton("simca_download","Download Plot"),
          br(),
          downloadButton("data_download","Download Data (tsv)")
      )
    }
  })
  
  simca_df <- reactive({
    comps <- meta_sub$data
    simca <- simca(comps[,12:ncol(comps)],"blah",scale = T)
    df <- data.frame(predict(simca, comps[,12:ncol(comps)])$scores[,1:2],stringsAsFactors = F)
    df$Grouping1 <- comps[,input$simca_grouping_1]
    df$Grouping2 <- comps[,input$simca_grouping_2]
    varexp <- signif(100*c(simca$eigenvals[1]/sum(simca$eigenvals),simca$eigenvals[2]/sum(simca$eigenvals)),4)
    list(df,varexp)
  })
  
  simca_data <- reactive({
    cbind(meta_sub$data[,1:11],simca_df())
  })
  
  
  simca_plot <- reactive({
    simca_df <- simca_df()
    ggplot(data=simca_df[[1]], aes(Comp.1,Comp.2))+
      stat_ellipse(color="gray40",fill="white",geom = "polygon",type = "norm")+
      geom_point(aes(color=Grouping2,shape=Grouping1),alpha=0.6,size=4)+
      xlab(paste("Component 1 (",simca_df[[2]][1],"%)",sep = ""))+
      ylab(paste("Component 2 (",simca_df[[2]][2],"%)",sep = ""))+
      geom_vline(xintercept = 0,linetype="dashed")+
      geom_hline(yintercept = 0,linetype="dashed")+
      theme_bw()+
      theme(panel.background = element_rect(fill = "gray95"))+
      theme(strip.background=element_rect(fill="gray50"),
            strip.text.x=element_text(size=14,color="white"),
            strip.text.y=element_text(size=14,color="white"))+
      theme(axis.text = element_text(size = 14),
            axis.title= element_text(size = 18))+
      theme(axis.ticks.length=unit(0.2,"cm"))
    })
  
  output$simca_plot_ind <- renderPlot({
    simca_plot()
  })
  
  output$simca_download <- downloadHandler(
    filename = function() {input$simca_plot_name},
    content=function(file){
      ggsave(file,simca_plot(),device = "png",width = 6,height = 5,dpi = 300)
    }
  )
  
  output$data_download <- downloadHandler(
    filename = function() {"metabolomics_simca_data.tsv"},
    content = function(file){
      head(simca_data())
      write.table(simca_data(),file,row.names = FALSE, quote = FALSE,sep = "\t")
    }
  )
  
}


shinyApp(ui, server)





