#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(dplyr)
library(ggplot2)
library(deSolve)
library(DEoptim)
library(nlsMicrobio)
library(tools)
library(shinydisconnect)



source("R/lags_helper.R")
available.methods = list("max growth acceleration",
                      "tangent",
                      "biomass increase",
                      "parameter fitting to a model")
parameters.default = Get.default.parameters()
my_theme = Get.Theme(text.size = 22)


# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
    withMathJax(),
    # Application title
    titlePanel("Microbial lag phase duration calculator"),
    mainPanel(
      tabsetPanel(type = "tabs",
          tabPanel("Home",
                   br(),
                    strong("The MICROBIAL LAG PHASE DURATION CALCULATOR"),
                    div("is an online tool designed for rapid analyses of growth curves (changes in microbial density monitored in time). 
                    It allows for an easy and fast determination of the lag phase duration with the four most commonly used methods (described below) for user-specified datasets. 
                      If you find our tool useful, please cite"),
                    strong("Opalek & Smug, 2022... []."),
                    br(),
                    div("
                      In the publication, you can find a more detailed description of lag duration determination methods, together with a discussion of their advantages and biases. 
                      We also discuss possible challenges related to each of the mrthods and point out where special attention should be paid to get the most reliable outcomes."),
          
                    br(),
                    strong("INSTRUCTION"),
                   br(),
                   strong("1.	Upload your dataset."),
                   div("The accepted file formats are csv and txt. The dataset must contain two columns as specified below:"),
                    div("The first column: time (preferably in hours)"),
                    div("The second column: biomass (preferably in CFU/mL)"),
                    div("We recommend using biomass values instead of raw absorbance measurements as the correlation between biomass and absorbance is rarely linear. 
                      However, if you’re unable to provide biomass values, the calculator will also work for absorbance data.
                      After uploading your dataset, specify the column and decimal separators."),
                   strong("2. Use some data pre-processing if needed"),
                   br(),
                   strong("3.	Chose lag duration calculation method"),
                   div("
                    Within this calculator we allow for a choice of one of the four most commonly used methods of lag duration calculation:"),
                    div("MAX GROWTH ACCELERATION - finds the point of the growth curve where the second derivative is maximal"),
                    div("TANGENT METHOD - the intersection of the initial density line and the  line tangent to the part of the curve where the growth rate is maximal"),
                    div("BIOMASS INCREASE - an increase of biomass (or absorbance) by the user specified value from the beginning of growth (or minimal biomass/absorbance value)"),
                    div("PARAMETER FITTING TO A MODEL - predicted values?"),
                    strong("4.	Adjust parameters of the selected model."),
                    ),
          tabPanel("Upload growth curve", 
          tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 0px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #CCFF66;
               z-index: 105;
             }
          ")),
          tags$head(
            tags$style(HTML("
              /* this will affect only the pre elements under the class myclass */
              .myclass pre {
                color: black;
                background-color: #d6860f;
                font-weight: bolder;
              }"))
          ),
          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                           tags$div("Calculating...",id="loadmessage")),
          fluidRow(column(12, align="center",
            disconnectMessage(),
            #actionButton("disconnect", "Disconnect the app"),
            div(class = "myclass",
              verbatimTextOutput("message"),
            ),
            br(),
            fileInput("growth.curve.file", 
                      "Choose CSV or TXT File (first column = time, second column = CFU/mL)",
                         multiple = FALSE,
                         accept = c("text/csv","text/comma-separated-values,text/plain", ".csv"),
                         placeholder = "No file selected, using example data:",
                         width = '80%'),
          )),
          fluidRow(
            column(6, align="center",
            selectInput("column.separator",
                        label = "column separator",
                        choices = c(",", "tab", ";"),
                        selected = ",")),
            column(6, align="center",
                   selectInput("dec",
                               label = "decimal separator",
                               choices = c(",", "."),
                               selected = ","))),
          br(),
          h5("Provided growth curve data"),
          DT::dataTableOutput("growth.curve.data.table", width = '80%'),
          br(),
          ),
          tabPanel("Pre-process the growth curve",
                   sidebarPanel(
                   selectInput("cut_data_flag",
                               label = "Cut the data at some time?",
                               choices = c("no", "yes"),
                               selected = "no"),
                   br(),
                   conditionalPanel(condition = "input.cut_data_flag == 'yes'",
                                    numericInput("cut_max_time",
                                                label = "Where to cut the data (max time)?",
                                                value = 24,
                                                min = 0, 
                                                max = 100,
                                                step = 1)),
                  br(),
                  selectInput("smooth_data_flag",
                              label = "Smooth the data?",
                              choices = c("no", "yes"),
                              selected = "no"),
                  br(),
                  conditionalPanel(condition = "input.smooth_data_flag == 'yes'",
                                   selectInput("smooth_method",
                                                label = "Choose the Tukey's smoother",
                                                choices = c('3RS3R', '3RSS', '3R'),
                                                selected = '3RS3R')),
                   ),
                  mainPanel(
                    br(),
                    h5("Growth curve data after pre-processing"),
                    plotOutput("growth.curve.data.processed"),
                    br(),
                  )),
          tabPanel("Lag calculation", 
          fluidRow(column(12, align="center",  
            conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                             tags$div("Calculating...",id="loadmessage")),
           br(),
           sidebarPanel(
             h4("General parameters"),
             selectInput("N0.method",
                         label = "How to get the initial biomass",
                         choices = c("first.observation", "minimal.observation"),
                         selected = "first.observation"),
             br(),
             conditionalPanel(condition = "input.method == 'parameter fitting to a model'",
                              h4("Parameter related to model fitting methods"),
                              selectInput("model",
                                          label = "Model to fit",
                                          choices = c("logistic", "baranyi"),
                                          selected = "logistic"),
                              br(),
                              selectInput("algorithm",
                                          label = "NLS fitting algorithm (defaults to Gauss-Newton)",
                                          choices = c("default", "plinear", "port", "nl2sol"),
                                          selected = "default"),
                              br(),
                              numericInput("max.iter",
                                           label = "Max number of iterations",
                                           value = 100,
                                           min = 100, 
                                           max = 1000,
                                           step = 100),
                              br(),
                              h4("Initial parameter values"),
                              numericInput("init.lag",
                                           label = "Initial guess for the lag [h]",
                                           value = 3.5,
                                           min = 0, 
                                           max = 24,
                                           step = 0.1),
                              numericInput("init.growth.rate",
                                           label = "Initial guess for the max. growth rate [1/h]",
                                           value = 0.1,
                                           min = 0, 
                                           max = 1,
                                           step = 0.01)
             ),
             conditionalPanel(condition = "input.method == 'biomass increase'",
                              h4("Parameters related to biomass increase method"),
                              numericInput("threshold",
                                           value = 100,
                                           label = "Min. meaningful threshold in biomass increase",
                                           min = 0)),
             br(),
           conditionalPanel(condition = "input.method == 'tangent'",
                            h4("Parameters related to the tangent method"),
                            selectInput("tangentmethod",
                                        label = "How to find the exponential phase",
                                        choices = c("local.regression", "to.point"),
                                        selected = "local.regression")),
                            br(),
           conditionalPanel(condition = "input.method == 'tangent' && input.tangentmethod == 'local.regression'",
                            h4("Choose the number of points to be included in the regression line"),
                            selectInput("n.points.in.curve",
                                        label = "The number of points to be included",
                                        choices = c(3,5,7,9,11,13,15),
                                        selected = 3)),
                            br(),
           
           
            width = 4),
           mainPanel(
           selectInput("method",
                       label = "Lag calculation method",
                       choices = available.methods,
                       selected = available.methods[[1]]),
           br(),
           h5("Lag fitting and growth curve data plot"),
           tableOutput("lag.value"),
           br(),
           plotOutput("growth.curve.plot"),
          ))),width = 8)
          )

          
        )
    )
)





# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  
    isValid_input <- reactive({ !is.null(input$method) & !is.null(growth.curve.file)})
    
    # tangent method is called 'exponential' within the code
    selected_method = reactive({if_else(input$method == "tangent","exponential",input$method)})
    
    model.params = reactive({
      pars = list(model = input$model,
                  N0.method = input$N0.method,
                  tangent.method = input$tangentmethod,
                  threshold = input$threshold,
                  n.points.in.curve = as.numeric(input$n.points.in.curve),
                  init.growth.rate = input$init.growth.rate, 
                  init.lag = input$init.lag,
                  algorithm = input$algorithm, 
                  max.iter = input$max.iter)
      
    })
    
    output$errortext <- renderText({
        if(isValid_input()){ invisible(NULL) 
        }else{ #
            "Parameter values are outside a reasonable range.
             Note: They must be numeric and non-negative"
        }
    })
    
    growth.curve.data = reactive({
      message = ""
      if (is.null(input$growth.curve.file)) {
          data.path = "R/example_data.csv"
        } else {
          data.path =  input$growth.curve.file$datapath
        }
        data.format = file_ext(data.path)
        
        if (input$column.separator == "tab") {
          column.separator = "\t"
        } else{
          column.separator = input$column.separator
        }
        
        if (!(data.format == "csv" | data.format == "txt")) {
          message = "data format not supported"
        } else {
        out =  tryCatch(
            {
              # in case there is a header ignore all non-numeric values
              # make sure the decimal is always "." so that we can convert it to string and back to numeric
              data =  read.table(data.path,
                                 header = FALSE,
                                 dec = input$dec, 
                                 sep = column.separator,
                                 stringsAsFactors = FALSE)
              data = mutate_all(data, function(x) as.numeric(sub(input$dec, ".", as.character(x), fixed = TRUE)))
              names(data) = c("time", "biomass") 
              data = data %>% select(time, biomass) %>%
                filter(!is.na(time) & !is.na(biomass))
              return(list(data=data, message = message))
            },
            error=function(cond) {
              message = "Can't read in the data. Are the separators specified correctly?"
              # Choose a return value in case of error
              return(list(data=data.frame(time = numeric(0), biomass = numeric(0)), 
                          message = message))
            }
          )   
          
          return(out)
        } 
    })
    

    
    output$growth.curve.data.table <- DT::renderDataTable(growth.curve.data()$data, 
                                                          rownames= FALSE,
                                                          options = list(searching = FALSE))
    
    growth.curve.data.cut = reactive({
      if (input$cut_data_flag == 'yes') {
        data = Cut.The.Data(growth.curve.data()$data, input$cut_max_time) 
      } else {
        data = growth.curve.data()$data
      }
      return(data)
    })
    
    
    growth.curve.data.processed = reactive({
      if (input$smooth_data_flag == 'yes') {
        data = Smooth.Data(growth.curve.data.cut(), input$smooth_method)
      } else {
        data = growth.curve.data.cut()
      }
        return(list(data = data,
                    message = growth.curve.data()$message))
    })
    
    growth.curve.data.with.lag =  reactive({ 
      if (nrow(growth.curve.data.processed()$data) > 0) {
        data.with.lag = Calculate.Lag(data = growth.curve.data.processed()$data, 
                      method = selected_method(), 
                     pars = model.params()) %>%
        mutate(lag.calculation.method = input$method) 
      } else {
        data.with.lag = growth.curve.data.processed()$data
      }
      return(data.with.lag)
    })
    

    output$message <- renderText({ 
      growth.curve.data()$message
    })
    
    lag.value = reactive({
      if (nrow(growth.curve.data.with.lag()) > 0) {
        lag.info =  paste(round(unique(growth.curve.data.with.lag()$lag),3), " [h].")
        t = data.frame(method = input$method, lag = lag.info)
      } else {
        t = data.frame(method = input$method, lag = NA)
      return(t)
      }
    })
    
    output$lag.value = renderTable(lag.value(), 
                                    rownames= FALSE)
    

      
    output$growth.curve.plot = renderPlot({
      if (nrow(growth.curve.data.with.lag()) > 0) {
        fig = Plot.Lag.Fit(growth.curve.data.with.lag()) + my_theme
      } else 
        fig = ggplot() + theme_void()
      return(fig)
    })
    
    
    output$growth.curve.data.processed = renderPlot({
      if (nrow(growth.curve.data.processed()$data) > 0) {
        fig = Plot.Data(growth.curve.data.processed()$data) + my_theme
      } else 
        fig = ggplot() + theme_void()
      return(fig)
    })
    
    
    observeEvent(input$disconnect, {
      session$close()
    })
    
    #observeEvent(input$login, {
    #  showModal(modalDialog(
    #    title = "You have logged in.",
    #    paste0("aaa"),
    #    easyClose = TRUE,
    #    footer = NULL
    #  ))
    #})
})

# Run the application 
shinyApp(ui = ui, server = server)
