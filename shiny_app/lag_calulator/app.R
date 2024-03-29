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
library(minpack.lm)
library(shinydisconnect)
library(remotes)
# miLAG package can be installed from here:
#devtools::install_github("https://github.com/bognabognabogna/microbial_lag_calulator", dependencies = TRUE, upgrade = "always", force = FALSE)
library(miLAG)
#source("~/MGG Dropbox/Bogna Smug/Projects/Quiesence/2022_Lags/GitHub/microbial_lag_calulator/R/milags_functions.R")


available.methods = list("max growth acceleration",
                         "tangent",
                         "biomass increase",
                         "parameter fitting to a model")
parameters.default = get_def_pars()
my_theme = get_theme(text_size = 22)


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
                         div("This calculator is a tool developed to estimate lag phase duration out of a growth curve data for microbial populations.
                        The user provides datasets of microbial density change in time and chooses a method, parameters, and pre-processing techniques, according to the instruction below.
                      If you find our tool useful, please cite:"),
                         tags$a(href="https://www.biorxiv.org/content/10.1101/2022.11.16.516631v1", "Opalek & Smug, 2022"),
                         br(),
                         strong("The datasets uploaded by users are not saved nor anyhow collected. "),
                         br(),
                         br(),
                         strong("INSTRUCTION"),
                         br(),
                         strong("1.	Upload your dataset."),
                         div("The accepted file formats are csv and txt. The dataset must contain two required columns as specified below:"),
                         div("(i) the first column: time (preferably in hours)"),
                         div("(ii) the second column: biomass i.e. population size (preferably in CFU/mL, otherwise in OD)"),
                         div("(iii) the third optional column: curve id (an identifier for a growth curve if dataset consists of multiple growth curves)"),
                         div("The population size is recommend to be measured by CFU/mL values instead of raw absorbance (OD),
                    This is because the correlation between CFU number and absorbance is rarely linear.
                      However, if you’re unable to provide CFU values, the calculator will also work for absorbance data (assuming OD is directly proportional to the CFU)."),
                      strong("If you decide to use OD values, remember to blank correct the data: 0 OD should mean no cells in the culture."),
                      div("After uploading your dataset, specify the column and decimal separators."),
                         div("If your dataset consists of multiple growth curves, they should be provided in long format i.e. one curve under another and identified by the curve_id provided in the third column.
                             If multiple growth curves are technical replicates, we advise to fit the lag to the averaged growth curve, 
                             otherwise lags can be fitted separately to each growth curve."),
                         strong("2. Use some data pre-processing if needed"),
                         br(),
                         strong("3.	Chose lag duration calculation method"),
                         div("
                    Within this calculator we allow for a choice of one of the four most commonly used methods of lag duration calculation:"),
                         div("BIOMASS INCREASE - an increase of biomass (CFU/ml or OD) by the user specified value from the beginning of growth (or minimal biomass value)"),
                         div("MAX GROWTH ACCELERATION - finds the point of the growth curve where the second derivative is maximal"),
                         div("TANGENT METHOD - the intersection of the initial density line and the  line tangent to the part of the curve where the growth rate is maximal"),
                         div("PARAMETER FITTING TO A MODEL - uses fitting procedures to simultaneously fit all growth curve parameters (e.g. lag phase length, maximal growth rate, and maximal population size). "),
                         strong("4.	Adjust parameters of the selected model or pre-process the data again.
                           For example, if parameter fitting is chosen it may be helpful to cut out the stationary phase data (within data pre-processing tab)"),
                         br(),
                         br(),
                         br(),
                         strong("WHICH METHOD SHOULD I CHOOSE?"),
                         imageOutput("Decission_Tree"),
                         #img(src='Decission_Tree.png'),
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
                                         conditionalPanel("input.manual == 'No'",
                                         fileInput("growth.curve.file",
                                                   "Choose CSV or TXT File (first column = time, second column = biomass i.e. CFU/ml or OD)",
                                                   multiple = FALSE,
                                                   accept = c("text/csv","text/comma-separated-values,text/plain", ".csv"),
                                                   placeholder = "No file selected, using example data:",
                                                   width = '80%')),
                         )),
                         conditionalPanel("input.manual == 'Yes'",
                          fluidRow(
                           column(12, align="center",
                                  textAreaInput("pasted.data", "Paste your tabular data from a spreadsheet:"),
                         ))),
                         conditionalPanel("input.manual == 'No'",
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
                          ),
                         fluidRow(
                           column(12, align="center",
                                  selectInput("manual",
                                              label = "Would you prefer to input data manually (paste from Excel)?",
                                              choices = c("Yes", "No"),
                                              selected = "No"))),
                         br(),
                         fluidRow(
                           column(12, align = "center",
                         selectInput("biomass_unit",
                                     label = "Biomass unit",
                                     choices = c("CFU/ml", "OD"),
                                     selected = ","))),
                         fluidRow(
                           column(12, align = "center",
                                  selectInput("separate._fit",
                                              label = "Do you want to fit lag to every curve separately or to the mean?",
                                              choices = c("Separete lags", "Fit lag to the mean"),
                                              selected = "Separate lags"))),
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
                         ),
                         mainPanel(
                           verbatimTextOutput("not_enough_many_points_message"),
                           br(),
                           h5("Growth curve data after pre-processing"),
                           uiOutput("growth.curve.data.processed.plot"),
                           textOutput("sd.explanation.message.1"),
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
                                                                        label = "NLS fitting algorithm (defaults to the best of:  bounded Levenberg-Marquardt, unbounded Levenberg-Marquardt or bounded port)",
                                                                        choices = c("auto","Levenberg-Marquardt","port", "plinear", "port", "nl2sol"),
                                                                        selected = "auto"),
                                                            br(),
                                                            numericInput("max.iter",
                                                                         label = "Max number of iterations",
                                                                         value = 100,
                                                                         min = 100,
                                                                         max = 1000,
                                                                         step = 100),
                                                            br(),
                                                            conditionalPanel(condition = "input.algorithm != 'auto'",
                                                                             selectInput("own_init_params",
                                                                                         label = "Do you want to specify \nyour own initial guesses n\for the parameters?",
                                                                                         choices = c("No", "Yes"),
                                                                                         selected = "No"),
                                                                             conditionalPanel(condition = "input.own_init_params == 'Yes'",
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
                                                                                                           step = 0.01)),
                                                            )),
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
                                           downloadButton("downloadData", "Download"),
                                           br(),
                                           tableOutput("lag.stats"),
                                           br(),
                                           h5("Note that you may be able improve the efficiency of the methods by varying the parameters or further pre-processing the data.
              For example, if parameter fitting method is chosen it may be helpful to cut out the stationary phase data (within data pre-processing tab)."),
                                           uiOutput("growth.curve.plot"),
                                           textOutput("sd.explanation.message.2"),
                                         ))),width = 8)
    )


  )
)
)





server <- shinyServer(function(input, output, session) {
  output$Decission_Tree <- renderImage({

    list(src = "www/Decission_Tree.png",
         #width = "120%",
         height = "200%"
    )

  }, deleteFile = F)


  isValid_input <- reactive({ !is.null(input$method) & !is.null(growth.curve.file)})

  model.params = reactive({
    pars = list(model = input$model,
                n0_method = input$N0.method,
                tangent_method = input$tangentmethod,
                threshold = input$threshold,
                curve_points = as.numeric(input$n.points.in.curve),
                init_gr_rate = NULL,
                init_lag = NULL,
                algorithm = input$algorithm,
                max_iter = input$max.iter)

    if (input$own_init_params == 'Yes') {
      pars$init_gr_rate = input$init.growth.rate
      pars$init_lag = input$init.lag
    }
    return(pars)

  })

  output$errortext <- renderText({
    if(isValid_input()){ invisible(NULL)
    }else{ #
      "Parameter values are outside a reasonable range.
             Note: They must be numeric and non-negative"
    }
  })
  
  

################################################## Tab 2: read data and display it ################################################## 

  growth.curve.data = reactive({
    message = ""
    if (is.null(input$growth.curve.file)) {
      data.path = "R/chosen_exampless_curves1to3_biomass_ml_for_app.csv" #"R/example_data.csv"
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
    }  else {
      out = tryCatch(
        {
          if (input$manual == "Yes") {
            data = data.table::fread(paste(input$pasted.data, collapse = "\n")) %>% as.data.frame()
          } else {
            # in case there is a header ignore all non-numeric values
            # make sure the decimal is always "." so that we can convert it to string and back to numeric
            data =  read.table(data.path,
                               header = FALSE,
                               dec = input$dec,
                               sep = column.separator,
                               stringsAsFactors = FALSE)
          }
          
          if (ncol(data) < 2) {
            error("At least two columns are needed (time and cell count)")
          } else if (ncol(data) == 2) {
            data = data %>% mutate(curve_id = "provided_growth_curve")
          }
          names(data) = c("time", "biomass", "curve_id")
          data = mutate_at(data, c("time", "biomass"), function(x) as.numeric(sub(input$dec, ".", as.character(x), fixed = TRUE)))
          
          data = data %>% 
            filter(!is.na(time) & !is.na(biomass) & !is.na(biomass)) %>%
            select(time, biomass, curve_id)
          
          return(list(data=data, message = message))
        },
        error=function(cond) {
          message = "Can't read in the data. Are the separators specified correctly?"
          # Choose a return value in case of error
          return(list(data=data.frame(time = numeric(0), biomass = numeric(0), curve_id = character(0)),
                      message = message))
        }
      )
      
      return(out)
    }
  })
  

  output$growth.curve.data.table <- DT::renderDataTable(growth.curve.data()$data,
                                                        rownames= FALSE,
                                                        options = list(searching = FALSE))
  

  
  
  ################################################## Tab 3: data plot and preprocessing ################################################## 
  log10_transform_param = reactive({
    if (input$biomass_unit == "CFU/ml") {
      return(TRUE)
    } else if (input$biomass_unit == "OD") {
      return(FALSE)
    } else {
      error("Incorrect biomass_unit")
    }
    })
    
    y.label = reactive({
      if (log10_transform_param()) {
        return(paste0("log10(", input$biomass_unit, ")"))
      } else {
        return("OD")
      } 
    })
    
    
  
    # for CFU/mL work on log scale
    growth.curve.data.sd = reactive({
      data = growth.curve.data()$data
      if (input$separate._fit == "Fit lag to the mean") {
        if (log10_transform_param()) { 
          data = data %>%
            mutate(biomass = log10(biomass))
        }
        sd.data = data %>% 
          group_by(time) %>%
          summarise(
            SD = sd(biomass, na.rm = TRUE)) %>%
          ungroup()
      } else {
        sd.data = data %>% distinct(time) %>% mutate(SD = NA)
      }
      return(sd.data)
    })
    
    output$sd.explanation.message.1 <- renderText({
      if (input$separate._fit == "Fit lag to the mean") {
        "Error bars represent standard devation."
      } else {
        ""
      }
    })
    output$sd.explanation.message.2 <- renderText({
      if (input$separate._fit == "Fit lag to the mean") {
      "Lag has been fitted to the mean growth curve.\n
      Error bars represent standard devation."
      } else {
        ""
      }
    })
    
    data.for.lag.calculation = reactive({
      if (input$separate._fit == "Fit lag to the mean") {
        data = growth.curve.data()$data %>%
          group_by(time) %>%
          summarise(biomass = mean(biomass, na.rm = TRUE)) %>%
          ungroup() %>%
          mutate(curve_id = "Mean growth curve")
      } else {
        data = growth.curve.data()$data
      }
      return(data)
    })
    
    num.curves = reactive({data.for.lag.calculation() %>% distinct(curve_id) %>% nrow()})
    
    
    
    growth.curve.data.cut = reactive({
      data = data.for.lag.calculation()
      if (input$cut_data_flag == 'yes') {
        data = cut_the_data(data, input$cut_max_time)
      } 
      return(data)
    })
    
    
    growth.curve.data.processed = reactive({
      if (input$smooth_data_flag == 'yes') {
        data.smooth = smooth_data(growth.curve.data.cut(), '3RS3R') #input$smooth_method)
      } else {
        data.smooth = growth.curve.data.cut()
      }
      return(data.smooth)
    })
    
    output$not_enough_many_points_message <- renderText({
      if (nrow(growth.curve.data.processed()) < 4) {
        return("Not enough many time points")
      } else {
        return("")
      }

    })
    
    output$growth.curve.data.processed.plot.raw = renderPlot({
      data.to.plot = growth.curve.data.processed() %>%
        left_join(growth.curve.data.sd())
      if (nrow(data.to.plot) > 0) {
        fig = plot_data(data_new = data.to.plot, log10_transform = log10_transform_param()) + 
          my_theme +
          scale_y_continuous() + 
          facet_grid(curve_id~ .) +
          ylab(y.label()) +
          geom_errorbar(aes(ymin = biomass - SD, ymax = biomass + SD))
      } else
        fig = ggplot() + theme_void() + 
          ylab(y.label())
      return(fig)
    })
    
    output$growth.curve.data.processed.plot =   renderUI({
      plotOutput("growth.curve.data.processed.plot.raw", height = 250*num.curves())
    })
    
  ################################################## Tab 4: Lag calculation and data plot ################################################# 
  growth.curve.data.with.lag =  reactive({
    if (nrow(growth.curve.data.processed()) > 0) {
      data.with.lag = calc_lag(data = growth.curve.data.processed(),
                                    method = input$method,
                                    pars = model.params())
    } else {
      data.with.lag = growth.curve.data.processed()
    }
    return(data.with.lag)
  })


  output$message <- renderText({
    growth.curve.data()$message
  })

  lag.value = reactive({
    if (nrow(growth.curve.data.with.lag()) > 0) {
      lag.data.table = growth.curve.data.with.lag() %>%
        distinct(lag, curve_id) %>%
        rowwise() %>%
        mutate(lag.info = round(lag,3),
               method = input$method) %>%
        ungroup() %>%
        select(curve_id, lag = lag.info, method)

    } else {
      lag.data.table = data.frame(method = input$method, lag = NA)
      return(lag.data.table)
    }
  })

  lag.stats = reactive({
    data.frame(
      lag.mean = lag.value()$lag %>% mean(na.rm = TRUE),
      lag.sd = lag.value()$lag %>% sd(na.rm = TRUE),
      lag.min = lag.value()$lag %>% min(na.rm = TRUE),
      lag.max = lag.value()$lag %>% max(na.rm = TRUE))
  })
    
  
  output$downloadData <- downloadHandler(
    filename = function() {
      "lag.values.csv"
    },
    content = function(file) {
      write.csv(lag.value(), file)
    }
  )
  
  
  output$lag.value = renderTable(lag.value(),
                                 rownames= FALSE)
  output$lag.stats = renderTable(lag.stats(),
                                 rownames= FALSE)



  
  output$plot = renderPlot({
    if (nrow(growth.curve.data.with.lag()) > 0) {
      data.to.plot = growth.curve.data.with.lag() %>%
        left_join(growth.curve.data.sd())
      # seems like the NA fail here, try without the unnecessary tangent lines etc
      #write.csv(data.to.plot, "/Users/bognasmug/MGG Dropbox/Bogna Smug/Projects/Quiesence/2022_Lags/GitHub/microbial_lag_calulator/shiny_app/lag_calulator/R/test_data/test.csv", row.names = FALSE)
      fig = plot_lag_fit(data.to.plot, log10_transform = log10_transform_param()) + 
        my_theme +
        ylab(y.label())  +
        geom_errorbar(aes(ymin = biomass - SD, ymax = biomass + SD))
    } else {
      fig = ggplot() + theme_void() + ylab(y.label())
    }
    return(fig)
  })

  output$growth.curve.plot =   renderUI({
    plotOutput("plot", height = 250*num.curves())
  })
  
  observeEvent(input$disconnect, {
    session$close()
  })
  session$onSessionEnded(function() {
    stopApp()
  })

})

# Run the application
shinyApp(ui = ui, server = server)
