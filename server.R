server <- function(input, output) {
  
  n.challenge <- reactive({
    n <- 0
    if (input$tab == "challenge"){
      n <- 1
    } else if (input$tab == "challenge_2"){
      n <- 2
    }
    return(n)
  })
  
  gt <- reactive({
    tmp <- all.gt[[1]]
    if (input$dataset == "Training set" && !input$switch){
      tmp <- all.gt[[1]]
    } else if (input$dataset == "Evaluation set" && !input$switch){
      tmp <- all.gt[[2]]
    } else if (input$dataset == "Training set" && input$switch){
      tmp <- all.gt.reg[[1]]$region.selected 
    } else if (input$dataset == "Evaluation set" && input$switch){
      tmp <- all.gt.reg[[2]]$region.selected 
    }
    return(tmp)
  })
  
  submission <- reactiveValues(x = NULL)
  
  observe({
    req(input$subm)
    submission$x <- scan(input$subm$datapath)
  })
  
  # submission <- reactive({
  #   if (is.null(input$subm)) return(NULL)
  #   scan(input$subm$datapath)
  # })
  
  source(file.path("server", "uiOutput.R"), local = TRUE)$value
  source(file.path("server", "users.R"), local = TRUE)$value
  source(file.path("server", "submission.R"), local = TRUE)$value
  source(file.path("server", "update.R"), local = TRUE)$value
  
}
