server <- function(input, output) {
  
  all.gt <- readRDS("anssnp.rds")
  all.gt.reg <- NULL
  all.gt.reg[[1]] <- readRDS("ansreg1a.rds")
  all.gt.reg[[2]] <- readRDS("ansreg2a.rds")
  
  n.challenge <- reactive({
    req(input$tab)
    if (input$tab == "challenge"){
      n <- 1
    } else if (input$tab == "challenge_2"){
      n <- 2
    } else {
      n <- NULL
    }
    return(n)
  })
  
  gt <- reactive({
    tmp <- NULL
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
    
    positions <- readRDS("positions.rds")
    if (input$tab == "challenge") {
      req(input$subm)
      if (input$dataset == "Training set") {
        submission$x <- positions[[1]][as.integer(scan(input$subm$datapath))]
      } else if (input$dataset == "Evaluation set") {
        submission$x <- positions[[2]][as.integer(scan(input$subm$datapath))]
      }
    } else if (input$tab == "challenge_2") {
      req(input$subm_2)
      submission$x <- as.integer(scan(input$subm_2$datapath))
    }
    
  })

  source(file.path("server", "uiOutput.R"), local = TRUE)$value
  source(file.path("server", "users.R"), local = TRUE)$value
  source(file.path("server", "submission.R"), local = TRUE)$value
  source(file.path("server", "update.R"), local = TRUE)$value

}
