server <- function(input, output) {
  n.challenge <- reactive({
    if (input$tab == "challenge"){
      n.challenge <- 1
    } else if (input$tab == "challenge_2"){
      n.challenge <- 2
    }
    return(n.challenge)
  })
  source(file.path("server", "uiOutput.R"), local = TRUE)$value
  source(file.path("server", "infobox.R"), local = TRUE)$value
  source(file.path("server", "users.R"), local = TRUE)$value
  source(file.path("server", "submission.R"), local = TRUE)$value
  source(file.path("server", "tab-barchart.R"), local = TRUE)$value
  source(file.path("server", "tab-table.R"), local = TRUE)$value
  
}
