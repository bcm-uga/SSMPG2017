server <- function(input, output) {
  
  r.challenge <- reactive({
    list(N = input$challenge)
  })
  
  r.data <- reactive({
    if (!is.null(file <- input$subm$datapath))
      parse_file(file)
  })

  source(file.path("server", "users.R"), local = TRUE)$value
  source(file.path("server", "submission.R"), local = TRUE)$value
  
  # output$score <- renderTable(score())
  source(file.path("server", "tab-barchart.R"), local = TRUE)$value
  source(file.path("server", "tab-table.R"), local = TRUE)$value
  
  # leaderboard
  # leaderboard <- eventReactive(input$leaderboard, {
  #   structure(update_leaderboard(numero_challenge()),
  #             time = lubridate::now())
  # })
  # output$leaderboard <- renderTable(leaderboard())
  # output$time <- renderText(as.character(attr(leaderboard(), "time")))
}
