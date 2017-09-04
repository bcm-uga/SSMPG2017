server <- function(input, output) {
  
  r.challenge <- reactive({
    list(N = input$challenge)
  })
  
  r.data <- reactive({
    if (!is.null(file <- input$subm$datapath))
      parse_file(file)
  })

  
  source(file.path("server", "submission.R"), local = TRUE)$value
  
  # score
  # score <- eventReactive(input$submit, {
  #   key <- getKey()
  #   if (key == "") {
  #     stop("Please provide your team key.")
  #   } else {
  #     if (is.null(getData())) {
  #       stop("Please provide a submission file.")
  #     } else {
  #       stop("Please provide a submission file.")
  #       # submit(submission = getData(), 
  #       #        opts = getOptions(),
  #       #        number_challenge = numero_challenge(), 
  #       #        key = key)
  #     }
  #   }
  # })
  # output$score <- renderTable(score())
  source(file.path("server", "tab-leaderboard.R"), local = TRUE)$value
  
  # leaderboard
  # leaderboard <- eventReactive(input$leaderboard, {
  #   structure(update_leaderboard(numero_challenge()),
  #             time = lubridate::now())
  # })
  # output$leaderboard <- renderTable(leaderboard())
  # output$time <- renderText(as.character(attr(leaderboard(), "time")))
}
