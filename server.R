server <- function(input, output) {
  
  getData <- reactive({
    if (!is.null(file <- input$submission$datapath))
      parse_file(file)
  })
  
  getOptions <- reactive(input$checkbox)
  numero_challenge <- reactive(match(input$challenge, challenges))
  getKey <- reactive(input$key)
  
  # score
  score <- eventReactive(input$submit, {
    key <- getKey()
    if (key == "") {
      stop("Please provide your team key.")
    } else {
      if (is.null(getData())) {
        stop("Please provide a submission file.")
      } else {
        stop("Please provide a submission file.")
        # submit(submission = getData(), 
        #        opts = getOptions(),
        #        number_challenge = numero_challenge(), 
        #        key = key)
      }
    }
  })
  output$score <- renderTable(score())
  
  # leaderboard
  leaderboard <- eventReactive(input$leaderboard, {
    structure(update_leaderboard(numero_challenge()),
              time = lubridate::now())
  })
  output$leaderboard <- renderTable(leaderboard())
  output$time <- renderText(as.character(attr(leaderboard(), "time")))
}
