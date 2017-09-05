typeModal <- function() {
  modalDialog(
    span('Leaderboard not available for the evaluation set.'),
    footer = tagList(
      actionButton("gotit", "Got it!")
    )
  )
}

# Show modal when button is clicked.
observeEvent(input$gotit, {
  removeModal()
})

observeEvent(input$display, {
  if (input$dataset == "Training set") {
  output$leaderboard <- renderPlotly({
    db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = "db.sqlite3")
    df <- RSQLite::dbReadTable(db, "submission")
    RSQLite::dbDisconnect(db)
    df <- df %>% 
      filter(challenge == 1) %>%
      mutate(date = as.POSIXct(date)) %>%
      group_by(name) %>%  
      summarise(rank = which.max(date),
                Power = power[rank],
                FDR = fdr[rank],
                G_score = score[rank],
                Date = date[rank],
                N = n()) %>%
      arrange(desc(G_score), N, Date) %>%
      mutate(Date = as.character(Date), rank = NULL) %>%
        plot_ly(x = ~G_score, 
                y = ~reorder(name, G_score),
                type = 'bar', 
                orientation = 'h',
                color = ~reorder(name, -G_score),
                marker = list(line = list(width = 1))) %>%
        layout(xaxis = list(title = " ",
                            range = c(0, 1)), 
               yaxis = list(title = " ", 
                            showticklabels = FALSE)) 
    })
  } else if (input$dataset == "Evaluation set") {
      showModal(typeModal())  
    } 
})
