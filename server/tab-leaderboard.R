observeEvent(input$display, {
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
      mutate(Date = as.character(Date), rank = NULL)
    
    type <- input$dataset
    if (input$dataset == "Training set") {
      df %>% 
        plot_ly(x = ~G_score, 
                y = ~reorder(name, G_score),
                type = 'bar', 
                orientation = 'h',
                color = ~reorder(name, -G_score),
                marker = list(line = list(width = 1))) %>%
        layout(xaxis = list(title = " "), yaxis = list(title = " ", showticklabels = FALSE)) 
    } else if (input$dataset == "Evaluation set") {
      
    } 
    
  })
})
