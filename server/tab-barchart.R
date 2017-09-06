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
  output$barchart <- renderPlotly({
    db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = "db.sqlite3")
    df <- RSQLite::dbReadTable(db, "submission")
    RSQLite::dbDisconnect(db)
    
    df <- df %>% 
      filter(challenge == n.challenge()) %>%
      mutate(date = as.POSIXct(date),
             candidates = NULL) %>%
      group_by(name) %>%  
      summarise(rank = which.max(date),
                Power = power[rank],
                FDR = fdr[rank],
                G_score = score[rank],
                Date = date[rank],
                N = n()) %>%
      arrange(desc(G_score), N, Date) %>%
      mutate(Date = as.character(Date), 
             rank = NULL) %>%
      as.data.frame()
    df %>%
      plot_ly(x = ~G_score, 
              y = ~reorder(name, G_score),
              type = "bar", 
              orientation = "h",
              marker = list(color = "rgba(50, 171, 96, 0.6)",
                            line = list(color = "rgba(50, 171, 96, 1.0)", width = 1)),
              hoverinfo = "text",
              text = ~paste("G-Score: ", round(G_score, digits = 2),
                            "<br> Power: ", round(Power, digits = 2),
                            "<br> FDR: ", round(FDR, digits = 2))) %>%
    layout(title = paste("<b> Challenge", n.challenge(), "-", input$dataset, "</b>"),
           xaxis = list(title = " ",
                        range = c(0, 1.1)),
           yaxis = list(title = " ",
                        showticklabels = FALSE)) %>%
    add_annotations(xref = "x1",
                    yref = "y",
                    x = pmax(0.05, df$G_score / 2),
                    y = reorder(df$name, df$G_score),
                    text = paste("<b>", reorder(df$name, df$G_score), "</b>"),
                    font = list(family = "Arial",
                                size = 14,
                                color = "rgba(255, 255, 255, 1.0)"),
                    showarrow = FALSE)
  })
  
})
