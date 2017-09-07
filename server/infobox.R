output$challengeBox <- renderValueBox({
  valueBox(paste0("#", n.challenge()),
           "Challenge",
           icon = icon("trophy"))
})

output$contestantBox <- renderValueBox({
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = "db.sqlite3")
  submission.df <- RSQLite::dbReadTable(db, "submission")
  RSQLite::dbDisconnect(db)
  submission.df %>%
    filter(challenge == n.challenge()) %>%
    dplyr::pull(name) %>%
    unique() %>%
    length() %>%
    valueBox("contestants",
             icon = icon("users"),
             color = "red")
})

observeEvent(input$tab, {
  output$leaderBox <- renderValueBox({
    db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = "db.sqlite3")
    submission.df <- RSQLite::dbReadTable(db, "submission")
    RSQLite::dbDisconnect(db)
    
    submission.df %>%
      filter(challenge == n.challenge()) %>%
      mutate(score = 1,
             date = as.POSIXct(date)) %>%
      group_by(name) %>%
      summarise(rank = which.max(date),
                G_score = score[rank],
                Date = date[rank],
                N = n()) %>%
      arrange(desc(G_score), N, Date) %>%
      pull(name) %>%
      head(n = 1) %>%
      valueBox("Current leader",
               icon = icon("heart"),
               color = "olive")
  })
})
