observeEvent(c(input$update, input$switch), {
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = "db.sqlite3")
  df <- RSQLite::dbReadTable(db, "submission")
  RSQLite::dbDisconnect(db)
  
  # Retrieve the latest submission for each team
  df <- df %>% 
    filter(challenge == n.challenge()) %>%
    group_by(name) %>% 
    summarise(rank = which.max(as.POSIXct(date)), 
              candidates = `if`(input$switch, regions[rank], candidates[rank]),
              date = date[rank],
              N = n()) %>%
    mutate(r = sapply(stringr::str_split(candidates, ", "), FUN = function(X) {mean(gt() %in% as.integer(X))}),
           p = sapply(stringr::str_split(candidates, ", "), FUN = function(X) {mean(as.integer(X) %in% gt())}),
           g = sqrt(r * p),
           FDR = 1 - p,
           Power = r,
           G_score = sapply(g, FUN = function(X) {`if`(is.nan(X), 0, X)}),
           candidates = NULL,
           regions = NULL,
           r = NULL,
           p = NULL,
           g = NULL) %>%
    arrange(desc(G_score), N, as.POSIXct(date)) %>%
    mutate(rank = (1:length(name))) %>%
    as.data.frame()
  
  output$summary <- DT::renderDataTable({
    df %>%
      mutate(Power = round(Power, digits = 5),
             FDR = round(FDR, digits = 5),
             G_score = round(G_score, digits = 5)) %>%
      select(rank, name, G_score, Power, FDR, N) %>%
      rename("G-score" = "G_score", "Team" = "name", "Rank" = "rank", "Submissions" = "N") %>%
      DT::datatable(escape = FALSE, rownames = FALSE, 
                    options = list(initComplete = DT::JS(
                      "function(settings, json) {",
                      "$(this.api().table().header()).css({'background-color': '#49c', 'color': '#fff'});",
                      "}"))) %>%
      formatStyle("Team", color = "steelblue")
  })
  
  output$barchart <- renderPlot({
    df %>% 
      ggplot(aes(x = reorder(name, -rank), 
                 y = G_score,
                 label = reorder(name, -rank))) +
      coord_flip() +
      ylim(0, 1) +
      geom_bar(stat = "identity", 
               fill = "#32ab60", 
               alpha = 0.6, 
               color = "#32ab60") +
      geom_label(fontface = "bold", color = "steelblue", hjust = 1.1, vjust = 0.5) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank())
  })
  
  output$challengeBox <- renderValueBox({
    valueBox(challenge_names[n.challenge()],
             "Challenge",
             icon = icon("trophy"))
  })
  
  output$contestantBox <- renderValueBox({
    df %>%
      dplyr::pull(name) %>%
      unique() %>%
      length() %>%
      valueBox("contestants",
               icon = icon("users"),
               color = "red")
  })
  
  output$leaderBox <- renderValueBox({
    df %>%
      pull(name) %>%
      head(n = 1) %>%
      valueBox("Current leader",
               icon = icon("heart"),
               color = "olive")
  })
  
  #reset("subm") 
  
})
