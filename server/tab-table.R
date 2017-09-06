observeEvent(input$summary, {
  if (input$dataset == "Training set") {
    db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = "db.sqlite3")
    df <- RSQLite::dbReadTable(db, "submission")
    RSQLite::dbDisconnect(db)
    output$summary <- DT::renderDataTable({
      df %>% 
        filter(challenge == n.challenge()) %>%
        mutate(date = as.POSIXct(date)) %>%
        group_by(name) %>%  
        summarise(rank = which.max(date),
                  Power = power[rank],
                  FDR = fdr[rank],
                  G_score = score[rank],
                  Date = date[rank],
                  Methods = methods[rank],
                  N = n()) %>%
        arrange(desc(G_score), N, Date) %>%
        mutate(Date = as.character(Date), 
               rank = NULL, 
               candidates = NULL,
               Power = round(Power, digits = 5),
               FDR = round(FDR, digits = 5),
               G_score = round(G_score, digits = 5)) %>%
        rename("G-score" = "G_score", 
               "Team" = "name") %>%
        DT::datatable(escape = FALSE, rownames = FALSE, 
                      options = list(initComplete = DT::JS(
                                       "function(settings, json) {",
                                       "$(this.api().table().header()).css({'background-color': '#49c', 'color': '#fff'});",
                                       "}")
                                     )
                      ) %>%
        formatStyle(
          "Team",
          color = "lightseagreen")
    })
  } else if (input$dataset == "Evaluation set") {
    showModal(typeModal())  
  } 
})
