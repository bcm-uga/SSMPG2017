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
    as.data.frame()
  
  output$summary <- DT::renderDataTable({
    df %>%
      mutate(Power = round(Power, digits = 5),
             FDR = round(FDR, digits = 5),
             G_score = round(G_score, digits = 5)) %>%
      rename("G-score" = "G_score", "Team" = "name") %>%
      DT::datatable(escape = FALSE, rownames = FALSE, 
                    options = list(initComplete = DT::JS(
                      "function(settings, json) {",
                      "$(this.api().table().header()).css({'background-color': '#49c', 'color': '#fff'});",
                      "}"))) %>%
      formatStyle("Team", color = "steelblue")
  })
  
  output$barchart <- renderPlot({
    df %>% 
      ggplot(aes(x = reorder(name, G_score), 
                 y = G_score,
                 label = reorder(name, G_score))) +
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
    # df %>%
    #   plot_ly(x = ~G_score, 
    #           y = ~reorder(name, G_score),
    #           type = "bar", 
    #           orientation = "h",
    #           marker = list(color = "rgba(50, 171, 96, 0.6)",
    #                         line = list(color = "rgba(50, 171, 96, 1.0)", width = 1)),
    #           hoverinfo = "text",
    #           text = ~paste("G-Score: ", round(G_score, digits = 2),
    #                         "<br> Power: ", round(Power, digits = 2),
    #                         "<br> FDR: ", round(FDR, digits = 2))) %>%
    #   layout(title = paste("<b> Challenge", n.challenge(), "-", input$dataset, "</b>"),
    #          xaxis = list(title = " ",
    #                       range = c(0, 1.1)),
    #          yaxis = list(title = " ",
    #                       showticklabels = FALSE)) %>%
    #   add_annotations(xref = "x1",
    #                   yref = "y",
    #                   x = pmax(0.05, df$G_score / 2),
    #                   y = reorder(df$name, df$G_score),
    #                   text = paste("<b>", reorder(df$name, df$G_score), "</b>"),
    #                   font = list(family = "Arial",
    #                               size = 14,
    #                               color = "rgba(255, 255, 255, 1.0)"),
    #                   showarrow = FALSE)
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
})
