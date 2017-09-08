observeEvent(input$display_region, {
  if (input$tab == "challenge") {
    output$barchart_region <- renderPlot({
      db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = "db.sqlite3")
      df <- RSQLite::dbReadTable(db, "submission")
      RSQLite::dbDisconnect(db)
      gt_reg <- all.gt.reg[[1]]$region.selected
      df %>% 
        filter(challenge == n.challenge()) %>%
        mutate(r = sapply(stringr::str_split(regions, ", "), FUN = function(X) {mean(gt_reg %in% as.integer(X))}),
               p = sapply(stringr::str_split(regions, ", "), FUN = function(X) {mean(as.integer(X) %in% gt_reg)}),
               g = sqrt(r * p),
               fdr = 1 - p,
               power = r,
               score = sapply(g, FUN = function(X) {`if`(is.nan(X), 0, X)}),
               date = as.POSIXct(date), 
               candidates = NULL,
               regions = NULL) %>%
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
        as.data.frame() %>%
        ggplot(aes(x = reorder(name, G_score), 
                   y = G_score,
                   label = reorder(name, G_score))) +
        coord_flip() +
        ylim(0, 1) +
        geom_bar(stat = "identity", 
                 fill = "#32ab60", 
                 alpha = 0.6, 
                 color = "#32ab60") +
        geom_label(fontface = "bold", color = "steelblue", hjust = 1, vjust = 0.5) +
        theme_bw() +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank())
    })
  }
  
})
