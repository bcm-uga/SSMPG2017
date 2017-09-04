output$leaderboard <- renderPlotly({
  n.chall <- r.challenge()$N
  df %>% 
    filter(challenge == n.chall) %>%
    plot_ly(x = ~score, 
            y = ~reorder(team, score),
            type = 'bar', 
            orientation = 'h',
            color = ~reorder(team, -score),
            marker = list(line = list(width = 1))) %>%
    layout(xaxis = list(title = " "), yaxis = list(title = " ", showticklabels = FALSE))
})
