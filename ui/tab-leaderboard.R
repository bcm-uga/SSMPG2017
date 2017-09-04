tabPanel(title = strong("Leaderboard"),
         fixedRow(
           column(10, 
                  plotlyOutput("leaderboard")
           )
         ), value = 2
)