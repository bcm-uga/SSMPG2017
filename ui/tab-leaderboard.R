tabPanel(title = strong("Leaderboard"),
         fixedRow(
           column(10, 
                  plotlyOutput("leaderboard"),
                  actionButton(inputId = "display",
                               label = "Display/Update")
           )
         ), value = 2
)