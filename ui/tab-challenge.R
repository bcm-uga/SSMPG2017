tabPanel("Challenge",
         
         fluidRow(
           box(width = 6, status = "primary", height = 200,
               selectInput("challenge", "Challenge",
                           choices = c(1, 2),
                           selected = 1),
               radioButtons("dataset", "Dataset",
                            choices = c("Training set", "Evaluation set"),
                            selected = "Training set")
           ),
           
           box(width = 6, status = "warning", height = 200,
               fileInput(
                 "subm",
                 div("Choose submission file",
                     div(tags$a(href = "https://drive.google.com/uc?export=download&id=0B9o4VIJJSodfOTNOelRKNm9MQ2M",
                                h6("download example file"),
                                target = "_blank")
                     )
                 ),
                 multiple = FALSE,
                 accept = c(
                   '.csv',
                   '.txt'
                 )
               )
           )
           
         ),
         fluidRow(
           box(width = 12, status = "success",
               tabsetPanel(
                 source(file.path("ui", "tab-leaderboard.R"), local = TRUE)$value,
                 id = "conditionedPanels"
               )  
           )
         )
)
