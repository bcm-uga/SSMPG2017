tabPanel("Challenge",
         fluidRow(
           box(title = strong("Challenge"),
               solidHeader = TRUE,
               width = 6, 
               status = "danger", 
               height = 250,
               selectInput("challenge", "Challenge",
                           choices = c(1, 2),
                           selected = 1),
               radioButtons("dataset", "Dataset",
                            choices = c("Training set", "Evaluation set"),
                            selected = "Training set")
           ),
           
           box(title = strong("Submission"),
               solidHeader = TRUE,
               width = 6, 
               status = "warning", 
               height = 250,
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
               ),
               div(actionButton(inputId = "submit",
                                label = "Submit"), 
                   align = "center")
           )
           
         ),
         fluidRow(
           box(title = strong("Leaderboard"),
               solidHeader = TRUE,
               width = 12, 
               status = "success",
               tabsetPanel(
                 source(file.path("ui", "tab-barchart.R"), local = TRUE)$value,
                 source(file.path("ui", "tab-table.R"), local = TRUE)$value,
                 id = "conditionedPanels"
               )  
           )
         )
)
