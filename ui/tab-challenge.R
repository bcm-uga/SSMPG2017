tabPanel("Challenge",
         fluidRow(
           valueBoxOutput("challengeBox"),
           valueBoxOutput("contestantBox"),
           conditionalPanel(condition = "input.tab == 'challenge' && input.dataset == 'Training set'",
                            valueBoxOutput("leaderBox")
           )
         ),
         fluidRow(
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
                                label = "Submit", 
                                icon = icon("upload")), 
                   align = "center")
           ),
           conditionalPanel(condition = "input.tab == 'challenge'",
                            box(title = strong("Challenge"),
                                solidHeader = TRUE,
                                width = 6, 
                                status = "danger", 
                                height = 250,
                                radioButtons("dataset", "Dataset",
                                             choices = c("Training set", "Evaluation set"),
                                             selected = "Training set")
                            )
           )
           
         ),
         fluidRow(
           conditionalPanel(condition = "(input.tab == 'challenge') && (input.dataset == 'Training set')",
                            box(title = strong("Leaderboard"),
                                solidHeader = TRUE,
                                width = 12, 
                                status = "success",
                                materialSwitch(inputId = "switch", label = "Switch between SNP and regions", status = "primary", right = TRUE),
                                conditionalPanel(condition = "input.switch == false",
                                                 tabsetPanel(
                                                   source(file.path("ui", "tab-barchart.R"), local = TRUE)$value,
                                                   source(file.path("ui", "tab-table.R"), local = TRUE)$value
                                                 )
                                ),
                                conditionalPanel(condition = "input.switch == true",
                                                 tabsetPanel(
                                                   source(file.path("ui", "tab-barchart-region.R"), local = TRUE)$value,
                                                   source(file.path("ui", "tab-table-region.R"), local = TRUE)$value
                                                 )
                                )
                            )
           )
         )
)
