tabPanel("Challenge",
         fluidRow(
           valueBoxOutput("challengeBox"),
           conditionalPanel(condition = "input.tab == 'challenge' && input.dataset == 'Training set'",
                            valueBoxOutput("contestantBox"),
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
                                fluidRow(
                                  column(width = 4,
                                         conditionalPanel(
                                           condition = "input.switch == true",
                                           h4(strong("Evaluation per region"))
                                         ),
                                         conditionalPanel(
                                           condition = "input.switch == false",
                                           h4(strong("Evaluation per SNP"))
                                         )
                                  ),
                                  column(width = 8,
                                         style = "margin-top: 12px;",
                                         materialSwitch(inputId = "switch", 
                                                        label = " ", 
                                                        status = "primary", 
                                                        right = TRUE)
                                         
                                  )
 
                                ),
                                
                                # Leaderboard panel
                                tabsetPanel(
                                  tabPanel(title = strong("Bar chart"),
                                           plotOutput("barchart")
                                  ),
                                  tabPanel(
                                    title = strong("Summary"), 
                                    DT::dataTableOutput("summary")
                                  )
                                ),
                                
                                # Update button
                                actionButton(inputId = "update",
                                             label = "Display/Update")
                            )
           )
         )
)
