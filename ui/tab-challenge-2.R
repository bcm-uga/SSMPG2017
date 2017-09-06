tabPanel("Challenge",
         fluidRow(
           box(title = strong("Submission"),
               solidHeader = TRUE,
               width = 6, 
               status = "warning", 
               height = 250,
               fileInput(
                 "subm_2",
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
               div(actionButton(inputId = "submit_2",
                                label = "Submit"), 
                   align = "center")
           )
           
         )
)
