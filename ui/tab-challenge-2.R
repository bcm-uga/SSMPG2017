tabPanel("Challenge 2",
         column(width = 4,
                fluidRow(
                  box(width = NULL, status = "warning",
                      fileInput(
                        "file_intrg_geno",
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
                  box(width = NULL, status = "primary",
                      radioButtons("dataset-2", "Dataset",
                                   choices = c("Training set", "Evaluation set"),
                                   selected = "Training set")
                  )
                ),
                fluidRow(
                  conditionalPanel(condition = "output.file_intrg_pop == true",
                                   box(width = NULL, status = "primary",
                                       selectInput("ancstrl1", "Ancestral 1", choices = c("2", "Tricho"), 
                                                   selected = NULL, 
                                                   multiple = FALSE,
                                                   selectize = TRUE,
                                                   width = NULL,
                                                   size = NULL),
                                       selectInput("ancstrl2", "Ancestral 2", choices = c("2", "Tricho"), 
                                                   selected = NULL, 
                                                   multiple = FALSE,
                                                   selectize = TRUE,
                                                   width = NULL,
                                                   size = NULL)
                                       
                                   )
                  )
                )
         ),
         mainPanel(width = 8
         )
)