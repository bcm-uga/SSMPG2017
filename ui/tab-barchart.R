tabPanel(title = strong("Bar chart"),
         plotlyOutput("barchart"),
         actionButton(inputId = "display",
                      label = "Display/Update"), 
         value = 1
)