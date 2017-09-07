tabPanel(title = strong("Bar chart - Region"),
         plotOutput("barchart_region"),
         actionButton(inputId = "display",
                      label = "Display/Update"), 
         value = 3
)