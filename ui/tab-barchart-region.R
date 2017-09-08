tabPanel(title = strong("Bar chart - Region"),
         plotOutput("barchart_region"),
         actionButton(inputId = "display_region",
                      label = "Display/Update"), 
         value = 3
)