tabPanel(
  strong("Summary"), 
  DT::dataTableOutput("summary"),
  actionButton(inputId = "summary",
               label = "Display/Update"),
  value = 2
)