tabPanel(
  strong("Summary - Region"), 
  DT::dataTableOutput("summary_region"),
  actionButton(inputId = "summary_region",
               label = "Display/Update"),
  value = 4
)