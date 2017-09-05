server <- function(input, output) {
  
  source(file.path("server", "users.R"), local = TRUE)$value
  source(file.path("server", "submission.R"), local = TRUE)$value
  source(file.path("server", "tab-barchart.R"), local = TRUE)$value
  source(file.path("server", "tab-table.R"), local = TRUE)$value
  
}
