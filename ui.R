header <- dashboardHeader(
  title = "SSMPG 2017"
)

body <- dashboardBody(
  uiOutput("ui")
)

sidebar <- dashboardSidebar(
  sidebarMenu(id = "tab",
              menuItem("Challenge", icon = icon("trophy"),
                       menuSubItem("Dahu", tabName = "challenge"),
                       menuSubItem("Cichlid", tabName = "challenge_2")), 
              menuItem("Vignette", tabName = "vignette", icon = icon("book")), 
              menuItem("GitHub", icon = icon("github"), href = "https://github.com/bcm-uga/SSMPG2017"),
              menuItem("Datasets", tabName = "download", icon = icon("download")),
              div(actionButton(inputId = "users", 
                               label = "Create team", 
                               icon = icon("users")), 
                  align = "center")
  )  
)

dashboardPage(header, 
              sidebar, 
              body)
