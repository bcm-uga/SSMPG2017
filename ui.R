header <- dashboardHeader(
  title = strong("SSMPG 2017")
)

body <- dashboardBody(
  fluidPage(
    theme = shinythemes::shinytheme("cerulean"),
    tabItems(
      tabItem("challenge",
              fluidPage(
                theme = shinythemes::shinytheme("flatly"),
                title = " ",
                source(file.path("ui", "tab-challenge.R"), local = TRUE)$value
              )
      ),
      tabItem("challenge_2",
              fluidPage(
                theme = shinythemes::shinytheme("flatly"),
                title = " ",
                source(file.path("ui", "tab-challenge-2.R"), local = TRUE)$value
              )
      ),
      tabItem("vignette",
              fluidPage(
                theme = shinythemes::shinytheme("flatly"),
                title = " ",
                source(file.path("ui", "tab-vignette.R"), local = TRUE)$value  
              )
      ),
      tabItem("download",
              fluidPage(
                theme = shinythemes::shinytheme("flatly"),
                title = " ",
                source(file.path("ui", "tab-download.R"), local = TRUE)$value  
              )
      )
    )
  )
)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Challenge 1", tabName = "challenge", icon = icon("trophy"), selected = TRUE),
    menuItem("Challenge 2", tabName = "challenge_2", icon = icon("trophy")),
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
