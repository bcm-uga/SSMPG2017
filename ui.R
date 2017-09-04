header <- dashboardHeader(
  title = "SSMPG 2017"
)

body <- dashboardBody(
  fluidPage(
    tagList(
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
        tags$script(type = "text/javascript", src = "busy.js")
      )
      
    ),
    theme = shinythemes::shinytheme("cerulean"),
    tabItems(
      tabItem("challenge",
              fluidPage(
                theme = shinythemes::shinytheme("flatly"),
                title = " ",
                source(file.path("ui", "tab-challenge.R"), local = TRUE)$value
              )
      ),
      tabItem("vignette",
              fluidPage(
                theme = shinythemes::shinytheme("flatly"),
                title = " ",
                source(file.path("ui", "tab-vignette.R"), local = TRUE)$value  
              )
      )
    )
  )
)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Challenge", tabName = "challenge", icon = icon("trophy"), selected = TRUE),
    menuItem("Vignette", tabName = "vignette", icon = icon("book")), 
    menuItem("Github", icon = icon("github"), href = "https://github.com/bcm-uga/SSMPG2017")
  )  
)

dashboardPage(header, 
              sidebar, 
              body)
