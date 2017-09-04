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
      tabItem("challenge-1",
              fluidPage(
                theme = shinythemes::shinytheme("flatly"),
                title = " ",
                source(file.path("ui", "tab-challenge-1.R"), local = TRUE)$value
              )
      ),
      tabItem("challenge-2",
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
      )
    )
  )
)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Challenge 1", tabName = "challenge-1", icon = icon("thumbs-o-up")),
    menuItem("Challenge 2", tabName = "challenge-2", icon = icon("hand-peace-o")),
    menuItem("Vignette", tabName = "vignette", icon = icon("book")), 
    menuItem("Github", icon = icon("github"), href = "https://github.com/BioShock38/pcadapt")
  )  
)

dashboardPage(header, 
              sidebar, 
              body)
