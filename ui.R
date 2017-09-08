header <- dashboardHeader(
  title = "SSMPG 2017"
)

body <- dashboardBody(
  shinyjs::useShinyjs(),
  tagList(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
      tags$script(type="text/javascript", src = "busy.js")
    )
    
  ), 
  div(class = "busy",
      p("Updating..."),
      img(src = "loading.gif")
  ),
  uiOutput("ui")
)

sidebar <- dashboardSidebar(
  sidebarMenu(id = "tab",
              menuItem("Challenge", icon = icon("trophy"),
                       menuSubItem("Dahu", tabName = "challenge"),
                       menuSubItem("Cichlid", tabName = "challenge_2")), 
              menuItem("Vignette", tabName = "vignette", icon = icon("book")), 
              menuItem("GitHub", icon = icon("github"), href = "https://github.com/bcm-uga/SSMPG2017"),
              #menuItem("Datasets", tabName = "download", icon = icon("download")),
              div(actionButton(inputId = "users", 
                               label = "Create team", 
                               icon = icon("users")), 
                  align = "center")
  )  
)

dashboardPage(header, 
              sidebar, 
              body)
