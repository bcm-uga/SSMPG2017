createModal <- function(failed = 0) {
  modalDialog(
    span("Enter your team name"),
    textInput("new_username", label = NULL, width = NULL),
    span("Enter your password"),
    passwordInput("new_password", label = NULL, width = NULL),
    
    if (failed == 1)
      div(tags$b("Username already exists", style = "color: red;")),
    
    footer = tagList(
      modalButton("Cancel"),
      actionButton("create", "Create")
    )
  )
}

# Show modal when button is clicked.
observeEvent(input$users, {
  showModal(createModal())
})

observeEvent(input$create, {
  hash.pwd <- digest::digest(paste0("SSMPG2017", input$new_password), algo = "md5")
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = "db.sqlite3")
  user.df <- RSQLite::dbReadTable(db, "user")
  if (input$new_username %in% user.df$name) {
    showModal(createModal(failed = 1))  
  } else {
    dplyr::db_insert_into(db, "user", tibble::tibble(name = input$new_username, password = hash.pwd))
    removeModal()
    showModal(modalDialog(
        span(paste(input$new_username, "has been successfully created.")),
        easyClose = TRUE,
        footer = NULL
        )
    )
  }
  RSQLite::dbDisconnect(db)
})



