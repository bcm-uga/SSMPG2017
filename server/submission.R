userModal <- function(failed = 0) {
  modalDialog(
    span('You are about to submit for challenge ',
         input$challenge, 
         'with respect to the',
         paste0(input$dataset, '.'),
         'If this is correct, please enter your team name and your password.'),
    
    textInput("username", "Team"),
    textInput("password", "Password"),
    
    if (failed == 1)
      div(tags$b("User not listed in the database", style = "color: red;")),
    
    if (failed == 2)
      div(tags$b("Wrong password", style = "color: red;")),
    
    if (failed == 3)
      div(tags$b("File not found", style = "color: red;")),
    
    footer = tagList(
      modalButton("Cancel"),
      actionButton("ok", "OK")
    )
  )
}

# Show modal when button is clicked.
observeEvent(input$submit, {
  showModal(userModal())
})

observeEvent(input$ok, {
  # Check that data object exists and is data frame.
  hash.pwd <- digest::digest(paste0("SSMPG2017", input$password), algo = "md5")
  if (!(input$username %in% user.df$name)) {
    showModal(userModal(failed = 1))  
  } else if (input$username %in% user.df$name && hash.pwd != user.df$password[user.df$name == input$username]){
    showModal(userModal(failed = 2))  
  } else if (is.null(file <- input$subm$datapath)){
    showModal(userModal(failed = 3))  
  } else {
    removeModal()
  }
  
  add_submission(user.name = input$username,
                 password = input$password,
                 challenge = as.character(input$challenge),
                 dataset = input$dataset,
                 fdr = 0,
                 power = 0,
                 score = 0)
})


