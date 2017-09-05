userModal <- function(failed = 0) {
  modalDialog(
    tags$div(align = "left", 
             class = "multicol",
             checkboxGroupInput(inputId = "methods", 
                                label = "Select the method(s) you have used for this submission:",  
                                choices = c("BAYPASS" = "BAYPASS",
                                            "HapFLK" = "HapFLK",
                                            "LEA" = "LEA",
                                            "OutFLANK" = "OutFLANK",
                                            "pcadapt" = "pcadapt",
                                            "REHH" = "REHH",
                                            "Selestim" = "SelEstim",
                                            "SweeD" = "SweeD",
                                            "Other" = "Other"),
                                width = NULL)),
    
    span('You are about to submit for challenge ',
         input$challenge, 
         'with respect to the',
         paste0(input$dataset, '.'),
         'If this is correct, please enter your team name and your password to validate your submission.'),
    
    textInput("username", "Team", width = NULL),
    passwordInput("password", "Password", width = NULL),
    
    if (failed == 1)
      div(tags$b("User not listed in the database", style = "color: red;")),
    
    if (failed == 2)
      div(tags$b("Wrong password", style = "color: red;")),
    
    if (failed == 3)
      div(tags$b("File not found", style = "color: red;")),
    
    if (failed == 4)
      div(tags$b("Check at least one box", style = "color: red;")),
    
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


# Show modal when button is clicked.
observeEvent(input$gotit, {
  removeModal()
})


observeEvent(input$ok, {
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = "db.sqlite3")
  user.df <- RSQLite::dbReadTable(db, "user")
  RSQLite::dbDisconnect(db)
  hash.pwd <- digest::digest(paste0("SSMPG2017", input$password), algo = "md5")
  if (!(input$username %in% user.df$name)) {
    showModal(userModal(failed = 1))  
  } else if (input$username %in% user.df$name && hash.pwd != user.df$password[user.df$name == input$username]) {
    showModal(userModal(failed = 2))  
  } else if (is.null(input$subm$datapath)) {
    showModal(userModal(failed = 3))  
  } else if (length(input$methods) == 0) {
    showModal(userModal(failed = 4))
  } else {
    removeModal()
    
    m <- paste(input$methods, collapse = "; ")
      
    submission <- scan(input$subm$datapath)
    all.gt <- readRDS("fake.rds")
    gt <- all.gt[[as.integer(input$challenge)]]
    r <- mean(gt %in% submission)
    p <- mean(submission %in% gt)
    g <- sqrt(r * p)
    add_submission(user.name = input$username,
                   password = input$password,
                   challenge = as.character(input$challenge),
                   dataset = input$dataset,
                   fdr = 1 - p,
                   power = r,
                   score = `if`(is.nan(g), 0, g),
                   methods = m)
    showModal(modalDialog(
      span(paste("Your entry for challenge", 
                 input$challenge, 
                 "has been successfully submitted.")),
      easyClose = TRUE,
      footer = NULL
    )
    )
  }
})


