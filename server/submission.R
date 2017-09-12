userModal <- function(failed = 0) {
  msg <- `if`(n.challenge() == 1, 
              paste0("You are about to submit for the ",
              challenge_names[n.challenge()],
              " challenge (",
              input$dataset,
              "). If this is correct, please enter your team name and your password to validate your submission."),
              paste0("You are about to submit for the ",
              challenge_names[n.challenge()],
              " challenge. If this is correct, please enter your team name and your password to validate your submission."))
              
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
                                width = "100%")),
    
    span(msg),
    
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

observeEvent(input$ok, {
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = db.file)
  user.df <- RSQLite::dbReadTable(db, "user")
  RSQLite::dbDisconnect(db)
  hash.pwd <- digest::digest(paste0("SSMPG2017", input$password), algo = "md5")
  if (!(input$username %in% user.df$name)) {
    showModal(userModal(failed = 1))  
  } else if (input$username %in% user.df$name && hash.pwd != user.df$password[user.df$name == input$username]) {
    showModal(userModal(failed = 2))  
  } else if (is.null(submission$x)) {
    showModal(userModal(failed = 3))  
  } else if (length(input$methods) == 0) {
    showModal(userModal(failed = 4))
  } else {
    removeModal()
    
    if (input$tab == "challenge" && input$dataset == "Training set") {
      fb <- all.gt.reg[[1]]$region.start
      fe <- all.gt.reg[[1]]$region.end
      reg <- paste(unique(sapply(submission$x, 
                                 FUN = function(X) {
                                   which(fb <= X & fe >= X)
                                 })), collapse = ", ")
    } else if (input$tab == "challenge" && input$dataset == "Evaluation set") {
      fb <- all.gt.reg[[2]]$region.start
      fe <- all.gt.reg[[2]]$region.end
      reg <- paste(unique(sapply(submission$x, 
                                 FUN = function(X) {
                                   which(fb <= X & fe >= X)
                                 })), collapse = ", ")
    } else if (input$tab == "challenge_2") {
      reg <- "None"
    }
    
    add_submission(user.name = input$username,
                   password = input$password,
                   challenge = as.character(n.challenge()),
                   dataset = `if`(input$tab == "challenge", 
                                  input$dataset,
                                  "Real"),
                   methods = paste(input$methods, collapse = "; "),
                   candidates = paste(submission$x, collapse = ", "),
                   regions = reg)
    
    closing_msg <- `if`(n.challenge() == 1, 
                        paste0("Your entry for the ",
                               challenge_names[n.challenge()],
                               " challenge (",
                               input$dataset,
                              ") has been successfully submitted."),
                        paste0("Your entry for the ",
                               challenge_names[n.challenge()],
                               " challenge has been successfully submitted."))
    
    showModal(modalDialog(
      span(closing_msg),
      easyClose = TRUE,
      footer = NULL)
    )
  
    # Reset fileInput
    submission$x <- NULL
    reset("subm") 
    reset("subm_2") 

  }
})


