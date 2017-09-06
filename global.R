### Packages

library(shiny)
library(googlesheets)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(dbplyr)
library(plotly)
library(DT)
library(RSQLite)

add_user <- function() {
  
  ## get username
  name <- readline(prompt="User name: ")
  
  ## get mdp
  password <- readline(prompt="User password: ")
  password.hash <- digest::digest(paste0("SSMPG2017", password), algo = "md5")
  
  ## write in table
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = "db.sqlite3")
  dplyr::db_insert_into(db, "user", tibble::tibble(name = name, password = password.hash))
  RSQLite::dbDisconnect(db)
}

add_submission <- function(user.name, 
                           password, 
                           challenge, 
                           dataset, 
                           fdr, 
                           power, 
                           score, 
                           methods,
                           candidates) {
  
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = "db.sqlite3")
  
  ## add submission
  dplyr::db_insert_into(db, "submission", tibble::tibble(name = user.name,
                                                         date = as.character(Sys.time()),
                                                         challenge = challenge,
                                                         dataset = dataset,
                                                         fdr = fdr,
                                                         power = power,
                                                         score = score,
                                                         methods = methods,
                                                         candidates = candidates)
                        )
  RSQLite::dbDisconnect(db)
}
