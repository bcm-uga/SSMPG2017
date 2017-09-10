### Packages
rm(list = ls())
library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(dbplyr)
library(DT)
library(RSQLite)
library(shinyWidgets)
library(stringr)
library(shinyjs)

db.file <- "db.sqlite3"

challenge_names <- c("Dahu", "Cichlid")

init_db <- function() {
  
  ## get username
  name <- readline(prompt = "User name: ")
  
  ## get mdp
  password <- readline(prompt = "User password: ")
  password.hash <- digest::digest(paste0("SSMPG2017", password), algo = "md5")
  
  ## write in table
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = db.file)
  dplyr::db_insert_into(db, "user", tibble::tibble(name = name, password = password.hash))
  dplyr::db_insert_into(con = db, 
                        table = "submission", 
                        values = tibble::tibble(name = name,
                                                date = as.character(Sys.time()),
                                                challenge = "1",
                                                dataset = "Training set",
                                                methods = "None",
                                                candidates = "0",
                                                regions = "0")
  )
  dplyr::db_insert_into(con = db, 
                        table = "submission", 
                        values = tibble::tibble(name = name,
                                                date = as.character(Sys.time()),
                                                challenge = "1",
                                                dataset = "Evaluation set",
                                                methods = "None",
                                                candidates = "0",
                                                regions = "0")
  )
  dplyr::db_insert_into(con = db, 
                        table = "submission", 
                        values = tibble::tibble(name = name,
                                                date = as.character(Sys.time()),
                                                challenge = "2",
                                                dataset = "Real",
                                                methods = "None",
                                                candidates = "0",
                                                regions = "None")
  )
  RSQLite::dbDisconnect(db)
  
}

add_submission <- function(user.name, 
                           password, 
                           challenge, 
                           dataset, 
                           methods,
                           candidates,
                           regions) {
  
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = db.file)
  
  ## add submission
  dplyr::db_insert_into(con = db, 
                        table = "submission", 
                        values = tibble::tibble(name = user.name,
                                                date = as.character(Sys.time()),
                                                challenge = challenge,
                                                dataset = dataset,
                                                methods = methods,
                                                candidates = candidates,
                                                regions = regions)
  )
  RSQLite::dbDisconnect(db)
}
