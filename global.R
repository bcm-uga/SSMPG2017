### Packages

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(dbplyr)
library(DT)
library(RSQLite)
library(shinyWidgets)

all.gt <- readRDS("anssnp.rds")
all.gt.reg <- readRDS("ansreg.rds")

add_submission <- function(user.name, 
                           password, 
                           challenge, 
                           dataset, 
                           fdr, 
                           power, 
                           score, 
                           methods,
                           candidates,
                           evaluation) {
  
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = "db.sqlite3")
  
  ## add submission
  dplyr::db_insert_into(con = db, 
                        table = "submission", 
                        values = tibble::tibble(name = user.name,
                                                date = as.character(Sys.time()),
                                                challenge = challenge,
                                                dataset = dataset,
                                                methods = methods,
                                                candidates = candidates)
  )
  RSQLite::dbDisconnect(db)
}
