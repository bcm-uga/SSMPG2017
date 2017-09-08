### Packages

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(dbplyr)
library(DT)
library(RSQLite)
library(shinyWidgets)
library(stringr)

all.gt <- readRDS("anssnp.rds")
all.gt.reg <- NULL
all.gt.reg[[1]] <- readRDS("ansreg1a.rds")
all.gt.reg[[2]] <- readRDS("ansreg2a.rds")

add_submission <- function(user.name, 
                           password, 
                           challenge, 
                           dataset, 
                           methods,
                           candidates,
                           regions) {
  
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = "db.sqlite3")
  
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
