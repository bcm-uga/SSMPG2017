### Packages

library(shiny)
library(googlesheets)
library(shinydashboard)


parse_file <- function(file) scan(file, what = "")  

challenges <- paste("Challenge", 1:2)
OPTIONS <- c("Local FDR", "Covariables", "Latent factors")