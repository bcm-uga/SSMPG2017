### Packages

library(shiny)
library(googlesheets)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(plotly)
library(RSQLite)


parse_file <- function(file) scan(file, what = "")  

challenges <- paste("Challenge", 1:2)
OPTIONS <- c("Local FDR", "Covariables", "Latent factors")

df <- data.frame(team = rep(c("bcm", "white walkers", "power rangers", "teletubbies"), 2),
                 score = c(0.8, 0.4, 0.7, 0.9, 0.2, 0.5, 0.1, 0.9),
                 challenge = c(1, 1, 1, 1, 2, 2, 2, 2))