create_fun <- function(w, WAIT = 0) {
  private <- function(i) w[[i]]
  obscure <- function(number_challenge) {
    true <- private(number_challenge)
    function(ens) {
      r <- mean(true %in% ens) # recall
      p <- mean(ens %in% true) # precision
      f <- 2 * p * r / (p + r) # f-score       
      data.frame(power = r, fdr = 1 - p, F_score = `if`(is.nan(f), 0, f)) 
    }
  }
  scoring <- function(ens, number_challenge) {
    obscure(number_challenge)(ens)
  }
  
  compiler::cmpfun(function(submission, opts, number_challenge, key) {
    
    now <- lubridate::now
    register <- function(x) googlesheets::gs_key(x, lookup = FALSE, 
                                                 visibility = "private")
    add_row <- googlesheets::gs_add_row
    
    minutes_from_now <- function(time) 
      as.double(lubridate::as.duration(now() - time - diff_time)) / 60
    
    my_doc <- register(key)
    diff_time <- now() - my_doc$reg_date
    
    wait <- WAIT - minutes_from_now(my_doc$updated)
    if (wait > 0) {
      stop(sprintf("You still have to wait %.1f minutes.", wait), call. = FALSE)
    } 
    
    add_row(my_doc, ws = number_challenge,
            input = c(as.character(now()), toString(opts), toString(submission)))
    
    score <- scoring(submission, number_challenge)
    
    results <- register("1YtuLCYNugnxdeDiTabIBhoX_C_aam62nokrPRoGkUr0")
    add_row(results, ws = number_challenge, 
            input = c(my_doc$sheet_title, as.character(now()), 
                      score[["power"]], score[["fdr"]],
                      score[["F_score"]], toString(opts)))
    
    score
  })
}

library(dplyr)

update_leaderboard <- compiler::cmpfun(function(number_challenge) {
  
  results <- googlesheets::gs_key(
    "1YtuLCYNugnxdeDiTabIBhoX_C_aam62nokrPRoGkUr0", 
    lookup = FALSE, visibility = "private")
  
  googlesheets::gs_read(results, ws = number_challenge) %>%
    group_by(team) %>%  
    summarise(rank = which.max(date),
              Power = power[rank],
              FDR = fdr[rank],
              F_score = score[rank],
              Date = date[rank],
              N = n()) %>%
    arrange(desc(F_score), N, Date) %>%
    mutate(Date = as.character(Date), rank = NULL) %>%
    rename_("Team" = "team")
})