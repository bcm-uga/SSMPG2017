db.file <- "db.sqlite3"

## create db

## create user
db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = db.file)
res <- RSQLite::dbSendQuery(conn = db,
                            "CREATE TABLE user(
                            name CHARACTER PRIMARY KEY,
                            password CHARACTER)"
                            )
RSQLite::dbClearResult(res)
RSQLite::dbDisconnect(db)

## create submission
db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = db.file)
res <- RSQLite::dbSendQuery(conn = db,
                            "CREATE TABLE submission(
                            name CHARACTER,
                            date CHARACTER,
                            challenge CHARACTER,
                            dataset CHARACTER,
                            methods CHARACTER,
                            candidates CHARACTER,
                            FOREIGN KEY (name) REFERENCES user(name))"
                            )
RSQLite::dbClearResult(res)
RSQLite::dbDisconnect(db)