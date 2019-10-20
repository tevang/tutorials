# You could also have several functions in a single R file and still document them separately. Simply put an identifier 
# starting with ## ---- before each function definition and then create empty chunks referring to each one of the identifiers.

library("ROCR")
library("hash")


## ---- common_molnames
common_molnames <- function(RESULTS_FILES) {
  "
  A function to find the molnames that are common in all score files in order to compare the scoring 
  functions properly.
  "
  x = read.table(RESULTS_FILES[1], header = TRUE)
  colnames(x)[1] = "molname"  ; # add a column header and operate on x$molname, otherwise 'inersect' fails!
  valid_molnames <- unique(sort(x$molname)) ; # unique molnames
  for (i in 2:length(RESULTS_FILES)) {
    x = read.table(RESULTS_FILES[i], header = TRUE)
    colnames(x)[1] = "molname"
    valid_molnames <- intersect(valid_molnames, x$molname)
  }
  return(valid_molnames)
}

## ---- count_actives_inactives
count_actives_inactives <- function(ACTIVITIES_FILE, valid_molnames) {
  a = read.table(ACTIVITIES_FILE)
  colnames(a)[1] = "molname"
  colnames(a)[2] = "label"
  actives <- a$molname[a$label==1]
  inactives <- a$molname[a$label==0]
  active_num <- length(actives[actives %in% valid_molnames])
  inactive_num <- length(inactives[actives %in% valid_molnames])
  paste("The molecules that have been scored by all scoring functions consist of", active_num, "actives", 
        " and", inactive_num, "inactives.")
  num <- list(actives=active_num, inactives=inactive_num)
  return(num)
}

## ---- read_scores
read_scores <- function(RESULTS_FILE, ACTIVITIES_FILE, valid_molnames) {
"
  The valid_molnames list ensures that only the molnames that were scored by all scoring functions will
  be considered.
"
  x = read.table(RESULTS_FILE, header = TRUE)
  colnames(x)[2] = "score"
  # ignore the other columns
  score_dict = hash()
  for (i in seq(1, nrow(x))) { score_dict[x[i,1]] <- x[i,2] }
  
  a = read.table(ACTIVITIES_FILE)
  colnames(a)[1] = "molname"
  colnames(a)[2] = "label"
  label_dict = hash()
  for (i in seq(1, nrow(a))) { label_dict[a[i,1]] <- a[i,2] }
  scores <- rep(0, length(valid_molnames))
  labels <- rep("", length(valid_molnames))  ; # initialize labels to a list of size(scores) but without contents
  i = 1
  for (molname in valid_molnames) {
    scores[i] <- score_dict[[molname]]
    labels[i] <- label_dict[[molname]]
    i <- i+1
  }
  pred = prediction(-1*scores, labels)  ;# IMPORTANT: -1* because the function prediction() assumes that the highest the score the better
  return(pred)
}
