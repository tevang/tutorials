# You could also have several functions in a single R file and still document them separately. Simply put an identifier 
# starting with ## ---- before each function definition and then create empty chunks referring to each one of the identifiers.

library("ROCR")
library("hash")

## ---- read_scores
read_scores <- function(RESULTS_FILE, ACTIVITIES_FILE) {
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
  scores <- rep(0, nrow(x))
  labels <- rep("", nrow(x))
  i = 1
  for (molname in names(score_dict)) {
    scores[i] <- score_dict[[molname]]
    labels[i] <- label_dict[[molname]]
    i <- i+1
  }
  pred = prediction(-1*scores, labels)  ;# IMPORTANT: -1* because the function prediction() assumes that the highest the score the better
  return(pred)
}