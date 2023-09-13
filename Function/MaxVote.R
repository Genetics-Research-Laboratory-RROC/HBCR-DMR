#================= Calculate max frequented entire in input vector =================
MaxVote <- function(VecLogical){
  return(names(which.max(table(VecLogical))))
}