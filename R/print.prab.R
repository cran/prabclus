"print.prab" <-
function(x, ...){
  cat("Presence-absence matrix object with\n")
  cat(x$n.species," species and ",x$n.regions," regions,\n")
  cat("including regions neighborhoods and \n")
  cat("between-species distance matrix of type ",x$distance,".\n")
}
