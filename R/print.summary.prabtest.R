"print.summary.prabtest" <-
function(x, ...){
  cat("* Monte Carlo test for presence-absence data *\n\n")
  cat("Test statistics: ",x$teststat,", Tuning constant=",x$tuning,"\n")
  cat("Distance: ",x$distance,"\n")
  cat("Simulation runs: ",x$times,"\n")
  cat("Disjunction parameter: ",x$pd,"\n")
  if (x$pdfnb)
    cat("Neighbor-based correction of region probabilities was used.\n")
  else 
    cat("Neighbor-based correction of region probabilities was not used.\n")
  cat("Statistics value for original data: ",x$datac,"\n")
  cat("Mean for null data: ",x$rmean,", range: ",x$rrange,"\n")
  cat("p= ",x$p.value,"\n")
}
