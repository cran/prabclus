"incmatrix" <-
function(regmat){
  nart <- ncol(regmat)
  nreg <- nrow(regmat)
  neq <- 0
  incmat <- diag(nart)
  for (i in 1:(nart-1))
    for (j in (i+1):nart){
#      cat (i," ",j," ",sum(regmat[,i]<regmat[,j]),sum(regmat[,j]<regmat[,i]),"\n")
      if (sum(regmat[,i]<regmat[,j])==0)
        incmat[i,j] <- 1
      if (sum(regmat[,j]<regmat[,i])==0)
        incmat[j,i] <- 1
      if (identical(regmat[,i],regmat[,j]))
#        incmat[i,j] <- 0
        neq <- neq+1
    }
  out <- list(m=incmat, ninc=sum(incmat)-2*neq-nart, neq=neq)
  out
}
