nb <- list()
for (i in 1:34)
  nb <- c(nb,list(scan(file="nb.dat",
                   skip=i-1,nlines=1, quiet=TRUE)))
remove(i)
