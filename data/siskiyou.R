siskiyou <- as.matrix(read.table("LeiMik1.dat"))
siskiyou.nb <- list()
for (i in 1:6)
  siskiyou.nb <- c(siskiyou.nb,list(scan(file="LeiMik1NB.dat",
                   skip=i-1,nlines=1, quiet=TRUE)))
remove(i)
siskiyou.groups <- scan("LeiMik1G.dat")
