gam <- function(N, A, numConvs, a, b){
  supp <- seq(from=0, to = A*2, length.out = N*2)
  
  f <- c(0,diff(pgamma(supp,a,b)))
  
  ft <- fft(f)
  
  newf <- Re(fft(ft^numConvs, inverse = TRUE))/(N*2)
  
  real <- pgamma(supp[1:N], a*numConvs,b)
  
  low <- cumsum(newf[1:N])-cumsum(newf[(N+1):(2*N)])
  upp <- cumsum(newf[1:(N+2)])[-(1:2)]
  
  f <- function(i) {
    sum(abs(low*(1-i)+upp*i-real))
  }
  
  min.error <- optimize(f, interval = c(0,1))$minimum
  
  return(list(cbind(low, real, upp), min.error))
}


gam(N=2^3, A=2^6, numConvs = 3, a = 2, b = 3)[[2]]


best.est2 <- matrix(NA, nrow = 8, ncol = 8)

for(i in 3:10){
  for(j in 3:10){
    best.est2[i-2,j-2] <- gam(N=2^i, A=2^j, numConvs = 3, a = 2, b = 3)[[2]]
  }
}

best.est2

