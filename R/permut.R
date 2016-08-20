permut <-
  function(file, silt){
  r2_sim_CV_2 <- matrix(0,nrow = as.numeric(silt), ncol = 2)
  Y = file$Y
  corrY <- c()
  for(i in 1:silt){
    sid <- sample(length(Y),length(Y))
    p1=pls::plsr(Y[sid]~X,data = file, validation = "CV",method ="oscorespls")
    r2=pls::R2(p1,estimate = "CV")
    corrY[i] <- cor(Y[sid],Y)
    #r2_cum <- cum(r2)
    #r2_cum <- r2_cum[1:slt]
    r2_sim_CV_2[i,] <- r2$val[1,1,2:3]
  }
  r2_sim_CV <- data.frame(r2_sim_CV_2,corrY)
}