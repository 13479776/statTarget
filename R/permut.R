permut <-
  function(file, silt){
    Y = file$Y
    file = file
    env <- environment()
    r2_sim_CV_2 <- matrix(0,nrow = as.numeric(silt), ncol = 2)
    corrY <- c()
    #pb <- txtProgressBar(min = 1, max = silt, style = 3)
    pCV <- function(x){
      Y <-  get("Y", envir=env) 
      file <- get("file", envir=env) 
      sid <- sample(length(Y),length(Y))
      p1=pls::plsr(Y[sid]~X,data = file, ncomp = 2, validation = "CV",
                   method ="oscorespls")
      q2r2=pls::R2(p1,estimate = "all")
      corrY <- cor(Y[sid],Y)
      #r2_cum <- cum(r2)
      #r2_cum <- r2_cum[1:slt]
      r2_sim_CV_2 <- cbind(t(as.data.frame(q2r2$val[5:6])),corrY)
    }
    output <- plyr::llply(1:silt,.fun= pCV,.progress = "text")
    output = do.call(rbind,output)
  }    
  