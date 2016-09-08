plsda_stat <-
  function (scaling,silt) {
    pwd.x = paste(getwd(), "/Preprocessing_Data_", 
                  scaling, "/ProcessedTable.csv", sep="")
    x = read.csv(pwd.x, header=TRUE)
    x.x = x[,2:ncol(x)]
    rownames(x.x) = x[,1]
    pwdK = paste(getwd(), "/Preprocessing_Data_", scaling, "/class.csv", sep="")
    k = read.csv(pwdK, header=TRUE)
    k.x = matrix(k[,-1], ncol=1)
    x.n = cbind(k.x, x.x)
    sorted = x.n[order(x.n[,1]),]
    g = c()
    for (i in 1:nrow(sorted)) { 
      if (any(g == sorted[i,1])) {g=g} 
      else {g=matrix(c(g,sorted[i,1]), ncol=1)}
    }
    dimB = nrow(g)*nrow(sorted)
    B = matrix(rep(NA, dimB), ncol=nrow(g))
    for (i in 1:nrow(sorted)) {
      for (j in 1:nrow(g)) {
        jn <- g[j,1] 
        if (sorted[i,1] == jn) { 
          B[i,j] = 1}
        else {
          B[i,j] = 0
        }
      }
    }
    
    #requireNamespace(pls)
    #library(pls)
    sorted.x = sorted[,-1]
    sorted.un = matrix(unlist(sorted.x), ncol=ncol(sorted.x))
    
    P = pls::plsr(B ~ sorted.un, 
                  method = c("oscorespls"),ncomp = 2,validation = "CV")
    
    message( "\nPLS(-DA) Two Component Model Summary\n")
    #print(summary(P))
    rownames(P$scores) =  rownames(sorted.x)
    rownames(P$loadings) =  colnames(sorted.x)
    dirout = paste(getwd(), "/PLS_DA_", scaling, "/", sep="")
    dir.create(dirout)
    out.score = paste(dirout, "PLSDA_Scores_", scaling, ".csv", sep="")
    write.csv(P$scores, out.score)
    out.load = paste(dirout, "PLSDA_Loadings_", scaling, ".csv", sep="")
    write.csv(P$loadings, out.load)
    k = matrix(sorted[,1], ncol=1)
    tutticolors=matrix(c(1,2,3,4,5,6,7,8,"rosybrown4", 
                         "green4", "navy", "purple2", "orange", "pink", 
                         "chocolate2", "coral3", "khaki3","thistle",
                         "turquoise3","palegreen1","moccasin","olivedrab3",
                         "azure4","gold3","deeppink"), ncol=1)
    col=c()
    for(i in 1:nrow(k)) {
      col=c(col, tutticolors[k[i,],])
    }
    if (ncol(P$scores) == 1) {
      xlab = "Samples"
      ylab = "Score values Component 1"
      graphics::plot(P$scores[,1], col=col, pch=19, xlab=c(xlab), 
                     ylab=c(ylab), main = paste(
                       "PLS-DA Score Plot (", scaling, ")", sep=""))
      lim1 = nrow(P$scores)*2
      axis(1, at=c(-lim1,lim1), col="grey", pos=c(0,0), labels=FALSE, lwd=1)
      text(P$scores[,1], col=col, cex=0.5, pos=1, labels=rownames(P$scores))
      pwdout=paste(dirout, "ScorePlot_PLSDA_1Component_",
                   scaling, ".pdf", sep="")
      dev.copy2pdf(file=pwdout)
      Max.pc2 = 1.1*(max(P$loadings[,1]))
      Min.pc2 = 1.1*(min(P$loadings[,1]))
      Mpc2=c(Min.pc2,Max.pc2)
      pwdout1=paste(dirout, "W.cPlot_PLSDA_1Component_", 
                    scaling, ".pdf", sep="")
      pdf(pwdout1)
      #dev.new()
      graphics::plot(P$loadings[,1], ylim=Mpc2, main = paste(
        "PLS-DA Loading Plot (", scaling, ")", sep=""), 
        xlab="Variables", ylab="W*c values (Component1)")
      text(P$loadings[,1], cex=0.7, pos=1, labels=rownames(P$loadings))
      pwdout1=paste(dirout, "W.cPlot_PLSDA_1Component_", 
                    scaling, ".pdf", sep="")
      dev.off()
    } else {
      pairs = paste(dirout, "Pairs_PLSDA_", scaling, ".pdf", sep="")
      pdf(pairs)
      pairs = c()
      if (ncol(P$scores) >= 3) {pairs = c(3)} else {pairs = c(ncol(P$scores))}
      
      pairs(P$scores[,1:pairs],col=col)
      #pairs = paste(dirout, "Pairs_PLSDA_", scaling, ".pdf", sep="")
      dev.off()
    }
    p.v1 = matrix(P$Xvar, ncol=1)
    p.v.csv1 = paste(dirout, "PLSDA_P_", scaling, ".csv", sep="")
    write.csv(p.v1, file=p.v.csv1)
    p.vtot = matrix(P$Xtotvar, ncol=1)
    p.vtot.csv = paste(dirout, "PLSDA_Ptot_", scaling, ".csv", sep="")
    write.csv(p.vtot, file=p.vtot.csv, row.names = FALSE)
    dp.vtot = as.numeric(p.vtot)
    p.v = p.v1/dp.vtot
    p.v = matrix(100*p.v1/dp.vtot, ncol=1)
    colnames(p.v) = c("R2X")
    p.v.csv = paste(dirout, "PLSDA_R2X_", scaling, ".csv", sep="")
    write.csv(p.v, file=p.v.csv)
    #message(date(),"\\PLS-DA Finished!")
   ######################## VIP#########
    corr=cor(sorted.un,P$scores[,1])
    splot=cbind(P$loadings[,1],corr)
    splotstat = paste(dirout, "PLSDA_SPlot_", scaling, ".csv", sep="")
    write.table(splot,splotstat)
    pwdsplot=paste(dirout, "PLS_DA_SPlot_", scaling, ".pdf", sep="")
    pdf(pwdsplot)
    graphics::plot(splot[,1],splot[,2],  main = paste(
      "PLS-DA SPlot (", scaling, ")", sep=""), 
      xlab="W*c values (Component1)", ylab="W*c values (Component1)")
    text(splot, cex=0.7, pos=1, labels=colnames(sorted.x))
    dev.off()
    # significant comp, when Q2S>0 R2_Q2........................................
    #R2X_pls=explvar(P)
    #R2X_pls = data.frame(P$"Xvar")
    Q2_all=as.data.frame(pls::R2(P,"all",intercept = 0)$"val") 
    #R2[which estimate,which response,which comp]
    unk <- t(Q2_all)
    unk[apply(unk, 1, function(x) !all(is.na(x))),] -> Q2_all
    #rownames(Q2_all) = rownames(R2X_pls)
    #stat=cbind(R2X_pls*0.01,Q2_all)
    colnames(Q2_all)<-c("R2Y(cum)","Q2(cum)")
    R2T <- pls::R2(P,"all",intercept = 0)$"val"
    R2_Q2 <- data.frame(R2T[,,2]) # two components
    rownames(R2_Q2) <- c("R2Y(cum)","Q2Y(cum)")
    #cat("PLS(-DA)")
    cat(nrow(x), "samples x", ncol(x)-1, "variables\n\n")
    cat("Cumulative Proportion of Variance Explained: R2X(cum) = ", 
        p.v[1] + p.v[2],"%",sep="")
    cat("\n\nCumulative Proportion of Response(s):\n")
    print(R2_Q2)
    plsdastat = paste(dirout, "PLSDA_R2Q2_", scaling, ".csv", sep="")
    write.table(Q2_all,plsdastat)
    
    # PERMUTATION ........................................
    
    leng = length(g)
    if (leng == 2) {
    message("\nPermutation of PLSDA Model START...!")
    #message("\\...........Tea Time! Take A Rest!...........")
    #Y=c()
    Dx=list(X=sorted.un,Y=B[,1])
    #attach(D)
    #Y <- B[,1]
    ##
    #slt <- 2      ## 
    #silt = silt  ## permutation time
    #r2_sim_CV <- with(Dx,permut(Dx, silt = silt))
    
    r2_sim_CV <- with(Dx,permut(Dx, silt = silt))
    
    colnames(r2_sim_CV) <- c("R2","Q2","correlation.coefficient")
    RQ_line = c(Q2_all[3,],1)
    r2_sim_P <- rbind(r2_sim_CV,RQ_line)
    permutstat = paste(dirout, "PLSDA_permut_", scaling, ".csv", sep="")
    write.table(r2_sim_P,permutstat)
    permuplot=paste(dirout, "Permutation_", scaling, ".pdf", sep="")
    par(cex.axis=1,cex.lab=1)
    pdf(permuplot,width=6, height=6)
    lim_R = 1.3*(max(abs(r2_sim_P[,1])))
    lim_Q = 1.3*(max(abs(r2_sim_P[,2])))
    if (lim_R > lim_Q) {lim = c(-lim_R)} else {lim = c(-lim_Q)}
    if (Q2_all[3,1] > Q2_all[3,2]) {lim_mY = c(Q2_all[3,1])
    } else {lim_mY = c(Q2_all[3,2])}
    graphics::plot(abs(r2_sim_P[,3]), r2_sim_P[,1], 
                   col="grey", xlab = "Correlation", ylab = "", 
                   xlim = c(-0.05,1.05), ylim = c(lim,1.3*lim_mY), 
                   pch=19,  main = paste("Permutation Plot (", scaling, ")", 
                                         sep=""))
    points(abs(r2_sim_P[,3]), r2_sim_P[,2], col="black",
           xlim = c(-0.05,1.05), ylim = lim, pch=19)
    graphics::abline(h=Q2_all[3,1],col="grey",lwd = 1.2)
    graphics::abline(h=Q2_all[3,2],col="black",lwd = 1.2)
    points(1,Q2_all[3,1],col="grey",lwd=1,pch=19)
    points(1,Q2_all[3,2],col="black",lwd=1,pch=19)
    #interceptR <- as.numeric(interR)
    #interQ <- lm(median(r2_sim_CV[,2])~median(r2_sim_CV[,3]))$fitted.values
    #interceptQ <- as.numeric(interQ)
    #segments(0,interceptR,1,Q2_all[3,1],col="#D55E00")
    #segments(0,interceptQ,1,Q2_all[3,2],col="#009E73")
    #abline(h=Q2_all[3,1],col="grey")
    #abline(h=Q2_all[3,1],col="gery")
    legend("bottomright",legend = c("R2Y","Q2Y"),bty="n",
           cex=1.2,pch = 19, col = c("grey","black"))
    #legend("bottomright",1.0, legend = c("Q2"),bty="n",
    #cex=1.5, pch = 19, col = "#009E73")
    grid(lwd = 0.8)
    dev.off()
    #message( "\nR2 and Q2!")
    #output <- Q2_all[1:6,]
    #colnames(output) <- c("R2 comp1", "R2 comp2")
    #print(output)
    #message(date(),"\\permutation time Finished!")
    ############################
    #permutstat = paste(dirout, "PLSDA_permut_", scaling, ".csv", sep="")
    #write.table(r2_sim_CV,permutstat)
    #permuplot=paste(dirout, "Permutation_", scaling, ".pdf", sep="")
      #pdf(permuplot)
      #x_lim_CV=min(r2_sim[,2]) - 1.2*sd(r2_sim[,2])
      #minV = min(r2_sim_CV[,2]) - 1.2*sd(r2_sim_CV[,2])
      #maxV = max(r2_sim_CV[,2]) + 5*sd(r2_sim_CV[,2])
      #maxV = max(r2_sim_CV[,2]) + abs(min(r2_sim_CV[,2]))
      
      #if (max(r2_sim_CV[,2]) > 1.2*Q2_all[3,2]) {
        #x_lim_CV = c(minV, maxV)
      #} else {
        #x_lim_CV = c(minV, 1.2*Q2_all[3,2])
      #}
      
      #hist(r2_sim_CV[,2],breaks=40,col='darkblue',border='white',
           #xlim=x_lim_CV,ylim=c(0,85),main="Permutation",axes = TRUE,
           #xlab= "Q2")
      #segments(Q2_all[3,2],38,Q2_all[3,2],0,col="#FF4F00",lwd=4)
      #points(Q2_all[3,2],42,col="red",lwd=1,bg="#FF4F00",pch=8)
      #d = density(r2_sim_CV[,2])
      #polygon(d,col = "lightgreen",border = NA)
      #dev.off()
      #message(date(),"\\Q2!")
      #output <- Q2_all[1:6,]
      #colnames(output) <- c("R2 comp1", "R2 comp2")
      #print(output)
      #message(date(),"\\permutation time Finished!")
  } 
    if (leng > 2) { 
        message( "\nWarning: More than two groups, Permutation Test Free!")
    }
    #detach()
    # PERMUTATION ..................................................
 #VIP
if (nrow(P$Yloadings) > 2) {
  message( "\nWarning: VIP was only implemented for single-response models!")
} else {
  if (P$method != "oscorespls")
    stop("\nOnly implemented for orthogonal scores algorithm.  
         Refit with 'method = \"oscorespls\"'")
  
  SS <- c(P$Yloadings[1,])^2 * colSums(P$scores^2)
  Wnorm2 <- colSums(P$loading.weights^2)
  SSW <- sweep(P$loading.weights^2, 2, SS / Wnorm2, "*")
  vip <- sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
  colnames(vip) <- colnames(x.x)
  vipstat = paste(dirout, "PLSDA_VIP_", scaling, ".csv", sep="")
  write.table(vip,vipstat)
  }
}


