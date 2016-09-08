volcano <-
function(file, upper.lim, lower.lim, sig.lim) {
 pwdfile=paste(getwd(), "/Univariate/DataTable.csv", sep="")
 file=pwdfile
 x <- read.csv(file, sep=",", header=TRUE)
 x.x = x[,3:ncol(x)]
 rownames(x.x) = x[,2]
 k = matrix(x[,1], ncol=1)
 slink = paste(getwd(), "/PreTable","/slink.csv", sep="")
 slink = read.csv(slink, header=TRUE)
 x.n = cbind(k, x.x)
 sorted = x.n[order(x.n[,1]),]
 sorted.x = as.matrix(sorted[,-1], ncol=ncol(sorted)-1)
 g = c()
 for (i in 1:nrow(sorted)) {
  if (any(g == sorted[i,1])) {g=g}
  else {g=matrix(c(g,sorted[i,1]), ncol=1)}
 }
NoF=nrow(g)
dirout.fc = paste(getwd(), "/Univariate/Fold_Changes/", sep="")
 dir.create(dirout.fc)
dirout.vol = paste(getwd(), "/Univariate/Volcano_Plots/", sep="")
 dir.create(dirout.vol)
for (i in 1:NoF) {
  for (j in 1:NoF) {
   if (i < j) {
    ni=paste("r.",ExcName(i,slink),".csv",sep="")
    nj=paste("r.",ExcName(j,slink),".csv",sep="")
    pwdi = paste(getwd(), "/Univariate/Groups/", ni, sep="")
    pwdj = paste(getwd(), "/Univariate/Groups/", nj, sep="")
    pv = paste(getwd(), "/Univariate/Pvalues/Pvalues_", 
               ExcName(i,slink), "vs", ExcName(j,slink), ".csv", sep="")
    I=read.csv(pwdi, header=TRUE)
    J=read.csv(pwdj, header=TRUE)
    I = I[,-1]
    J = J[,-1]
meanI = matrix(colMeans(I), ncol=ncol(I))
meanJ = matrix(colMeans(J), ncol=ncol(J))
MeanI = matrix(rep(NA, ncol(I)), nrow=1)
MeanJ = matrix(rep(NA, ncol(I)), nrow=1)
  for (m in 1:ncol(I)) {
    if (meanI[,m] < 0 | meanJ[,m] < 0) {
      MeanI[,m] = 1
      MeanJ[,m] = 1
    } else {
	MeanI[,m] = meanI[,m]
	MeanJ[,m] = meanJ[,m]
      }
  }
	
FC = matrix(MeanI/MeanJ, nrow=ncol(I))
rownames(FC) = colnames(I)
fc.csvfile = paste("Fold_Change_", ExcName(i,slink), "vs",
                   ExcName(j,slink), ".csv", sep="")
write.csv(FC, paste(dirout.fc, fc.csvfile, sep=""))
PV = read.csv(pv, header=TRUE)
PV = matrix(PV[,-1], ncol=1)
logfc = log2(FC)
logpv = -log10(PV)

upper <- log2(upper.lim)
lower <- log2(lower.lim)
sig <- -log10(sig.lim)
colpv=matrix(rep(NA, nrow(PV)), ncol=1)
  for (p in 1:nrow(PV)) {
    if (logfc[p,] < lower | logfc[p,] > upper) {
      if (logpv[p,] > sig) {
      colpv[p,] = "navy"
      } else {
      colpv[p,] = "dark grey"
	}
    } else { 
      if (logpv[p,] > sig) {
        colpv[p,] = "#D55E00"
      } else {
        colpv[p,] = "dark grey"
      }
      }
    }
  
max.fc = 1.3*(max(abs(logfc)))
V = paste(dirout.vol, "VolcanoPlot_", ExcName(i,slink), "vs", 
          ExcName(j,slink), ".pdf", sep="")
pospv=matrix(rep(NA, nrow(PV)), ncol=1)
    for (p in 1:nrow(PV)) {
        if (logfc[p,] < 0) {
            pospv[p,] = 2
        } else { 
        pospv[p,] = 4
    }
   }
pdf(V)
graphics::plot(logfc, logpv, col=colpv, cex= 1.2,pch = 19, 
               xlim=c(-max.fc,max.fc), xlab = "Log2 (Fold Change)", 
               ylab = "Log10 (Pvalue)", 
               main = paste("Volcano Plot ",ExcName(i,slink), 
                            " vs ", ExcName(j,slink), sep=""), 
               sub = paste("(Variables in Blue are significant 
                           (Pvalue<",sig.lim, ") and showed Fold Changes>",
                           upper.lim," or <",lower.lim,")", sep = ""))
if(length(colnames(sorted.x)[grep("navy",colpv)]) > 0){
text(logfc[grep("navy",colpv)], logpv[grep("navy",colpv)], 
     labels=colnames(sorted.x)[grep("navy",colpv)], cex=0.6,
     pos=pospv, col=colpv[grep("navy",colpv)])
   } else {
 }
axis(2, at = c(-1,150), pos=c(lower,0), col="blue", lwd=1,lty= 2)
axis(2, at = c(-1,150), pos=c(upper,0), col="blue", lwd=1,lty= 2)
axis(1, at = c(-150,150), pos=c(sig,0), col="blue", lwd=1,lty= 2)
dev.off()

}}}
}
