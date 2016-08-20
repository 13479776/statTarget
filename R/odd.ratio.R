#' this function provide the odd.ratio value for metabolomics data or others.
#' @param file is the data imput
#' @details 
#' An odds ratio (OR) is a measure of association between an exposure and an outcome. OR=1 Exposure does not affect odds of outcome; OR>1 Exposure associated with higher odds of outcome; OR<1 Exposure associated with lower odds of outcome.The 95% confidence interval (CI) is used to estimate the precision of the OR. A large CI indicates a low level of precision of the OR, whereas a small CI indicates a higher precision of the OR. It is important to note however, that unlike the p value, the 95% CI does not report a measure???s statistical significance. 
#' @references 
#' J Can Acad Child Adolesc Psychiatry. 2010 Aug; 19(3): 227???229. Magdalena Szumilas
#' @usage odd.ratio(file)
#' @export
odd.ratio <- function(file) {
  pwdfile = paste(getwd(), "/Univariate/DataTable.csv", sep = "")
  file = pwdfile
  x <- read.csv(file, sep = ",", header = TRUE)
  x.x = x[, 3:ncol(x)]
  rownames(x.x) = x[, 2]
  k = matrix(x[, 1], ncol = 1)
  x.n = cbind(k, x.x)
  sorted = x.n[order(x.n[, 1]), ]
  g = c()
  for (i in 1:nrow(sorted)) {
    if (any(g == sorted[i, 1])) {
      g = g
    }
    else {
      g = matrix(c(g, sorted[i, 1]), ncol = 1)
    }
  }
  NoF = nrow(g)
  dirout.w = paste(getwd(), "/Univariate/oddratio", sep = "")
  dir.create(dirout.w)
  for (i in 1:NoF) {
    for (j in 1:NoF) {
      if (i < j) {
        ni = paste("r.", i, ".csv", sep = "")
        nj = paste("r.", j, ".csv", sep = "")
        pwdi = paste(getwd(), "/Univariate/Groups/", ni, sep = "")
        pwdj = paste(getwd(), "/Univariate/Groups/", nj, sep = "")
        I = read.csv(pwdi, header = TRUE)
        J = read.csv(pwdj, header = TRUE)
        I = I[, -1]
        J = J[, -1]
        Ilf = matrix(rep(i),nrow(I))
        Jlf = matrix(rep(j),nrow(J))    
        colnames(Ilf) = c("lf")
        colnames(Jlf) = c("lf")
        I = cbind(Ilf,I)
        J = cbind(Jlf,J)
        IJ = rbind(I,J)
        #fin = ncol(sorted) - 1
        #we <- matrix(rep(NA, fin))
        #odds.radio..........................................................................................................................
        IJM <- as.matrix(IJ[,2:ncol(IJ)])
        outFactor <-as.factor(IJ[,1])
        logit.or <- function(outcome,a){
                    glm.log <- glm(outcome~scale(a),family=binomial(link="logit"))
                    res <- summary(glm.log)$coefficients[2,1:2]
                    p <- summary(glm.log)$coefficients[2,4]
                    out <- c(res,p)
                    Odd.radio <- exp(out[1])
                    Std.error1 <- exp(out[1]-1.96*out[2])
                    Std.error2 <- exp(out[1]+1.96*out[2])
                    Out <- c(Odd.radio,Std.error1,Std.error2,p)
                    names(Out)=c("Odd.radio","2.5%CI","97.5%CI","p.value")
                    return(Out)
        }
        or <- apply(t(IJM),1,function(x){logit.or(outFactor,x)})
        or <- as.data.frame(t(or))
        #write.table(or,"odds_radio.txt",sep="\t",quote=F)
        or.ij = paste("odds_radio_", i, "vs", j, ".csv", sep = "")
        assign(or.ij, or)
        write.csv(or, paste(dirout.w, or.ij, sep = "/"))
      }
    }
  }
}
        