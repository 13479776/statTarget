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