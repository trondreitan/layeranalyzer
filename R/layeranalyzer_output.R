
#####################################################
#
# Functions for summarizing analysis results and
# printing such  summaries.
# 
# Trond Reitan, 17. Aug. 2018
#
#####################################################


# Summary data object:

summary.layered=function(object, ...)
{
  input=list(...)

  #show("summary called")

  if(is.null(object$ML.loglik))
  {
    mat=matrix(NA,nrow=length(object$parameter.names),ncol=4)
    for(i in 1:length(object$parameter.names))
     for(j in 1:4)
      mat[i,j]=as.numeric(object[[i+6]][j])
    dimnames(mat)=list(object$parameter.names, c("Mean","Median",
     "Lower 95%", 
     "Upper 95%"))
    
    ans=list(coefficients=round(mat,6), model.log.likelihood=object$model.log.lik)
    class(ans)="summary.layered"
    
    return(ans)
  }  
  if(!is.null(object$ML.loglik))
  {
    mat=matrix(NA,nrow=length(object$parameter.names),ncol=3)
    for(i in 1:length(object$parameter.names))
     for(j in 1:3)
      mat[i,j]=as.numeric(object[[i+10]][c(1,4,5)[j]])
    dimnames(mat)=list(object$parameter.names, c("ML estimate",
     "Bayesian Lower 95%", "Bayesian Upper 95%"))
    
    ans=list(coefficients=round(mat,6), ML.loglik=object$ML.loglik, 
      AIC=object$AIC,AICc=object$AICc,BIC=object$BIC)
    class(ans)="summary.layered"
    
    return(ans)
  }  
}


# Function for printing summary:

print.summary.layered=function(x, ...)
{
  input=list(...)
  
  if(is.null(x$ML.loglik))
  {
    print.srcref("Coefficients:")
    print(x$coefficients)
    if(!is.na(x$model.log.likelihood))
    {
      print.srcref("")
      print.srcref(sprintf("Model log-likelihood: %9.3f", x$model.log.likelihood))
    }
  }
  
  if(!is.null(x$ML.loglik))
  {
    print.srcref("Coefficients:")
    print(x$coefficients)
    
    print.srcref("")
    print.srcref(sprintf("Log-likelihood: %9.3f", x$ML.loglik))
    print.srcref(sprintf("AIC : %9.3f", x$AIC))
    print.srcref(sprintf("AICc: %9.3f", x$AICc))
    print.srcref(sprintf("BIC : %9.3f", x$BIC))
  }
}




