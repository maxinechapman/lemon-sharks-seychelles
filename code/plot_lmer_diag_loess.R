plot_lmer_diag<-function(x, ...){
  if (all(!inherits(x,"lmerMod"),!inherits(x,"glmerMod"))){
    stop("this is not an (g)lmerMod object")
  }
  df <- suppressWarnings(broom.mixed::augment(x))
  df$id <- 1:nrow(df)
  df$std.resid<-scale(df$.resid)
  df$.sqrt.abs.std.resid<-sqrt(abs(scale(df$.resid)))
  # based on https://github.com/lme4/lme4/issues/693
  df$cooks_dist <- cooks.distance(influence(x))
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = .fitted, y = .resid, ...)) +
    ggplot2::geom_point() +
    ggplot2::stat_smooth(data=df, ggplot2::aes(x = .fitted, y = .resid), method="loess", se=FALSE) +
    ggplot2::labs(x="Fitted values",y="Residuals", title="Residuals vs Fitted") +
    ggplot2::theme(legend.position="none")+
    ggplot2::geom_text(data=df[max_n(abs(df$.resid),n=3),],ggplot2::aes(x = .fitted, y = .resid,label=id),
              col="red",hjust="inward",vjust="inward")
  if (inherits(x,"lmerMod")) {
    p2 <- ggplot2::ggplot(df, ggplot2::aes(x = .fitted, y = .sqrt.abs.std.resid, ...)) + 
    ggplot2::geom_point() +
    ggplot2::stat_smooth(method="loess", se=FALSE) +
    ggplot2::labs(x="fitted values",y="sqrt(abs(std. residuals))", title="Scale-Location") +
    ggplot2::theme(legend.position="none")+
    ggplot2::geom_text(data=df[max_n(df$.sqrt.abs.std.resid,n=3),],
                       ggplot2::aes(x = .fitted, y = .sqrt.abs.std.resid,label=id),col="red",hjust="inward",vjust="inward")

  # to colour code points, we first create a temp data frame with the qq values
  df2<-qqnorm(scale(df$.resid),plot=FALSE)
  this_qqline<- qqline_eq(scale(df$.resid))
  df2<-data.frame(x=df2$x, y=df2$y)
  df2$delta<-abs(df2$y-(this_qqline$int+df2$x*this_qqline$slope)) # get the vertical distance between points and lines
  df2<-cbind(df2,df) ## in case we need other variable to colour code points
  p3 <- ggplot2::ggplot(df2, ggplot2::aes(x=x, y=y,...))  +
    ggplot2::geom_point() + ggplot2::geom_abline(slope=this_qqline$slope,intercept=this_qqline$int)+
    #ggplot2::stat_qq() + ggplot2::stat_qq_line() +
    ggplot2::labs(x="Theoretical quantiles",y="Standardized residuals", title="Normal Q-Q") +
    ggplot2::theme(legend.position="none")+
    ggplot2::geom_text(data=df2[max_n(df2$delta,n=3),],
                       ggplot2::aes(x = x, y = y,label=id),col="red",hjust="inward",vjust="inward")
  }
  p4 <- ggplot2::ggplot(df, ggplot2::aes(x =id, y = cooks_dist, ...))+
    ggplot2::geom_point() +
    ggplot2::labs(x="Sample id",y="Cook's distance", title="Cook's distance") +
    ggplot2::theme(legend.position="none")+
    ggplot2::geom_text(data=df[max_n(df$cooks_dist,n=3),],
                       ggplot2::aes(x =id, y = cooks_dist,label=id),col="red",hjust="inward",vjust="inward")
  if (any(df$cooks_dist>=0.5)){
    p4 <- p4 + ggplot2::geom_hline(yintercept=0.5,linetype="dashed", color="red")
  }
  if (any(df$cooks_dist>=1)){
    p4 <- p4 + ggplot2::geom_hline(yintercept=1, color="red")
  }
  # return all plots for lmer
  if (inherits(x,"lmerMod")) {
  return(gridExtra::grid.arrange(p1, p3, p2, p4, nrow = 2))
  } else { #and only a subset for glmer
    return(gridExtra::grid.arrange(p1, p4, nrow = 1))
  }
    
}

# retun index of max n values in vector
max_n <-function(x,n){
  x >= sort(x, decreasing=T)[n]
  #which(x >= sort(x, decreasing=T)[n], arr.ind=TRUE)
}

# get equation for the line in a qqplot
qqline_eq<-function (y, datax = FALSE, distribution = qnorm, probs = c(0.25, 
                                                            0.75), qtype = 7, ...) 
{
  stopifnot(length(probs) == 2, is.function(distribution))
  y <- as.vector(quantile(y, probs, names = FALSE, type = qtype, 
                          na.rm = TRUE))
  x <- distribution(probs)
  if (datax) {
    slope <- diff(x)/diff(y)
    int <- x[[1L]] - slope * y[[1L]]
  }
  else {
    slope <- diff(y)/diff(x)
    int <- y[[1L]] - slope * x[[1L]]
  }
  return(list(slope=slope,int=int))
}

plot_lmer_rand_diag<-function(x, random, random_label=NULL){
  if (all(!inherits(x,"lmerMod"),!inherits(x,"glmerMod"))){
    stop("this is not an (g)lmerMod object")
  }
  if (is.null(random_label)){
    random_label <- random
  }
  n_params <- ncol(ranef(x)[[random]])
  # random intercept
  par(mfrow=c(n_params,2))
  hist(ranef(x)[[random]][,1],
       main=paste("Random intercepts by",random_label),
       xlab="random intercept")
  qqnorm(ranef(x)[[random]][,1],
         main = paste("Normal Q-Q Plot by",random_label))
  qqline(ranef(x)[[random]][,1])
  # random slopes
  if(ncol(ranef(x)[[random]])>1){
    hist(ranef(x)[[random]][,2],
         main=paste("Random slopes by",random_label),
         xlab="random slopes")
    qqnorm(ranef(x)[[random]][,2],
           main = paste("Normal Q-Q Plot by",random_label))
    qqline(ranef(x)[[random]][,2])    
  }
} 

