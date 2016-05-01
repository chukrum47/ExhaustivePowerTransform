ExhaustivePowerTransform = function(
  x, 
  y = NA, 
  xmin = -3, 
  xmax = 3, 
  xstep = 0.25, 
  ymin = -3, 
  ymax = 3, 
  ystep = 0.25, 
  type = "rsquared", 
  makeHeatmap = T,
  xlab = "X", 
  ylab = "Y",
  main = "", 
  roundDigits = 3) {
  # Check for required packages
  library(nortest)
  library(fBasics)
  library(gplots)
  library(car)
  library(lmtest)
  # Set xrange to be a vector of powers to try
  if (is.factor(x)) {
    xrange=c(1)
    xmin=1
  } else if (xstep==0) {
    xrange=xmin
  } else {
    xrange = seq(xmin, xmax, xstep)
  }
  # Set yrange to be a vector of powers to try
  if (is.factor(y) | length(y) == 1) {
    yrange=c(1)
    ymin=1
  } else if (ystep==0 | type == "norm" | type == "normCVM" | type == "normLillie" | type == "normPearson") {
    yrange=ymin
  } else {
    yrange = seq(ymin, ymax, ystep)
  }
  # Total number of different combinations to try
  total = length(xrange) * length(yrange)
  i = 1
  # Initialize space for results
  res = data.frame(X=numeric(total), Y=numeric(total), Score=numeric(total))
  for (xt in xrange) {
    # Set tx to be the transformed vector of x
    if (is.factor(x)) { # If x is a factor, don't transform
      tx = x
    } else if (xt == 0) { # If x is zero, use log transformation
      tx = log(x)
    } else { # Otherwise use power transformation
      tx = x^xt
    }
    for (yt in yrange) {
      # Set ty to be the transformed vector of y
      if (is.factor(y)) { # If y is a factor, don't transform
        ty = y
      } else if (yt == 0) { # If y is a zero, use log transformation
        ty = log(y)
      } else { # Otherwise use power transformation
        ty = (y)^yt
      }
      # Choose which test to run based on the specified type
      if (type=="rsquared") {
        score = summary(lm(ty ~ tx))$r.squared
      } else if (type=="variance") {
        score = leveneTest(ty~tx)$"Pr(>F)"[1]
      } else if (type=="norm") {
        score = ad.test(tx)$p.value
      } else if (type=="normCVM") {
        score = cvm.test(tx)$p.value
      } else if (type=="normLillie") {
        score = lillie.test(tx)$p.value
      } else if (type=="normPearson") {
        score = pearson.test(tx)$p.value
      } else if (type=="normShapiro") {
        score = sf.test(tx)$p.value
      } else if (type=="residNorm") {
        score = ad.test(lm(ty~tx)$residuals)$p.value
      } else if (type=="residNormCVM") {
        score = cvm.test(lm(ty~tx)$residuals)$p.value
      } else if (type=="residNormLillie") {
        score = lillie.test(lm(ty~tx)$residuals)$p.value
      } else if (type=="residNormPearson") {
        score = pearson.test(lm(ty~tx)$residuals)$p.value
      } else if (type=="residNormShapiro") {
        score = sf.test(lm(ty~tx)$residuals)$p.value
      } else if (type=="residVar") {
        score = bptest(lm(ty~tx))$p.value
      } else if (type=="residVarGQ") {
        score = gqtest(lm(ty~tx))$p.value
      } else if (type=="residVarHMC") {
        score = hmctest(lm(ty~tx))$p.value
        #To add a new possible test, simply add an additional 'else if' statement here, and just be sure to set the score equal to the resulting value of interest.
      } else {
        score = summary(lm(ty~tx))$r.squared
      }
      # Add values to results matrix
      res[i,1] = xt
      res[i,2] = yt
      res[i,3] = score
      i = i+1
    }
  }
  if (makeHeatmap) {
    # Generate heatmap matrix
    if (xstep == 0 | length(xrange) == 1) {
      mat = matrix(rep(res$Score, 2), byrow=F, ncol=2)
      colnames(mat) = c(xmin, xmin) 
      rownames(mat) = yrange
    } else if (ystep == 0 | length(yrange) == 1) {
      mat = matrix(rep(res$Score, 2), byrow=T, nrow=2)
      colnames(mat) = xrange
      rownames(mat) = c(ymin, ymin)
    } else {
      mat = matrix(res$Score, byrow=F, ncol = length(xrange), nrow=length(yrange))
      colnames(mat) = xrange
      rownames(mat) = yrange
    }
    # Find max score to add lines on the heatmap at those positions
    max.row = which.max(res$Score)
    maxx = which(xrange==res[max.row,1])
    maxy = length(yrange)-which(yrange==res[max.row,2])+1
    #return (mat)
    # Must assign the variables as global variables due to a bug in the 'add.expr' parameter of heatmap which fails to check callers's environment and only checks the global environment. Names are chosen to be suitably long/unique so as to minimize chance of overwritign a global variable.
    assign(x="ExhaustivePowerTransformXTemp", value=maxx, envir=.GlobalEnv)
    assign(x="ExhaustivePowerTransformYTemp", value=maxy, envir=.GlobalEnv)
    # Only create space for title if one is requested
    if (main=="") {
      heatmap.2(mat, 
                Rowv=NA, 
                Colv=NA, 
                dendrogram="none", 
                symkey=F, 
                symbreaks=F, 
                density.info="none", 
                keysize=0.5, 
                key.title=NA, 
                key.xlab="Score", 
                key.par=list(mar=c(3.5,0,3,0)), 
                lmat=rbind(c(5, 4, 2), c(6, 1, 3)), 
                lhei=c(1.2, 5), lwid=c(1, 10, 1), 
                margins=c(4,4),
                col=rev(rainbow(max( length(table(round(mat,roundDigits))), 5), 
                                start = 0/6, 
                                end = 5/6)), 
                scale="none", 
                revC=F, 
                xlab = xlab, 
                ylab = ylab,  
                add.expr=abline(v=ExhaustivePowerTransformXTemp, 
                                h=ExhaustivePowerTransformYTemp), 
                trace="none")
    } else {
      heatmap.2(mat, 
                Rowv=NA, 
                Colv=NA, 
                dendrogram="none", 
                symkey=F, 
                symbreaks=F, 
                density.info="none", 
                keysize=0.5, 
                key.title=NA, 
                key.xlab="Score", 
                key.par=list(mar=c(3.5,0,3,0)), 
                lmat=rbind(c(7,8,9),c(5, 4, 2), c(6, 1, 3)), 
                lhei=c(0.2,1.2, 5), 
                lwid=c(1, 10, 1), 
                margins=c(4,4), 
                col=rev(rainbow(max( length(table(round(mat,roundDigits))), 5), 
                                start = 0/6, 
                                end = 5/6)), 
                scale="none", 
                revC=F, 
                xlab = xlab, 
                ylab = ylab,  
                add.expr=abline(v=ExhaustivePowerTransformXTemp, 
                                h=ExhaustivePowerTransformYTemp), 
                trace="none")
      # Must separately create title due to peculiarities about how the heatmap is laid out
      title(main)
    }
    # Remove global environment variables after assignment
    with(.GlobalEnv, rm(ExhaustivePowerTransformXTemp))
    with(.GlobalEnv, rm(ExhaustivePowerTransformYTemp))
  }
  # Return the result matrix sorted high-to-low
  return(res[sort(res$Score, decreasing=T, index.return=T)$ix,])
}