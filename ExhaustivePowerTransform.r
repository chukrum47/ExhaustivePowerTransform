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
  if (is.factor(x)) { # Don't transform if x is a factor
    xrange=c(1)
    xmin=1
  } else if (xstep==0) {
    xrange=xmin
  } else {
    xrange = seq(xmin, xmax, xstep)
  }
  # Set yrange to be a vector of powers to try
  if (is.factor(y) | length(y) == 1) { # Don't transform is y is a factor or if y is not specified (ie. is.na(y))
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
    if (xstep == 0 | length(xrange) == 1) { # If x has only a single value, we have to make two columns because heatmap.2 will not work with a single column
      mat = matrix(rep(res$Score, 2), byrow=F, ncol=2)
      colnames(mat) = c(xmin, xmin) 
      rownames(mat) = yrange
    } else if (ystep == 0 | length(yrange) == 1) { # If y has only a single value, we have to make two rows because heatmap.2 will not work with a single row
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
    # Must convert the coordinates to the heatmap.2 scale
    maxx = which(xrange==res[max.row,1]) 
    maxy = length(yrange)-which(yrange==res[max.row,2])+1
    # Must assign the variables as global variables due to a bug in the 'add.expr' parameter of heatmap which fails to check callers's environment and only checks the global environment. Names are chosen to be suitably long/unique so as to minimize chance of overwriting a global variable.
    assign(x="ExhaustivePowerTransformXTemp", value=maxx, envir=.GlobalEnv)
    assign(x="ExhaustivePowerTransformYTemp", value=maxy, envir=.GlobalEnv)
    # Only create space for title if one is requested
    if (main=="") {
      heatmap.2(mat, 
                Rowv=NA, # Don't reorder rows
                Colv=NA, # Don't reorder columns
                dendrogram="none", # Don't create dendrograms
                symkey=F, # Don't draw extra stuff on heatmap
                symbreaks=F, # Don't draw extra stuff on heatmap
                density.info="none", # Don't draw density lines on color key
                keysize=0.5, # Reduce side of color key
                key.title=NA, # Don't add title to color key
                key.xlab="Score", # Add x-axis label for color key
                key.par=list(mar=c(3.5,0,3,0)), # Set margins on color key element
                lmat=rbind(c(5, 4, 2), c(6, 1, 3)), # Create extra buffering to left of heatmap/color key to make it appear centered
                lhei=c(1.2, 5), # Set the relative heights of the color key and heatmap
                lwid=c(1, 10, 1), # Set the width of the heatmap/color key in respect to the left and right edges
                margins=c(4,4), # Set the margins of the graph
                col=rev(rainbow(max( length(table(round(mat,roundDigits))), 5), # Set the number of different colors 
                                start = 0/6, 
                                end = 5/6)), # Set the color palette to be rainbow
                scale="none", # Don't rescale plot
                revC=F, # Don't reverse columns' order
                xlab = xlab, # Set x-axis label
                ylab = ylab,  # Set y-axis label
                add.expr=abline(v=ExhaustivePowerTransformXTemp, 
                                h=ExhaustivePowerTransformYTemp), # Add black lines to mark highest location
                trace="none") # Don't draw extra stuff on heatmap
    } else {
      heatmap.2(mat, 
                Rowv=NA, # Don't reorder rows
                Colv=NA, # Don't reorder columns
                dendrogram="none", # Don't create dendrograms
                symkey=F, # Don't draw extra stuff on heatmap 
                symbreaks=F, # Don't draw extra stuff on heatmap 
                density.info="none", # Don't draw density lines on color key
                keysize=0.5, # Reduce side of color key
                key.title=NA, # Don't add title to color key
                key.xlab="Score", # Add x-axis label for color key
                key.par=list(mar=c(3.5,0,3,0)), # Set margins on color key element
                lmat=rbind(c(7,8,9),c(5, 4, 2), c(6, 1, 3)), # Create extra buffering to left of heatmap/color key to make it appear centered and extra room above for title
                lhei=c(0.2,1.2, 5), # Set the relative heights of the title, color key, and heatmap
                lwid=c(1, 10, 1), # Set the width of the heatmap/color key in respect to the left and right edges
                margins=c(4,4), # Set the margins of the graph
                col=rev(rainbow(max( length(table(round(mat,roundDigits))), 5), # Set the number of different colors  
                                start = 0/6, 
                                end = 5/6)), # Set the color palette to be rainbow
                scale="none", # Don't rescale plot
                revC=F, # Don't reverse columns' order
                xlab = xlab, # Set x-axis label
                ylab = ylab, # Set y-axis label
                add.expr=abline(v=ExhaustivePowerTransformXTemp, 
                                h=ExhaustivePowerTransformYTemp), # Add black lines to mark highest location
                trace="none")# Don't draw extra stuff on heatmap
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