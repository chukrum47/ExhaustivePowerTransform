<center><h2>Exhaustive Power Transformation</h2></center>


### Description
Exhaustively try all power transformations, within a given range, for a data set.

### Usage
```
ExhaustivePowerTransform( x, 
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
                          roundDigits = 3 )
```

### Arguments
`x` A vector of values. If a model is created, this vector will become the independent variable. Must be specified.  

`y` A vector of values. If a model is created, this vector will become the dependent variable. Default is NA.  

`xmin` The minimum power transformation to be applied to the `x` vector. If `xstep=0`, then `xmin` is the chosen power. Default is -3.  

`xmax` The maximum power transformation to be applied to the `x` vector. Default is 3.  

`xstep` The increment between successive powers to be applied to `x`. If set to 0, then `xmin` is the only power that is used. Default is 0.25.  

`ymin` The minimum power transformation to be applied to the `y` vector. If `ystep=0`, then `ymin` is the chosen power. Does nothing if no `y` vector is specified. Default is -3.  

`ymax` The maximum power transformation to be applied to the `y` vector. Does nothing if no `y` vector is specified. Default is 3.  

`ystep` The increment between successive powers to be applied to `y`. If set to 0, then `ymin` is the only power that is used. Does nothing if no `y` vector is specified. Default is 0.25.  

`type` The criteria by which the attempted powers should be ranked in descending order. Available options are: `"rsquared"`, `"variance"`, `"norm"`, `"normCVM"`, `"normLillie"`, `"normPearson"`, `"normShapiro"`, `"residNorm"`, `"residNormCVM"`, `"residNormLillie"`, `"residNormPearson"`, `"residNormShapiro"`, `"residVar"`, `"residVarGQ"`, and `"residVarHMC"`. Information about each of these types can be found in the details. Default is `"rsquared"`.  

`makeHeatmap` Determines whether a heatmap should be drawn. Default is `True`.

`xlab` The label for the x-axis of the heatmap. Default is "X".  

`ylab` The label for the y-axis of the heatmap. Default is "Y".  

`main` The title of the heatmap. Default is "" (blank).  

`roundDigits` The number of decimals to which each score (as determined by the `type`) is rounded, for purposes of choosing colors for the heatmap. Default is 3.  

### Details

The following packages are required for the function to work: `nortest`, `fBasics`, `gplots`, `car`, and `lmtest`.

For each attempted power transformation, the resulting transformation is graded based on the `type` parameter. The following are the possible options for `type`:  

* `"rsquared"` - Creates a linear model of the form `lm(y~x)` and grades the model based on the r^2 value.  

* `"variance"` - Tests for equal variance with Levene's Test (`leveneTest(y~x)`), and grades based on the p-value of the test.   Note: `x` must be categorical; no transformations are applied to `x`.  

* `"norm"` - Tests for normality of `x` with the Anderson-Darling test for normality, and grades based on the p-value of the test. Note: if this type is chosen, `y` is ignored.  

* `"normCVM"` - Tests for normality of `x` with the Cramer von-Mises test for normality, and grades based on the p-value of the test. Note: if this type is chosen, `y` is ignored.  

* `"normLillie"` - Tests for normality of `x` with the Lilliefors (Kolmogorov-Smirnov) test for normality, and grades based on the p-value of the test. Note: if this type is chosen, `y` is ignored.  
* `"normPearson"` - Tests for normality of `x` with the Pearson chi-square test for normality, and grades based on the p-value of the test. Note: if this type is chosen, `y` is ignored.  
 
* `"normShapiro"` - Tests for normality of `x` with the Shapiro-Francia test for normality, and grades based on the p-value of the test. Note: if this type is chosen, `y` is ignored.  

* `"residNorm"` - Creates a linear model of the form `lm(y~x)` and tests the residuals for normality with the Anderson-Darling test for normality, and grades based on the p-value of the test.  

* `"residNormCVM"` - Creates a linear model of the form `lm(y~x)` and tests the residuals for normality with the Cramer von-Mises test for normality, and grades based on the p-value of the test.  

* `"residNormLillie"` - Creates a linear model of the form `lm(y~x)` and tests the residuals for normality with the Lilliefors (Kolmogorov-Smirnov) test for normality, and grades based on the p-value of the test.  

* `"residNormPearson"` - Creates a linear model of the form `lm(y~x)` and tests the residuals for normality with the Pearson chi-square test for normality, and grades based on the p-value of the test.  

* `"residNormShapiro"` - Creates a linear model of the form `lm(y~x)` and tests the residuals for normality with the Shapiro-Francia test for normality, and grades based on the p-value of the test.  

* `"residVar"` - Creates a linear model of the form `lm(y~x)` and tests the residuals for constant variance with the Breusch-Pagan test for heteroskedasticity, and grades based on the p-value of the test.  

* `"residVarGQ"` - Creates a linear model of the form `lm(y~x)` and tests the residuals for constant variance with the Goldfeld-Quandt test for heteroskedasticity, and grades based on the p-value of the test.  

* `"residVarHMC"` - Creates a linear model of the form `lm(y~x)` and tests the residuals for constant variance with the Harrison-McCabe test for heteroskedasticity, and grades based on the p-value of the test. Note: this test takes significantly longer than other tests, and so caution should be used when running on a large number of transformations.  

If `x` or `y` are detected to be factors, they are not transformed. If only `x` should be transformed, `ystep` should be set to 0 and `ymin` to 1; similarly so in reverse if only `y` should be transformed. Whenever only one vector is being transformed, the resultant heatmap will still appear to show two values for that transformation: this is a limitation of the `heatmap.2` function which does not allow a heatmap to be created with only a single row or a single column, and thus the same row or column is simply repeated. This does not prevent the function from working as intended, and beyond appearance has no effect.  

If `y` is not specified, only `"norm"`, `"normCVM"`, `"normLillie"`, and `"normPearson"` types may be used.  

If an error occurs reading "figure margins too small", delete the partially created plot, increase the size of the plot window, and re-run the function.  

Whenever the power transformation is 0, a natural log transformation is applied instead. 

Disclaimer: this function is designed to be a tool to help automate the process of finding a useful power for transformations, and is not meant to be a replacement for manual analysis: all suggested power transformations should be manually checked as well to ensure. The onus is on the user to ensure that a given transformation is genuinely useful, as opposed to simply passing a test.

### Value

A data frame is returned with three columns: `X`, `Y`, and `Score`, where `Score` is the resulting value (either r^2 or p-value) of the test when run on vectors `x` and `y` using the powers given in column `X` and column `Y`. The data frame is sorted in descending order by the `Score` column, and so the 'optimal' transformation is in the first row.  

If `makeHeatmap` is set to `True`, then a heatmap is also created, with a horizontal and vertical line added to indicate the location of the maximum score.  

### Examples

```{r eval=F, tidy=T}
#Load Boston dataset from the MASS package
boston = MASS::Boston

library(mosaic)
#Look at raw data
xyplot(boston$nox ~ boston$dis, xlab=list(label="Distance", cex=2), ylab=list(label="NOX", cex=2), cex=2, cex.axis=2, scales=list(cex=2))

#Clear transformation necessary, so lets find which works best
result = ExhaustivePowerTransform(boston$dis, boston$nox, main="Heatmap of Transforms (r^2)")
head(result)

#It appears that with x=0.25 and y=-2.5, r^2 = 0.7784, so let's graph the transformation
nox=boston$nox^-2.5
dis = boston$dis^0.25
xyplot(-nox ~ dis, xlab=list(label="Distance^0.25", cex=2), ylab=list(label="-NOX^-2.5", cex=2), cex=2, cex.axis=2, scales=list(cex=2), type=c('p','r'))
#Data now appears linear

#Larger range and larger step
result = ExhaustivePowerTransform(boston$dis, boston$nox, main="Heatmap of Transforms (r^2)", xmin=-4, xmax=4, ymin=-4, ymax=4, xstep=0.5, ystep=0.5)

#Smaller range and smaller step
result = ExhaustivePowerTransform(boston$dis, boston$nox, main="Heatmap of Transforms (r^2)", xmin=-2, xmax=2, ymin=-2, ymax=2, xstep=0.1, ystep=0.1)

#Rounding the r^2 to one decimal
result = ExhaustivePowerTransform(boston$dis, boston$nox, main="Heatmap of Transforms (r^2)", roundDigits=1)

#Rounding the r^2 to two decimals
result = ExhaustivePowerTransform(boston$dis, boston$nox, main="Heatmap of Transforms (r^2)", roundDigits=2)

#Only transform the 'nox' variable looking for normality
result = ExhaustivePowerTransform(boston$nox, main="Heatmap of Transforms (norm)", type="norm")

#Only transform the 'nox' variable looking for normality with the Shapiro-Francia test
result = ExhaustivePowerTransform(boston$nox, main="Heatmap of Transforms (normShapiro)", type="normShapiro")

#Look for normal residuals with Lilliefors test
result = ExhaustivePowerTransform(boston$dis, boston$nox, main="Heatmap of Transforms (residNormLillie)", xstep=0.1, ymin=-5, ymax=0, ystep=0.1, type="residNormLillie")

#Look for homoskedastic residuals with Breusch-Pagan test
result = ExhaustivePowerTransform(boston$dis, boston$nox, main="Heatmap of Transforms (residVar)", xstep=0.1, ystep=0.1, type="residVar")
```
