---
title: "Estimation using the Propensity Score"
author: "Luis Alvarez"
date: "17/09/2020"
output: 
  html_document:
    keep_md: true
---



In this session, we will discuss three different estimation methods for the propensity score: subclassification, matching on the estimated propensity score and inverse probability weighting. In what follows, we assume a random sample $\{Y_i,X_i,D_i\}_{i=1}^n$ from the population of interest.

## Subclassification on the estimated propensity score

This method consists in dividing the sample in blocks using the algorithm in the previous session (though with option TRIM set to FALSE), and then, within each block, to compute an estimate of the treatment effect. The intuition here is that, since there is much less variation in propensity scores within each block, comparisons between treatment and control group are more credible (if unconfoundedness given the propensity score holds).

If our sample is divided in $B$ blocks, then Imbens and Rubin suggest estimating treatment effects within each block, $\tau_b$, by fitting a linear model for each block $b=1,2\ldots B$:

$$Y_{ib} = \alpha_b + \tau_b D_{ib} + \gamma_b' A_{ib} + \epsilon_{ib}$$

where $Y_{ib}$ is the outcome of observation $i$ in block $b$; and  $A_{ib}$ is a subset of covariates included for bias reduction purposes (we may not include any covariates if we wish; in this case, the estimator reduces to a comparison of means within each block). The idea here is that any remaining imbalances in covariates within each block can be accounted for if we ``control'' for some covariates. Since observations within each block are rather homogeneous in terms of covariates, including these in a linear fashion should not lead to the extrapolation biases that a linear regression using the whole sample usually incurs.

Finally, after estimating each treatment effect within each block, we aggreate these as:

\begin{equation}
\hat{\tau} = \sum_{b=1}^B \omega_b \hat{\tau}_b \quad (*)
\end{equation}
where the weights depend on the estimand of interest. If interest lies in the average treatment effect (ATE), then $\omega_b = n_b/n$, where $n_b$ is the number of observations in group $b$. If interest is in the average treatment effect on the treated (ATT), then $\omega_b = n_b^t/n^t$, where $n_b^t$ is the number of treated units in block $b$; and $n^t$ is the overall number of treated units.

Notice that, since the partition in blocks depends only on covariates and the estimated propensity score, and since each observation lies within a single block; a conditional variance formula for the estimator in $(*)$ is readily available. Indeed, defining $\mathbf{X} := [X_1, X_2, \ldots X_n]$ and $\mathbf{D} := [D_1, D_2, \ldots D_n]$, we have:

\begin{equation}
\mathbb{V}[\hat{\tau}|\mathbf{D}, \mathbf{X}] = \sum_{b=1}^B \omega_b^2  \mathbb{V}[\hat{\tau}_b|\mathbf{D}, \mathbf{X}]
\end{equation}

Estimates for $\mathbb{V}[\hat{\tau}_b|\mathbf{D}, \mathbf{X}]$ are readily available from each "within-block" regression output. Imbens and Rubin suggest using heteroskedasticity-robust standard errors in estimating $\mathbb{V}[\hat{\tau}|\mathbf{D}, \mathbf{X}]$. Importantly, since the variance formula conditions on covariates and treatment indicators, it abstracts from variability due to estimation error in the estimated propensity score (since we are conditioning, this error only generates bias). Moreover, because it conditions on covariates, it also abstracts from variability in the conditional treatment effect $\tau(x) = \mathbb{E}[Y_i(1)-Y_i(0)|X_i=x]$, which should be accounted for if we are interested in estimating the treatment effect in the population from which our sample was drawn (intuitively, we need to take into account the difference between the distribution of the control variables in our sample and its distribution in the population, as it may well be that our weights $\omega$ differ from their population counterparts). The second source of variability always increases the variance, whereas the first one may actually decrease it (more on this later). These two sources of variability are actually not a problem if interest lies in the sample treatment effect (cf. Chapter 19 in Imbens and Rubin for a discussion), though we will discuss them in detail for the next two estimators.

Let's implement the subclassification method in our example from the previous session. First, we just repeat some crucial steps from the previous session.




```r
#Loading package
library("Matching")
```

```
## Loading required package: MASS
```

```
## Warning: package 'MASS' was built under R version 3.6.2
```

```
## ## 
## ##  Matching (Version 4.9-7, Build Date: 2020-02-05)
## ##  See http://sekhon.berkeley.edu/matching for additional documentation.
## ##  Please cite software as:
## ##   Jasjeet S. Sekhon. 2011. ``Multivariate and Propensity Score Matching
## ##   Software with Automated Balance Optimization: The Matching package for R.''
## ##   Journal of Statistical Software, 42(7): 1-52. 
## ##
```

```r
data(lalonde)

#Estimating the propensity score model
formula.model = imbens.rubin.stepwise("treat",c("nodegr","black","educ"),c("age","re74","re75","u74","u75","married"), lalonde)
ps.imbens <- glm(formula.model, family = "binomial", data = lalonde)
lalonde$ps.linear = predict(ps.imbens, type = "link")

#Trimming to improve balance
keep.indicators = trimming.imbens(lalonde$ps.linear)
```

```
## [1] "Trimming threshold alpha is  0.13308191816384"
```

```r
lalonde.trimmed = lalonde[keep.indicators,]
```

Next, we construct a function that estimates our effects using blocking with regression. This function uses the block construction function presented in the previous session:


```r
#Robust standard errors package
library("sandwich")

#Blocking with regresison
#Outcome - character with outcome name
#Treatment - character with treatment indicator
#contorls = vector with names of variables to be included for covariate adjustment. You may leave it as "c()" if you do not wish to include any controls.
#data = data.frame with observations
#lin.score = vector with estimates of the linearized propensity scores
#estimand of interest: either ATT or ATE
#... you may pass additional arguments to the block construction function if you wish (e.g. change T+max)
blocking.regression <- function(outcome, treatment, controls, data, lin.score, estimand = "ATE", ...){
  K = length(controls)
  
  if(K==0)
    formula.est = paste(outcome, treatment, sep = " ~ ") else formula.est = paste(outcome, paste(c(treatment, controls), collapse=" + "), sep = " ~ ")
  
  treat.vec = data[,treatment]
  #Constructing blocks for subclassification
  block.var = propensity.score.blocks(treat.vec,lin.score, K, trim =F, ...)
  
  print(paste("Sample was divided in", length(unique(block.var)), "blocks"))
  
  #print(block.var)
  est = c()
  var = c()
  weights = c()
  block.results = c()
  
  for(j in unique(block.var))
  {
    #print(block.var == j)
    model.block = lm(formula.est, data = data[block.var == j,])
    
    estimate.block = model.block$coefficients[treatment]
    est = c(est, estimate.block)
    
    vcov = vcovHC(model.block)
    var.block = vcov[treatment, treatment]
    var = c(var, var.block)
    
    if(estimand== "ATE")
      weights = c(weights, sum(block.var==j)/length(block.var))
      else if(estimand == "ATT")
        weights = c(weights, sum(treat.vec[block.var==j])/sum(treat.vec))
    
    block.results[[j]] = list("Block no" = j, "Nobs" = sum(block.var==j), "Ntreat" = sum(treat.vec[block.var==j]), "Est" = estimate.block, "SE" = sqrt(var.block))
  }
  
  eff.whole = sum(weights*est)
  var.whole = sum((weights^2)*var)
  
  
  return(list("Estimand" = estimand, "Estimate" = eff.whole, "SE" = sqrt(var.whole), "Details per block" = block.results))
    
}
```

Next, we apply the above function in our sample


```r
#Estimating ATE - no adjustment for covariates
blocking.regression("re78", "treat", c(), lalonde.trimmed, lalonde.trimmed$ps.linear)
```

```
## [1] "Sample was divided in 2 blocks"
```

```
## $Estimand
## [1] "ATE"
## 
## $Estimate
## [1] 1617.683
## 
## $SE
## [1] 656.5996
## 
## $`Details per block`
## $`Details per block`[[1]]
## $`Details per block`[[1]]$`Block no`
## [1] 1
## 
## $`Details per block`[[1]]$Nobs
## [1] 219
## 
## $`Details per block`[[1]]$Ntreat
## [1] 69
## 
## $`Details per block`[[1]]$Est
##    treat 
## 548.8325 
## 
## $`Details per block`[[1]]$SE
## [1] 863.9399
## 
## 
## $`Details per block`[[2]]
## $`Details per block`[[2]]$`Block no`
## [1] 2
## 
## $`Details per block`[[2]]$Nobs
## [1] 225
## 
## $`Details per block`[[2]]$Ntreat
## [1] 115
## 
## $`Details per block`[[2]]$Est
##    treat 
## 2658.031 
## 
## $`Details per block`[[2]]$SE
## [1] 985.7469
```

```r
#Estimating ATT - adjusting for earnings in 75
blocking.regression("re78", "treat", c("re75"), lalonde.trimmed, lalonde.trimmed$ps.linear, estimand = "ATT")
```

```
## [1] "Sample was divided in 2 blocks"
```

```
## $Estimand
## [1] "ATT"
## 
## $Estimate
## [1] 1815.821
## 
## $SE
## [1] 697.5131
## 
## $`Details per block`
## $`Details per block`[[1]]
## $`Details per block`[[1]]$`Block no`
## [1] 1
## 
## $`Details per block`[[1]]$Nobs
## [1] 219
## 
## $`Details per block`[[1]]$Ntreat
## [1] 69
## 
## $`Details per block`[[1]]$Est
##    treat 
## 581.1108 
## 
## $`Details per block`[[1]]$SE
## [1] 864.4793
## 
## 
## $`Details per block`[[2]]
## $`Details per block`[[2]]$`Block no`
## [1] 2
## 
## $`Details per block`[[2]]$Nobs
## [1] 225
## 
## $`Details per block`[[2]]$Ntreat
## [1] 115
## 
## $`Details per block`[[2]]$Est
##    treat 
## 2556.648 
## 
## $`Details per block`[[2]]$SE
## [1] 988.163
```

## Matching on the estimated propensity score

Next, we consider matching on the estimated propensity score. The idea here is to compare the outcome of each unit (all units if interest lies in ATE; treated units if interest lies in ATT) with the average outcome of the closest $M$ units, in terms of the propensity score, in the opposite group. Typically $M$ is small, often $M=1$ (Abadie and Imbens, 2016). We may also match units with or without replacement; in the latter, an observation can serve as a counterfactual to a *single* observation in the opposite group.

Package Matching of Sekhon (2011) implements propensity score matching (and other forms of matching) under a variety of options. See the author's page (and his JSS paper in the references) for several examples: http://sekhon.berkeley.edu/matching/ . Below, we will implement matching on our estimated propensity score:


```r
#Estimation
match  = Match(Y = lalonde.trimmed$re78, Tr = lalonde.trimmed$treat, X = lalonde.trimmed$ps.linear, estimand = "ATT", M=1, replace= T)
summary(match)
```

```
## 
## Estimate...  1746.8 
## AI SE......  830.11 
## T-stat.....  2.1043 
## p.val......  0.035351 
## 
## Original number of observations..............  444 
## Original number of treated obs...............  184 
## Matched number of observations...............  184 
## Matched number of observations  (unweighted).  2464
```
By default, the matching package reports Abadie and Imbens (2006) standard errors, which account for variability in the conditional treatment effect $\tau(x)$ (i.e. it implicitly assumes the estimand of interest is the ATE/ATT in the population from which our sample was drawn). If we want to report conditional standard errors like in the previous session (i.e. interest lies in the sample ATE/ATT), we just make:


```r
summary(match, full = T)
```

```
## 
## Estimate...  1746.8 
## AI SE......  830.11 
## T-stat.....  2.1043 
## p.val......  0.035351 
## 
## Est noAdj..  1746.8 
## SE.........  709.89 
## T-stat.....  2.4607 
## p.val......  0.013867 
## 
## Original number of observations..............  444 
## Original number of treated obs...............  184 
## Matched number of observations...............  184 
## Matched number of observations  (unweighted).  2464
```
Importantly, though, the Abadie and Imbens (2006) standard errors reported by the package do *not* account for variability due to estimation error in the propensity score. This is not much of a problem if interest lies in the population ATE, as Abadie and Imbens (2016) actually show that estimating the propensity score *reduces* the variance, so we may think of our variance formula as conservative. However, in the ATT case, we may be under- or overestimating the variance, so you may want to check their paper to see how to account for this error.

Package Matching further allows for bias-reduction techniques by specifying a set of covariates as argument $Z$. See the documentation previously referenced and Chapter 18 of the book for further discussions.

## Inverse probability weighting

Finally, we consider inverse probability weighting. This method uses the observation that:

\begin{equation}
{\tau}^{ATE} = \mathbb{E}[Y_i(1) - Y_i(0)] = \mathbb{E}\left[\frac{Y_i \cdot D_i}{\mathbb{P}[D_i=1|X_i]} - \frac{Y_i \cdot (1-D_i)}{1-\mathbb{P}[D_i=1|X_i]}\right] 
\end{equation}

\begin{equation}
{\tau}^{ATT} = \mathbb{E}[Y_i(1) - Y_i(0)|D_i=1] = \mathbb{E}\left[\frac{Y_i \cdot D_i}{\mathbb{P}[D_i=1]} - \frac{Y_i \cdot (1-D_i)\cdot\mathbb{P}[D_i=1|X_i]}{\mathbb{P}[D_i=1](1-\mathbb{P}[D_i=1|X_i])}\right] 
\end{equation}

to construct analog estimators such as 

\begin{equation}
\hat{\tau}^{ATE} = \frac{1}{n}\sum_{i=1}^n \frac{Y_i \cdot D_i}{\hat{p}_i} - \frac{1}{n}\sum_{i=1}^n \frac{Y_i \cdot (1-D_i)}{1-\hat{p}_i}
\end{equation}
where $\hat{p}_i$ are the predicted conditional probabilities of treatment from the propensity score model. A similar formula holds for the ATT, with $\hat{p}$ -- the proportion of treated units in the sample -- entering the formula.

Notice that $\mathbb{E}[D_i/\mathbb{P}[D_i=1|X_i]] = \mathbb{E}[(1-D_i)/(1-\mathbb{P}[D_i=1|X_i])]=1$. Thus, to force estimated weights to sum up to unity, one usually considers normalizing them, i.e. we consider

\begin{equation}
\tilde{\tau}^{ATE} =  \frac{1}{\sum_{i=1}^N D_i \hat{p}_i}\sum_{i=1}^n \frac{Y_i \cdot D_i}{\hat{p}_i} - \frac{1}{\sum_{i=1}^n (1-D_i)(1-\hat{p}_i)}\sum_{i=1}^n \frac{Y_i \cdot (1-D_i)}{1-\hat{p}_i}
\end{equation}

This estimator can be easily implemented by running a *weighted* regression:

\begin{equation}
Y_i = \alpha + \beta D_i + \epsilon_i
\end{equation}
with weights $w_i = (\hat{p}_i)^{-D_i} \cdot (1-\hat{p}_i)^{-(1-D_i)}$. This formulation is useful, because it allows further covariate adjustment by considering a weighted regression:

\begin{equation}
Y_i = \alpha + \beta D_i + \gamma' A_i + \epsilon_i
\end{equation}
with weights  $w_i = (\hat{p}_i)^{-D_i} \cdot (1-\hat{p}_i)^{-(1-D_i)}$. This approach is quite common in applied practice. It enjoys a property known as *double robustness*: the estimated treatment effect will be consistent for the parameter of interest if either the propensity score is correctly specified, or the linear specification is correct. 

The regression formulation of the IPW is also useful because it yields immediate variance estimators from the regression output. Importantly, it can be shown that using robust standard errors in the formulation above already yields an estimator of the variance that accounts for variability in $\tau(x)$ (Hirano, Imbens and Ridder, 2003). As for variability due to estimation of the propensity score, we note that not accounting for the estimation of the propensity score actually leads to *conservative* inference on population ATT/ATE (Hirano, Imbens and Ridder, 2003), since estimation of the propensity score only increases efficiency.


```r
#IPW for ATE estimation
weights = (lalonde.trimmed$treat/plogis(lalonde.trimmed$ps.linear)) + ((1-lalonde.trimmed$treat)/(1-plogis(lalonde.trimmed$ps.linear)))

ipw.model = lm(re78~treat, data = lalonde.trimmed, weights = weights)

#Package to produce coefficient tests
library(lmtest)
```

```
## Loading required package: zoo
```

```
## Warning: package 'zoo' was built under R version 3.6.2
```

```
## 
## Attaching package: 'zoo'
```

```
## The following objects are masked from 'package:base':
## 
##     as.Date, as.Date.numeric
```

```r
#Inference using HC standard errors
coeftest(ipw.model, vcov. = vcovHC)
```

```
## 
## t test of coefficients:
## 
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  4549.53     345.94 13.1514  < 2e-16 ***
## treat        1551.38     656.95  2.3615  0.01863 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Finally, we note that the inverse in IPW may actually lead to problems: indeed, small biases in the propensity score may lead to large biases in estimated effects. Imbens and Rubin argue that subclassification is prefferable because it "smooths" inverse weights by assigning the same weight to observations in the same block. Cf. Section 17.8 in the book for a discussion.


## References

### Book

Imbens, G., & Rubin, D. (2015). Causal Inference for Statistics, Social, and Biomedical Sciences: An Introduction. Cambridge: Cambridge University Press. doi:10.1017/CBO9781139025751

### Articles

Abadie, A. and Imbens, G.W. (2006), Large Sample Properties of Matching Estimators for Average Treatment Effects. Econometrica, 74: 235-267. doi:10.1111/j.1468-0262.2006.00655.x

Abadie, A. and Imbens, G.W. (2016), Matching on the Estimated Propensity Score. Econometrica, 84: 781-807. doi:10.3982/ECTA11293

Hirano, K., Imbens, G.W. and Ridder, G. (2003), Efficient Estimation of Average Treatment Effects Using the Estimated Propensity Score. Econometrica, 71: 1161-1189. doi:10.1111/1468-0262.00442

Sekhon, J. (2011). Multivariate and Propensity Score Matching Software with Automated Balance Optimization: The Matching package for R. Journal of Statistical Software, 42(7), 1 - 52. doi:http://dx.doi.org/10.18637/jss.v042.i07
