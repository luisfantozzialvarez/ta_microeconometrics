---
title: "Estimating the propensity score"
author: "Luis Alvarez"
date: "9/8/2020"
output:
  html_document:
    keep_md: true
---



## The selection on observables framework

We consider a binary treatment $D_i$ assigned to units $i$ in the population of interest. Let $Y_i(0)$ and $Y_i(1)$ denote the potential outcomes associated with this treatment. The observed outcome is $Y_i = Y_i(0) (1-D_i)  + Y_i(1)D_i$. We will work under the assumption that there exists a vector covariates $X_i$ such that, in the population from which we draw our samples:

\begin{equation}
Y_i(0),Y_i(1) \perp D_i | X_i \quad \text{(*)}
\end{equation}

and, for some $0 < \underline{p} < \bar{p} < 1$ 

\begin{equation}
\underline{p} \leq e(x) \leq \bar{p} \quad \forall x \in \mathcal{X} \quad \text{(**)}
\end{equation}
where $\mathcal{X}$ denotes the support of $X$; and $e(x) = \mathbb{P}[D_i=1|X_i=x]$.

You have seen in class that assumptions $(*)$ and $(**)$ enable identification of both ATT and ATE. Crucially, under the above assumptions:

\begin{equation}
Y_i(0), Y_i(1) \perp D_i | e(X_i) \quad \text{(***)}
\end{equation}

i.e. it sufficies to "control", in an appropriate manner, for the propensity score in order to identify ATE and ATT. Indeed, in light of this observation, most estimation methods you have seen in class involve estimating the propensity score and adjusting for it in a "smart" manner.

In this note, we will show how to perform analyses under assumptions $\text{(*)}$ and $\text{(**)}$. We will follow Imbens and Rubin's approach, which proceeds in three steps. In the first step, we try to estimate a "good" propensity score; and try to ensure that $\text{(**)}$ holds, possibly by restricting our sample to ensure we are drawing from a subpopulation where overlap seems more plausible. Importantly, we do not deal with outcomes $Y_i$ in this section; all we care is between the relationship of $X_i$ and $D_i$. In the second step, if possible, we try to assess the credibility of $\text{(*)}$ in our (possibly restricted) sample, by performing some placebo checks. Finally, in the third step, we estimate our causal effects of interest. In this note, we deal with the first step. In the next TA session, we will discuss steps 2 and 3.

Throughout, we work with the Lalonde dataset available on the R package Matching. The dataset store data on several variables and on participation in a job training program.



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

#Outcome of interest is re78
summary(lalonde)
```

```
##       age             educ          black             hisp        
##  Min.   :17.00   Min.   : 3.0   Min.   :0.0000   Min.   :0.00000  
##  1st Qu.:20.00   1st Qu.: 9.0   1st Qu.:1.0000   1st Qu.:0.00000  
##  Median :24.00   Median :10.0   Median :1.0000   Median :0.00000  
##  Mean   :25.37   Mean   :10.2   Mean   :0.8337   Mean   :0.08764  
##  3rd Qu.:28.00   3rd Qu.:11.0   3rd Qu.:1.0000   3rd Qu.:0.00000  
##  Max.   :55.00   Max.   :16.0   Max.   :1.0000   Max.   :1.00000  
##     married           nodegr           re74              re75      
##  Min.   :0.0000   Min.   :0.000   Min.   :    0.0   Min.   :    0  
##  1st Qu.:0.0000   1st Qu.:1.000   1st Qu.:    0.0   1st Qu.:    0  
##  Median :0.0000   Median :1.000   Median :    0.0   Median :    0  
##  Mean   :0.1685   Mean   :0.782   Mean   : 2102.3   Mean   : 1377  
##  3rd Qu.:0.0000   3rd Qu.:1.000   3rd Qu.:  824.4   3rd Qu.: 1221  
##  Max.   :1.0000   Max.   :1.000   Max.   :39570.7   Max.   :25142  
##       re78            u74              u75             treat       
##  Min.   :    0   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
##  1st Qu.:    0   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000  
##  Median : 3702   Median :1.0000   Median :1.0000   Median :0.0000  
##  Mean   : 5301   Mean   :0.7326   Mean   :0.6494   Mean   :0.4157  
##  3rd Qu.: 8125   3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000  
##  Max.   :60308   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000
```

## Initial balancing tests

As an initial step, it is good to assess the balance in some covariates we think may be relevant for identification. This will give us an idea on how difficult it will be to estimate a good propensity score; and whether such task will be fruitful or not.

We can do so by using the MatchBalance function in the Matching package.


```r
MatchBalance(treat~age+educ+black+nodegr+re74+re75+u74+u75+married, data = lalonde)
```

```
## 
## ***** (V1) age *****
## before matching:
## mean treatment........ 25.816 
## mean control.......... 25.054 
## std mean diff......... 10.655 
## 
## mean raw eQQ diff..... 0.94054 
## med  raw eQQ diff..... 1 
## max  raw eQQ diff..... 7 
## 
## mean eCDF diff........ 0.025364 
## med  eCDF diff........ 0.022193 
## max  eCDF diff........ 0.065177 
## 
## var ratio (Tr/Co)..... 1.0278 
## T-test p-value........ 0.26594 
## KS Bootstrap p-value.. 0.504 
## KS Naive p-value...... 0.7481 
## KS Statistic.......... 0.065177 
## 
## 
## ***** (V2) educ *****
## before matching:
## mean treatment........ 10.346 
## mean control.......... 10.088 
## std mean diff......... 12.806 
## 
## mean raw eQQ diff..... 0.40541 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 2 
## 
## mean eCDF diff........ 0.028698 
## med  eCDF diff........ 0.012682 
## max  eCDF diff........ 0.12651 
## 
## var ratio (Tr/Co)..... 1.5513 
## T-test p-value........ 0.15017 
## KS Bootstrap p-value.. 0.014 
## KS Naive p-value...... 0.062873 
## KS Statistic.......... 0.12651 
## 
## 
## ***** (V3) black *****
## before matching:
## mean treatment........ 0.84324 
## mean control.......... 0.82692 
## std mean diff......... 4.4767 
## 
## mean raw eQQ diff..... 0.016216 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 1 
## 
## mean eCDF diff........ 0.0081601 
## med  eCDF diff........ 0.0081601 
## max  eCDF diff........ 0.01632 
## 
## var ratio (Tr/Co)..... 0.92503 
## T-test p-value........ 0.64736 
## 
## 
## ***** (V4) nodegr *****
## before matching:
## mean treatment........ 0.70811 
## mean control.......... 0.83462 
## std mean diff......... -27.751 
## 
## mean raw eQQ diff..... 0.12432 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 1 
## 
## mean eCDF diff........ 0.063254 
## med  eCDF diff........ 0.063254 
## max  eCDF diff........ 0.12651 
## 
## var ratio (Tr/Co)..... 1.4998 
## T-test p-value........ 0.0020368 
## 
## 
## ***** (V5) re74 *****
## before matching:
## mean treatment........ 2095.6 
## mean control.......... 2107 
## std mean diff......... -0.23437 
## 
## mean raw eQQ diff..... 487.98 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 8413 
## 
## mean eCDF diff........ 0.019223 
## med  eCDF diff........ 0.0158 
## max  eCDF diff........ 0.047089 
## 
## var ratio (Tr/Co)..... 0.7381 
## T-test p-value........ 0.98186 
## KS Bootstrap p-value.. 0.612 
## KS Naive p-value...... 0.97023 
## KS Statistic.......... 0.047089 
## 
## 
## ***** (V6) re75 *****
## before matching:
## mean treatment........ 1532.1 
## mean control.......... 1266.9 
## std mean diff......... 8.2363 
## 
## mean raw eQQ diff..... 367.61 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 2110.2 
## 
## mean eCDF diff........ 0.050834 
## med  eCDF diff........ 0.061954 
## max  eCDF diff........ 0.10748 
## 
## var ratio (Tr/Co)..... 1.0763 
## T-test p-value........ 0.38527 
## KS Bootstrap p-value.. 0.046 
## KS Naive p-value...... 0.16449 
## KS Statistic.......... 0.10748 
## 
## 
## ***** (V7) u74 *****
## before matching:
## mean treatment........ 0.70811 
## mean control.......... 0.75 
## std mean diff......... -9.1895 
## 
## mean raw eQQ diff..... 0.037838 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 1 
## 
## mean eCDF diff........ 0.020946 
## med  eCDF diff........ 0.020946 
## max  eCDF diff........ 0.041892 
## 
## var ratio (Tr/Co)..... 1.1041 
## T-test p-value........ 0.33033 
## 
## 
## ***** (V8) u75 *****
## before matching:
## mean treatment........ 0.6 
## mean control.......... 0.68462 
## std mean diff......... -17.225 
## 
## mean raw eQQ diff..... 0.081081 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 1 
## 
## mean eCDF diff........ 0.042308 
## med  eCDF diff........ 0.042308 
## max  eCDF diff........ 0.084615 
## 
## var ratio (Tr/Co)..... 1.1133 
## T-test p-value........ 0.068031 
## 
## 
## ***** (V9) married *****
## before matching:
## mean treatment........ 0.18919 
## mean control.......... 0.15385 
## std mean diff......... 8.9995 
## 
## mean raw eQQ diff..... 0.037838 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 1 
## 
## mean eCDF diff........ 0.017672 
## med  eCDF diff........ 0.017672 
## max  eCDF diff........ 0.035343 
## 
## var ratio (Tr/Co)..... 1.1802 
## T-test p-value........ 0.33425 
## 
## 
## Before Matching Minimum p.value: 0.0020368 
## Variable Name(s): nodegr  Number(s): 4
```
Imbens and Rubins also suggest reporting the normalized differences:


$$  \frac{\bar{X}_{\text{treat}} - \bar{X}_{\text{control}}}{\sqrt{\frac{S^2_{X,\text{treat}} + S^2_{X,\text{control}}}{2}}} $$
where $S^2_{X,\text{treat}}$ and $S^2_{X,\text{control}}$ are the sample variances in the treatment and control group. Their point is that, with a large enough sample size, a t-test will always have enough power to detect even the smallest differences in means between groups, and that may not always be helpful. In this sense, the measure above, which is in standard deviations, may be more useful in assessing balance.

Let's construct a function that reports both a t-test and the normalized differences.

```r
table.test <- function(covariates, treat_var, data)
{
  treat_vec = data[,treat_var]
  table = c()
  
  for(lab in covariates)
  {
    cov_vec = data[,lab]
    control = cov_vec[treat_vec==0&!is.na(treat_vec)]
    treatment = cov_vec[treat_vec==1&!is.na(treat_vec)]
    
    normalized_diff = (mean(treatment) - mean(control))/sqrt((var(treatment)+var(control))/2)
    
    table.line = cbind( "Mean control" = mean(control), "Mean treatment" = mean(treatment), 
                            "t-stat" = tryCatch({t.test(control, x = treatment)$statistic}, error =function(e){NaN}), "Normalized diff" = normalized_diff)
    
    table = rbind(table, table.line)
  }
  
  rownames(table) = covariates
  
  return(table)
}
```



```r
table.test(c("age","educ","black","nodegr","re74","re75","u74","u75","married"), "treat", lalonde)
```

```
##         Mean control Mean treatment     t-stat Normalized diff
## age       25.0538462     25.8162162  1.1140361     0.107277121
## educ      10.0884615     10.3459459  1.4421840     0.141219821
## black      0.8269231      0.8432432  0.4577777     0.043886611
## nodegr     0.8346154      0.7081081 -3.1084981    -0.303986439
## re74    2107.0268154   2095.5740000 -0.0227466    -0.002159921
## re75    1266.9092408   1532.0556297  0.8692061     0.083863254
## u74        0.7500000      0.7081081 -0.9746890    -0.094140477
## u75        0.6846154      0.6000000 -1.8299743    -0.176809436
## married    0.1538462      0.1891892  0.9668363     0.093640701
```

## Estimating the propensity score

Next, we will estimate the propensity score. The usual method is via logistic regression. Let $Z$ denote a transformation of the vector of covariates $X$. This may include linear interactions, logs and powers of the covariates in $X$. We then estimate the propensity score via MLE:

\begin{equation}
\hat{\gamma} \in \text{argmax}_{b \in \mathbb{R}^d} \sum_{i=1}^N D_i \log(\ell(Z_i'b))  + (1-D_i) \log (1 - \ell(Z_i'b)) \quad \text{(****)}
\end{equation}
where $\ell(\cdot)$ is the logistic link function:

$$\ell(s) = \frac{\exp(s)}{1+\exp(s)}$$
Once we estimate $\gamma$, the predicted conditional probability of treatment for observation $i$ is given by:

\begin{equation}
\hat{p}_i = \ell(\hat{\gamma}'z_i)
\end{equation}
where $\hat{\gamma}'z_i$ is the latent linear index.

The crucial step above is to select which covariates/transformations enter $Z_i$. Recall that the goal is to have a good approximation of $e(X_i) = \mathbb{P}[D_i=1|X_i]$, so we must be flexible enough. However, we don't to be *too* flexible, otherwise we incur in *overfitting*. Indeed, note that, if $Z_i$ is flexible enough, then an upper bound to the maximization problem $\text{(****)}$, which is to set $\hat{p}_i = 1$ if $D_i=1$ and $\hat{p}_i = 0$ if $D_i=0$; may become feasible. Clearly, this is an undesirable solution, as it violates overlap in the estimated propensity score  -- indeed, most of our estimators won't even work under this solution. So we must have a principled method to select variables to be included in the propensity score that is able to balance between the two forces discussed.

In this session, we will discuss two selection methods:

1. Imbens and Rubin's stepwise selction method.
2. The Lasso.

### Imbens and Rubin's stepwise selection method
Imbens and Rubin's method lets us define a ``basic'' set of covariates $X^b$, which we think should always be kept in the specification; and an additional set of covariates, $X^t$, which we will then test for inclusion. After testing for inclusion of $X^t$, we arrive at a vector $X^* = [X^{b'}, X^{a'}]'$, where $X^a$ is the selected subset of $X^t$ (it may be empty). We then test for the inclusion of linear interactions and quadratic terms of variables in $X^*$. All tests for inclusion are based on likelihood ratio statistics; and in defining thresholds in terms of these statistics for inclusion (cf. their book, or Appendix A of Imbens, 2015, for the algorithm).

In the code below, we present a function that implements Imbens and Rubin's approach.


```r
# Imbens and Rubin's stepwise selection algorithm
# treatment: character variable for treatment indicator variable
# Xb: character vector with names of basic covariates: you may pass it as  c() if you do not want any basic covariate
# Xt: character vector with names for covariates to be tested for inclusion
# data: dataframe with variables
# Clinear: threshold, in terms of likelihood ratio statistics, for inclusion of linear terms
# Cquadratic: threshold, in terms of likelihood ratio statistics, for inclusion of quadratic/interaction terms
# Intercept: does model include intercept?
imbens.rubin.stepwise <- function(treatment, Xb, Xt, data, Clinear = 1, Cquadratic = 2.71, intercept = T)
{
  #Add or not intercept
  if(intercept)
    inter.add = "1" else inter.add = "-1"
  
  
  #Formula for model
  if(length(Xb)==0)
    formula = paste(treatment, inter.add, sep = " ~ ") else formula = paste(treatment, paste(c(inter.add,Xb), collapse = " + "), sep = " ~ ")
    
  continue = T
  
  Xt_left = Xt
  
  while(continue){
    null.model = glm(as.formula(formula), data, family = "binomial")
    
    null.lkl = logLik(null.model)
    
    test.stats = c()
    for(covariate in Xt_left)
    {
      formula.test = paste(formula, covariate, sep = " + ")
      test.model = glm(as.formula(formula.test), data, family = "binomial")
      
      lkl.ratio = 2*(as.numeric(logLik(test.model))-as.numeric(null.lkl))
      test.stats = c(test.stats, lkl.ratio)
    }
    
    if(max(test.stats,na.rm = T)<Clinear)
      continue = F else {
        
        add.coef = Xt_left[which.max(test.stats)]
        
        formula = paste(formula, add.coef, sep = " + ")
        
        Xt_left = Xt_left[-which.max(test.stats)]
      }
      
  }
  
  #Defining Xstar set
  Xstar = c(Xb, Xt[!(Xt%in%Xt_left)])
  
  #Creating all combinations of Xstar interactions
  combinations = expand.grid(Xstar, Xstar)
  Xcomb = paste(combinations[,1],combinations[,2],sep=":")

  continue = T
  
  Xcomb_left = Xcomb
  
  while(continue){
    null.model = glm(as.formula(formula), data, family = "binomial")
    
    null.lkl = logLik(null.model)
    
    test.stats = c()
    for(covariate in Xcomb_left)
    {
      formula.test = paste(formula, covariate, sep = " + ")
      test.model = glm(as.formula(formula.test), data, family = "binomial")
      
      lkl.ratio = 2*(as.numeric(logLik(test.model))-as.numeric(null.lkl))
      test.stats = c(test.stats, lkl.ratio)
    }

    if(max(test.stats,na.rm = T)<Cquadratic)
      continue = F else {
        
        add.coef = Xcomb_left[which.max(test.stats)]
        
        formula = paste(formula, add.coef, sep = " + ")
        
        Xcomb_left = Xcomb_left[-which.max(test.stats)]
      }
      
  }
  
  return(formula)
}
```

Using the function above in our data


```r
formula.model = imbens.rubin.stepwise("treat",c("nodegr","black","educ"),c("age","re74","re75","u74","u75","married"), lalonde)

print(formula.model)
```

```
## [1] "treat ~ 1 + nodegr + black + educ + u75 + re74 + educ:nodegr"
```

```r
#Estimating the model
ps.imbens <- glm(formula.model, family = "binomial", data = lalonde)

#Loading Robust (to misspecification) SEs
library(sandwich)
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
coeftest(ps.imbens, vcov. = vcovHC)
```

```
## 
## z test of coefficients:
## 
##                Estimate  Std. Error z value Pr(>|z|)  
## (Intercept) -6.8600e+00  4.1091e+00 -1.6695  0.09502 .
## nodegr       7.7225e+00  4.1974e+00  1.8398  0.06579 .
## black        2.0183e-01  2.7171e-01  0.7428  0.45759  
## educ         5.9601e-01  3.3561e-01  1.7759  0.07575 .
## u75         -5.7152e-01  2.4307e-01 -2.3513  0.01871 *
## re74        -3.0169e-05  2.0414e-05 -1.4778  0.13945  
## nodegr:educ -7.1123e-01  3.4525e-01 -2.0600  0.03939 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
In this step, it is always good to check the robustness of our propensity score model to different choices of $X^b$ and $X^t$; $C_{\text{linear}}$ and $C_{\text{quadratic}}$.

### The Lasso

Lasso works by modifying $\text{(****)}$ by including a *cost-function*:

$$\begin{equation}
\hat{\gamma} \in \text{argmax}_{b \in \mathbb{R}^d} \sum_{i=1}^N D_i \log(\ell(Z_i'b))  + (1-D_i) \log (1 - \ell(Z_i'b))  + \lambda \sum_{i=1}^d |b_i|
\end{equation}$$

for some penalty constant $\lambda > 0$. The idea here is that the penalty introduces a cost for including variable $i$ in the model (i.e. setting $b_i \neq 0$). Since the marginal cost of increasing the value of $|b_i|$ is constant, the Lasso leads to solutions with *exact* zeros, making it a method of variable selection. This is in stark constrast with other "smooth" penalty functions (e.g. the Ridge). The Lasso is known to work well in "sparse" data-generating models, i.e. models where only a few of the relevant population coefficients are different from zero (i.e. only "a few", but *unknown*, transformations of $X$ truly determine the probability of selection).

The package HDM of R implements post-Lasso, where we first select the variables using Lasso, and then run standard logistic regression on the selected components. This is known to ameliorate the attenuation biases induced by Lasso. The package implements data-dependent choices of penalty constant. Let's use it to estimate a propensity score model in our data.



```r
library(hdm)

#I will construct a vector with powers up to order 2 of the covariates:
X = as.matrix(lalonde[,c("educ","nodegr","age","black","re74","re75","u74","u75","married")])

Z = poly(X, degree = 2, raw = T, simple = T)
D = as.vector(lalonde[,"treat"])
head(Z)
```

```
##      1.0.0.0.0.0.0.0.0 2.0.0.0.0.0.0.0.0 0.1.0.0.0.0.0.0.0 1.1.0.0.0.0.0.0.0
## [1,]                11               121                 1                11
## [2,]                 9                81                 1                 9
## [3,]                12               144                 0                 0
## [4,]                11               121                 1                11
## [5,]                 8                64                 1                 8
## [6,]                 9                81                 1                 9
##      0.2.0.0.0.0.0.0.0 0.0.1.0.0.0.0.0.0 1.0.1.0.0.0.0.0.0 0.1.1.0.0.0.0.0.0
## [1,]                 1                37               407                37
## [2,]                 1                22               198                22
## [3,]                 0                30               360                 0
## [4,]                 1                27               297                27
## [5,]                 1                33               264                33
## [6,]                 1                22               198                22
##      0.0.2.0.0.0.0.0.0 0.0.0.1.0.0.0.0.0 1.0.0.1.0.0.0.0.0 0.1.0.1.0.0.0.0.0
## [1,]              1369                 1                11                 1
## [2,]               484                 0                 0                 0
## [3,]               900                 1                12                 0
## [4,]               729                 1                11                 1
## [5,]              1089                 1                 8                 1
## [6,]               484                 1                 9                 1
##      0.0.1.1.0.0.0.0.0 0.0.0.2.0.0.0.0.0 0.0.0.0.1.0.0.0.0 1.0.0.0.1.0.0.0.0
## [1,]                37                 1                 0                 0
## [2,]                 0                 0                 0                 0
## [3,]                30                 1                 0                 0
## [4,]                27                 1                 0                 0
## [5,]                33                 1                 0                 0
## [6,]                22                 1                 0                 0
##      0.1.0.0.1.0.0.0.0 0.0.1.0.1.0.0.0.0 0.0.0.1.1.0.0.0.0 0.0.0.0.2.0.0.0.0
## [1,]                 0                 0                 0                 0
## [2,]                 0                 0                 0                 0
## [3,]                 0                 0                 0                 0
## [4,]                 0                 0                 0                 0
## [5,]                 0                 0                 0                 0
## [6,]                 0                 0                 0                 0
##      0.0.0.0.0.1.0.0.0 1.0.0.0.0.1.0.0.0 0.1.0.0.0.1.0.0.0 0.0.1.0.0.1.0.0.0
## [1,]                 0                 0                 0                 0
## [2,]                 0                 0                 0                 0
## [3,]                 0                 0                 0                 0
## [4,]                 0                 0                 0                 0
## [5,]                 0                 0                 0                 0
## [6,]                 0                 0                 0                 0
##      0.0.0.1.0.1.0.0.0 0.0.0.0.1.1.0.0.0 0.0.0.0.0.2.0.0.0 0.0.0.0.0.0.1.0.0
## [1,]                 0                 0                 0                 1
## [2,]                 0                 0                 0                 1
## [3,]                 0                 0                 0                 1
## [4,]                 0                 0                 0                 1
## [5,]                 0                 0                 0                 1
## [6,]                 0                 0                 0                 1
##      1.0.0.0.0.0.1.0.0 0.1.0.0.0.0.1.0.0 0.0.1.0.0.0.1.0.0 0.0.0.1.0.0.1.0.0
## [1,]                11                 1                37                 1
## [2,]                 9                 1                22                 0
## [3,]                12                 0                30                 1
## [4,]                11                 1                27                 1
## [5,]                 8                 1                33                 1
## [6,]                 9                 1                22                 1
##      0.0.0.0.1.0.1.0.0 0.0.0.0.0.1.1.0.0 0.0.0.0.0.0.2.0.0 0.0.0.0.0.0.0.1.0
## [1,]                 0                 0                 1                 1
## [2,]                 0                 0                 1                 1
## [3,]                 0                 0                 1                 1
## [4,]                 0                 0                 1                 1
## [5,]                 0                 0                 1                 1
## [6,]                 0                 0                 1                 1
##      1.0.0.0.0.0.0.1.0 0.1.0.0.0.0.0.1.0 0.0.1.0.0.0.0.1.0 0.0.0.1.0.0.0.1.0
## [1,]                11                 1                37                 1
## [2,]                 9                 1                22                 0
## [3,]                12                 0                30                 1
## [4,]                11                 1                27                 1
## [5,]                 8                 1                33                 1
## [6,]                 9                 1                22                 1
##      0.0.0.0.1.0.0.1.0 0.0.0.0.0.1.0.1.0 0.0.0.0.0.0.1.1.0 0.0.0.0.0.0.0.2.0
## [1,]                 0                 0                 1                 1
## [2,]                 0                 0                 1                 1
## [3,]                 0                 0                 1                 1
## [4,]                 0                 0                 1                 1
## [5,]                 0                 0                 1                 1
## [6,]                 0                 0                 1                 1
##      0.0.0.0.0.0.0.0.1 1.0.0.0.0.0.0.0.1 0.1.0.0.0.0.0.0.1 0.0.1.0.0.0.0.0.1
## [1,]                 1                11                 1                37
## [2,]                 0                 0                 0                 0
## [3,]                 0                 0                 0                 0
## [4,]                 0                 0                 0                 0
## [5,]                 0                 0                 0                 0
## [6,]                 0                 0                 0                 0
##      0.0.0.1.0.0.0.0.1 0.0.0.0.1.0.0.0.1 0.0.0.0.0.1.0.0.1 0.0.0.0.0.0.1.0.1
## [1,]                 1                 0                 0                 1
## [2,]                 0                 0                 0                 0
## [3,]                 0                 0                 0                 0
## [4,]                 0                 0                 0                 0
## [5,]                 0                 0                 0                 0
## [6,]                 0                 0                 0                 0
##      0.0.0.0.0.0.0.1.1 0.0.0.0.0.0.0.0.2
## [1,]                 1                 1
## [2,]                 0                 0
## [3,]                 0                 0
## [4,]                 0                 0
## [5,]                 0                 0
## [6,]                 0                 0
```

```r
ps.lasso = rlassologit(x = Z, y =  D)

summary(ps.lasso)
```

```
## 
## Call:
## rlassologit.default(x = Z, y = D)
## 
## Post-Lasso Estimation:  TRUE 
## 
## Total number of variables: 54
## Number of selected variables: 2 
## 
##                   Estimate
## (Intercept)          0.269
## 1.0.0.0.0.0.0.0.0    0.000
## 2.0.0.0.0.0.0.0.0    0.000
## 0.1.0.0.0.0.0.0.0    0.000
## 1.1.0.0.0.0.0.0.0   -0.054
## 0.2.0.0.0.0.0.0.0    0.000
## 0.0.1.0.0.0.0.0.0    0.000
## 1.0.1.0.0.0.0.0.0    0.000
## 0.1.1.0.0.0.0.0.0    0.000
## 0.0.2.0.0.0.0.0.0    0.000
## 0.0.0.1.0.0.0.0.0    0.000
## 1.0.0.1.0.0.0.0.0    0.000
## 0.1.0.1.0.0.0.0.0    0.000
## 0.0.1.1.0.0.0.0.0    0.000
## 0.0.0.2.0.0.0.0.0    0.000
## 0.0.0.0.1.0.0.0.0    0.000
## 1.0.0.0.1.0.0.0.0    0.000
## 0.1.0.0.1.0.0.0.0    0.000
## 0.0.1.0.1.0.0.0.0    0.000
## 0.0.0.1.1.0.0.0.0    0.000
## 0.0.0.0.2.0.0.0.0    0.000
## 0.0.0.0.0.1.0.0.0    0.000
## 1.0.0.0.0.1.0.0.0    0.000
## 0.1.0.0.0.1.0.0.0    0.000
## 0.0.1.0.0.1.0.0.0    0.000
## 0.0.0.1.0.1.0.0.0    0.000
## 0.0.0.0.1.1.0.0.0    0.000
## 0.0.0.0.0.2.0.0.0    0.000
## 0.0.0.0.0.0.1.0.0    0.000
## 1.0.0.0.0.0.1.0.0    0.000
## 0.1.0.0.0.0.1.0.0    0.000
## 0.0.1.0.0.0.1.0.0    0.000
## 0.0.0.1.0.0.1.0.0    0.000
## 0.0.0.0.1.0.1.0.0    0.000
## 0.0.0.0.0.1.1.0.0    0.000
## 0.0.0.0.0.0.2.0.0    0.000
## 0.0.0.0.0.0.0.1.0    0.000
## 1.0.0.0.0.0.0.1.0    0.000
## 0.1.0.0.0.0.0.1.0   -0.412
## 0.0.1.0.0.0.0.1.0    0.000
## 0.0.0.1.0.0.0.1.0    0.000
## 0.0.0.0.1.0.0.1.0    0.000
## 0.0.0.0.0.1.0.1.0    0.000
## 0.0.0.0.0.0.1.1.0    0.000
## 0.0.0.0.0.0.0.2.0    0.000
## 0.0.0.0.0.0.0.0.1    0.000
## 1.0.0.0.0.0.0.0.1    0.000
## 0.1.0.0.0.0.0.0.1    0.000
## 0.0.1.0.0.0.0.0.1    0.000
## 0.0.0.1.0.0.0.0.1    0.000
## 0.0.0.0.1.0.0.0.1    0.000
## 0.0.0.0.0.1.0.0.1    0.000
## 0.0.0.0.0.0.1.0.1    0.000
## 0.0.0.0.0.0.0.1.1    0.000
## 0.0.0.0.0.0.0.0.2    0.000
```


## Assessing the quality of the estimated propensity score

The next step is to assess the quality of the estimated propensity score. Recall that the propensity score is a *balancing score*, in the sense that:

\begin{equation}
X_i \perp W_i | e(X_i)
\end{equation}

So we could try to test this assumption using our estimates. Indeed, Imbens and Rubin suggest partitioning the sample in *sub-blocks*, where the estimated propensity score is approximately constant; and to conduct balancing tests within each block. They provide us with an algorithm to construct the subblocks, where we iteratively subdivide the sample based on the median propensity score. Cf the algorithm in Chapter 13 of their book. 

Below, we implement their algorithm, given a sequence of estimated linear indices $Z_i'\hat{\gamma}$ of the propensity score.


```r
#Function that subdivides a given propensity score vector in subblocks
#treat = vector with treatment assignments
#lin.psm = vector with linearized PSs
#K = how many covariates will we want to test/use in bias correction of estimates later on? 
#t.max = threshold for tstat in making a further subdivide 
#trim = should we discard extreme observations so there is overlap?
propensity.score.blocks <- function(treat,lin.psm, K, t.max = 1.96,  trim =T)
{
  if(trim){
    b0 = min(plogis(lin.psm[treat==1]))
    b1 = max(plogis(lin.psm[treat==0]))
  } else
  {
    b0 = 0
    b1 = 1
  }
  b_vec =c(b0,b1)
  while(TRUE)
  {
    J = length(b_vec)-1
    b_vec_new = do.call(c,lapply(1:J, function(j){
      sample = (b_vec[j] <= plogis(lin.psm)) & (plogis(lin.psm) < b_vec[j+1])
      
      ps.treat = lin.psm[sample&treat==1]
      ps.control = lin.psm[sample&treat==0]
      
      #print(length(ps.control))
      #print(length(ps.treat))
      
      t.test.pass = tryCatch({abs(t.test(ps.control, ps.treat)$statistic) > t.max}, error = function(e){return(F)})
      
      med.val = median(c(ps.treat, ps.control))
      
      Nt.below = sum(ps.treat < med.val)
      Nt.above = sum(ps.treat >= med.val)
      Nc.below = sum(ps.control < med.val)
      Nc.above = sum(ps.control >= med.val)
      
      sample.crit = min(Nt.below, Nt.above, Nc.below, Nc.above) >= max(3, K+2)
      
      if(t.test.pass&sample.crit)
        return(c(b_vec[j], plogis(med.val), b_vec[j+1])) else return(c(b_vec[j], b_vec[j+1]))
      
    }))
    b_vec_new = unique(b_vec_new)
    
    #print(length(b_vec_new))
    if(length(b_vec_new)==length(b_vec))
      break else b_vec = b_vec_new
  }
  
  #Constructing blocking variable now
  block_var = rep(NA, length(treat))
  
  for(j in 1:(length(b_vec)-1))
    block_var[(b_vec[j] <= plogis(lin.psm)) & (plogis(lin.psm) < b_vec[j+1])] = j
  
  return(block_var)
}
```

Applying it to our data:



```r
lalonde$linear.ps.stepwise = predict(ps.imbens, type = "link")
lalonde$linear.ps.lasso = predict(ps.lasso, type = "link")

lalonde$blocks.stepwise = propensity.score.blocks(lalonde$treat, lalonde$linear.ps.stepwise, K =  length(c("age","educ","black","nodegr","re74","re75","u74","u75","married")))

lalonde$blocks.lasso = propensity.score.blocks(lalonde$treat, lalonde$linear.ps.lasso, K =  length(c("age","educ","black","nodegr","re74","re75","u74","u75","married")))

head(lalonde)
```

```
##   age educ black hisp married nodegr re74 re75     re78 u74 u75 treat
## 1  37   11     1    0       1      1    0    0  9930.05   1   1     1
## 2  22    9     0    1       0      1    0    0  3595.89   1   1     1
## 3  30   12     1    0       0      0    0    0 24909.50   1   1     1
## 4  27   11     1    0       0      1    0    0  7506.15   1   1     1
## 5  33    8     1    0       0      1    0    0   289.79   1   1     1
## 6  22    9     1    0       0      1    0    0  4056.49   1   1     1
##   linear.ps.stepwise linear.ps.lasso blocks.stepwise blocks.lasso
## 1        -0.77469498      -0.7330867               1            1
## 2        -0.74607979      -0.6258434               1            2
## 3        -0.07757216       0.2690294               2           NA
## 4        -0.77469498      -0.7330867               1            1
## 5        -0.42903088      -0.5722217               2            2
## 6        -0.54425225      -0.6258434               1            2
```

Let's now calculate balance within each blocL


```r
for(j in unique(lalonde$blocks.stepwise)[order(unique(lalonde$blocks.stepwise))])
 if(!is.na(j))
 {
   print(paste("Block",j,sep= " "))
   MatchBalance(treat~age+educ+black+nodegr+re74+re75+u74+u75+married, data = lalonde[lalonde$blocks.stepwise==j&!is.na(lalonde$blocks.stepwise),])
 }
```

```
## [1] "Block 1"
## 
## ***** (V1) age *****
## before matching:
## mean treatment........ 25.13 
## mean control.......... 24.703 
## std mean diff......... 6.4394 
## 
## mean raw eQQ diff..... 1.058 
## med  raw eQQ diff..... 1 
## max  raw eQQ diff..... 10 
## 
## mean eCDF diff........ 0.030325 
## med  eCDF diff........ 0.01994 
## max  eCDF diff........ 0.094753 
## 
## var ratio (Tr/Co)..... 0.80817 
## T-test p-value........ 0.67192 
## KS Bootstrap p-value.. 0.528 
## KS Naive p-value...... 0.79529 
## KS Statistic.......... 0.094753 
## 
## 
## ***** (V2) educ *****
## before matching:
## mean treatment........ 10.145 
## mean control.......... 10.034 
## std mean diff......... 11.336 
## 
## mean raw eQQ diff..... 0.14493 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 1 
## 
## mean eCDF diff........ 0.027686 
## med  eCDF diff........ 0.0047976 
## max  eCDF diff........ 0.11964 
## 
## var ratio (Tr/Co)..... 1.1313 
## T-test p-value........ 0.431 
## KS Bootstrap p-value.. 0.152 
## KS Naive p-value...... 0.51508 
## KS Statistic.......... 0.11964 
## 
## 
## ***** (V3) black *****
## before matching:
## mean treatment........ 0.81159 
## mean control.......... 0.81379 
## std mean diff......... -0.55824 
## 
## mean raw eQQ diff..... 0 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 0 
## 
## mean eCDF diff........ 0.0010995 
## med  eCDF diff........ 0.0010995 
## max  eCDF diff........ 0.0021989 
## 
## var ratio (Tr/Co)..... 1.0169 
## T-test p-value........ 0.96953 
## 
## 
## ***** (V4) nodegr *****
## before matching:
## mean treatment........ 1 
## mean control.......... 1 
## std mean diff......... 0 
## 
## mean raw eQQ diff..... 0 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 0 
## 
## mean eCDF diff........ 0 
## med  eCDF diff........ 0 
## max  eCDF diff........ 0 
## 
## var ratio (Tr/Co)..... NaN 
## T-test p-value........ 1 
## 
## 
## ***** (V5) re74 *****
## before matching:
## mean treatment........ 1421.6 
## mean control.......... 1444.7 
## std mean diff......... -0.51172 
## 
## mean raw eQQ diff..... 281.84 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 5956.7 
## 
## mean eCDF diff........ 0.0090142 
## med  eCDF diff........ 0.0075462 
## max  eCDF diff........ 0.031084 
## 
## var ratio (Tr/Co)..... 0.97492 
## T-test p-value........ 0.97225 
## KS Bootstrap p-value.. 0.918 
## KS Naive p-value...... 1 
## KS Statistic.......... 0.031084 
## 
## 
## ***** (V6) re75 *****
## before matching:
## mean treatment........ 721.05 
## mean control.......... 811.93 
## std mean diff......... -2.7804 
## 
## mean raw eQQ diff..... 338.31 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 8195.6 
## 
## mean eCDF diff........ 0.01556 
## med  eCDF diff........ 0.014493 
## max  eCDF diff........ 0.04068 
## 
## var ratio (Tr/Co)..... 1.3247 
## T-test p-value........ 0.8433 
## KS Bootstrap p-value.. 0.664 
## KS Naive p-value...... 1 
## KS Statistic.......... 0.04068 
## 
## 
## ***** (V7) u74 *****
## before matching:
## mean treatment........ 0.86957 
## mean control.......... 0.84828 
## std mean diff......... 6.2754 
## 
## mean raw eQQ diff..... 0.028986 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 1 
## 
## mean eCDF diff........ 0.010645 
## med  eCDF diff........ 0.010645 
## max  eCDF diff........ 0.021289 
## 
## var ratio (Tr/Co)..... 0.88805 
## T-test p-value........ 0.67467 
## 
## 
## ***** (V8) u75 *****
## before matching:
## mean treatment........ 0.89855 
## mean control.......... 0.88276 
## std mean diff......... 5.1925 
## 
## mean raw eQQ diff..... 0.028986 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 1 
## 
## mean eCDF diff........ 0.0078961 
## med  eCDF diff........ 0.0078961 
## max  eCDF diff........ 0.015792 
## 
## var ratio (Tr/Co)..... 0.88757 
## T-test p-value........ 0.72836 
## 
## 
## ***** (V9) married *****
## before matching:
## mean treatment........ 0.18841 
## mean control.......... 0.13103 
## std mean diff......... 14.565 
## 
## mean raw eQQ diff..... 0.057971 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 1 
## 
## mean eCDF diff........ 0.028686 
## med  eCDF diff........ 0.028686 
## max  eCDF diff........ 0.057371 
## 
## var ratio (Tr/Co)..... 1.3533 
## T-test p-value........ 0.30018 
## 
## 
## Before Matching Minimum p.value: 0.152 
## Variable Name(s): educ  Number(s): 2 
## 
## [1] "Block 2"
## 
## ***** (V1) age *****
## before matching:
## mean treatment........ 25.937 
## mean control.......... 25.62 
## std mean diff......... 4.3437 
## 
## mean raw eQQ diff..... 0.75 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 7 
## 
## mean eCDF diff........ 0.024508 
## med  eCDF diff........ 0.025025 
## max  eCDF diff........ 0.067568 
## 
## var ratio (Tr/Co)..... 1.1567 
## T-test p-value........ 0.73945 
## KS Bootstrap p-value.. 0.82 
## KS Naive p-value...... 0.964 
## KS Statistic.......... 0.067568 
## 
## 
## ***** (V2) educ *****
## before matching:
## mean treatment........ 10.279 
## mean control.......... 10.13 
## std mean diff......... 6.504 
## 
## mean raw eQQ diff..... 0.2037 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 1 
## 
## mean eCDF diff........ 0.019728 
## med  eCDF diff........ 0.017768 
## max  eCDF diff........ 0.061812 
## 
## var ratio (Tr/Co)..... 1.1147 
## T-test p-value........ 0.62162 
## KS Bootstrap p-value.. 0.712 
## KS Naive p-value...... 0.98497 
## KS Statistic.......... 0.061812 
## 
## 
## ***** (V3) black *****
## before matching:
## mean treatment........ 0.85586 
## mean control.......... 0.86111 
## std mean diff......... -1.4895 
## 
## mean raw eQQ diff..... 0.0092593 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 1 
## 
## mean eCDF diff........ 0.0026276 
## med  eCDF diff........ 0.0026276 
## max  eCDF diff........ 0.0052553 
## 
## var ratio (Tr/Co)..... 1.0312 
## T-test p-value........ 0.91167 
## 
## 
## ***** (V4) nodegr *****
## before matching:
## mean treatment........ 0.55856 
## mean control.......... 0.62037 
## std mean diff......... -12.392 
## 
## mean raw eQQ diff..... 0.064815 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 1 
## 
## mean eCDF diff........ 0.030906 
## med  eCDF diff........ 0.030906 
## max  eCDF diff........ 0.061812 
## 
## var ratio (Tr/Co)..... 1.0467 
## T-test p-value........ 0.35476 
## 
## 
## ***** (V5) re74 *****
## before matching:
## mean treatment........ 2609 
## mean control.......... 1925.2 
## std mean diff......... 13.258 
## 
## mean raw eQQ diff..... 655.71 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 12181 
## 
## mean eCDF diff........ 0.041307 
## med  eCDF diff........ 0.042167 
## max  eCDF diff........ 0.074575 
## 
## var ratio (Tr/Co)..... 1.7115 
## T-test p-value........ 0.27084 
## KS Bootstrap p-value.. 0.606 
## KS Naive p-value...... 0.92105 
## KS Statistic.......... 0.074575 
## 
## 
## ***** (V6) re75 *****
## before matching:
## mean treatment........ 2099.1 
## mean control.......... 1587.3 
## std mean diff......... 16.263 
## 
## mean raw eQQ diff..... 459.51 
## med  raw eQQ diff..... 192.2 
## max  raw eQQ diff..... 7119 
## 
## mean eCDF diff........ 0.062193 
## med  eCDF diff........ 0.049424 
## max  eCDF diff........ 0.13814 
## 
## var ratio (Tr/Co)..... 1.5231 
## T-test p-value........ 0.18694 
## KS Bootstrap p-value.. 0.16 
## KS Naive p-value...... 0.24713 
## KS Statistic.......... 0.13814 
## 
## 
## ***** (V7) u74 *****
## before matching:
## mean treatment........ 0.59459 
## mean control.......... 0.64815 
## std mean diff......... -10.858 
## 
## mean raw eQQ diff..... 0.055556 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 1 
## 
## mean eCDF diff........ 0.026777 
## med  eCDF diff........ 0.026777 
## max  eCDF diff........ 0.053554 
## 
## var ratio (Tr/Co)..... 1.0567 
## T-test p-value........ 0.41623 
## 
## 
## ***** (V8) u75 *****
## before matching:
## mean treatment........ 0.40541 
## mean control.......... 0.42593 
## std mean diff......... -4.1607 
## 
## mean raw eQQ diff..... 0.027778 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 1 
## 
## mean eCDF diff........ 0.01026 
## med  eCDF diff........ 0.01026 
## max  eCDF diff........ 0.020521 
## 
## var ratio (Tr/Co)..... 0.9856 
## T-test p-value........ 0.75936 
## 
## 
## ***** (V9) married *****
## before matching:
## mean treatment........ 0.18919 
## mean control.......... 0.16667 
## std mean diff......... 5.7246 
## 
## mean raw eQQ diff..... 0.018519 
## med  raw eQQ diff..... 0 
## max  raw eQQ diff..... 1 
## 
## mean eCDF diff........ 0.011261 
## med  eCDF diff........ 0.011261 
## max  eCDF diff........ 0.022523 
## 
## var ratio (Tr/Co)..... 1.1042 
## T-test p-value........ 0.66469 
## 
## 
## Before Matching Minimum p.value: 0.16 
## Variable Name(s): re75  Number(s): 6
```

```r
for(j in unique(lalonde$blocks.stepwise)[order(unique(lalonde$blocks.stepwise))])
 if(!is.na(j))
 {
   print(paste("Block",j,sep= " "))
   print(table.test(c("age","educ","black","nodegr","re74","re75","u74","u75","married"), "treat", lalonde[lalonde$blocks.stepwise==j&!is.na(lalonde$blocks.stepwise),]))
 }
```

```
## [1] "Block 1"
##         Mean control Mean treatment      t-stat Normalized diff
## age       24.7034483     25.1304348  0.42435910     0.060882464
## educ      10.0344828     10.1449275  0.79001308     0.116797411
## black      0.8137931      0.8115942 -0.03827226    -0.005605648
## nodegr     1.0000000      1.0000000         NaN             NaN
## re74    1444.6594414   1421.5610145 -0.03484484    -0.005084582
## re75     811.9289655    721.0537391 -0.19810286    -0.029682531
## u74        0.8482759      0.8695652  0.42062481     0.060865595
## u75        0.8827586      0.8985507  0.34800329     0.050354593
## married    0.1310345      0.1884058  1.04064293     0.156199073
## [1] "Block 2"
##         Mean control Mean treatment     t-stat Normalized diff
## age       25.6203704     25.9369369  0.3330083      0.04498706
## educ      10.1296296     10.2792793  0.4942682      0.06678055
## black      0.8611111      0.8558559 -0.1110559     -0.01500875
## nodegr     0.6203704      0.5585586 -0.9273706     -0.12532397
## re74    1925.2074352   2608.9502703  1.1041041      0.14896243
## re75    1587.2935426   2099.0802748  1.3240043      0.17869834
## u74        0.6481481      0.5945946 -0.8145360     -0.11007202
## u75        0.4259259      0.4054054 -0.3067021     -0.04145596
## married    0.1666667      0.1891892  0.4340460      0.05864583
```

Since there is a lot of information, you can aggregate the statistics above to perform a global analysis. See Chapter 13 of Imben's and Rubin's book.

## Trimming in order to improve overlap

As argued in Imben's and Rubin's book, it suffices to test for equality of means of (a monotonic transformation) of the propensity score to assess if covariate distributions are equal in the treatment and control group (Theorem 14.1). We do so by including the linearised propensity score in our test table:


```r
table.test(c("age","educ","black","nodegr","re74","re75","u74","u75","married","linear.ps.stepwise"), "treat", lalonde)
```

```
##                    Mean control Mean treatment     t-stat Normalized diff
## age                  25.0538462     25.8162162  1.1140361     0.107277121
## educ                 10.0884615     10.3459459  1.4421840     0.141219821
## black                 0.8269231      0.8432432  0.4577777     0.043886611
## nodegr                0.8346154      0.7081081 -3.1084981    -0.303986439
## re74               2107.0268154   2095.5740000 -0.0227466    -0.002159921
## re75               1266.9092408   1532.0556297  0.8692061     0.083863254
## u74                   0.7500000      0.7081081 -0.9746890    -0.094140477
## u75                   0.6846154      0.6000000 -1.8299743    -0.176809436
## married               0.1538462      0.1891892  0.9668363     0.093640701
## linear.ps.stepwise   -0.4289470     -0.2323553  4.4024141     0.432278474
```

From what we see above, there appears to be some, though arguably not much, imbalance in the estimated propensity score. Unbalance is undesirable, because it probably means our estimation will hinge on some unwanted level of extrapolation between individuals with different covariates, i.e. there is lack of overlap. In this case, one possibility is to *trim* the sample, so as to improve overlap. The approach in Imbens and Rubin's is based on *efficiency bound* calculations. Using semiparametric theory, we can calculate the best asymptotic variance a consistent estimator of an average treatment effect *in a subpopulation* can achieve. The idea of Imbens and Rubin is then to trim our sample so as to achieve the smallest variance (using an estimate of the bound). Note that in this approach, we change the estimand of interest, as we are looking at a more restrictive subpopulation.

We implement this approach in a function below. Given a linearized propensity score, we calculate the trimming constant and return a vector which indicates which observations should be kept. See Chapter 16 in Imbens and Rubin for details.


```r
trimming.imbens <- function(lin.ps,step.gamma.grid = 1e-5)
{
 inv.vec = 1/(plogis(lin.ps)*(1-plogis(lin.ps)))
 
 if(max(inv.vec)<=2*mean(inv.vec))
 {
   print("No trimming")
   return(rep(T,length(lin.ps))) 
  }else {
     gamma.grid = seq(min(inv.vec), max(inv.vec), by = step.gamma.grid)
     
     values = sapply(gamma.grid, function(gamma){
        (2*sum(inv.vec[inv.vec<=gamma])  - gamma*sum(inv.vec<=gamma))
     })
       
     values[values<0] = Inf
     
     gamma = gamma.grid[which.min(values)]
     
     alpha.trim = 1/2 - sqrt(1/4 - 1/gamma)
     print(paste("Trimming threshold alpha is ",alpha.trim))
     return(plogis(lin.ps)<= 1-alpha.trim &  plogis(lin.ps)>=alpha.trim)
   }
  
}
```


```r
lalonde$keep.trim.stepwise = trimming.imbens(lalonde$linear.ps.stepwise)
```

```
## [1] "Trimming threshold alpha is  0.13308191816384"
```

```r
lalonde.rest = lalonde[lalonde$keep.trim.stepwise,]

table.test(c("age","educ","black","nodegr","re74","re75","u74","u75","linear.ps.stepwise"), "treat", lalonde.rest)
```

```
##                    Mean control Mean treatment        t-stat Normalized diff
## age                  25.0538462     25.8369565  1.1411928382    1.100868e-01
## educ                 10.0884615     10.3152174  1.2845755594    1.258279e-01
## black                 0.8269231      0.8423913  0.4326711464    4.155127e-02
## nodegr                0.8346154      0.7119565 -3.0163914923   -2.954304e-01
## re74               2107.0268154   2106.9629891 -0.0001264475   -1.202582e-05
## re75               1266.9092408   1540.3820190  0.8939014756    8.640235e-02
## u74                   0.7500000      0.7065217 -1.0089443024   -9.762317e-02
## u75                   0.6846154      0.5978261 -1.8728429643   -1.812650e-01
## linear.ps.stepwise   -0.4289470     -0.2461534  4.2883854951    4.195512e-01
```

Trimming improved balance a bit. If, in your case, you think it is not enough, you may want to plot the linearised propensity score in the treatment and control groups (two histograms) to see if you can visually determine where overlap fails. The important thing here is that we are *never* looking at oucomes, so we do not introduce systematic biases in our analysis. In the end, one may have to live with a bit of extrapolation.

When we are interested in ATT, an alternative to trimming is to construct a matched sample, where each treated unit is matched (without replacement) to the ``closest'' (in terms of covariates) control units. This creates a sample with $2*N_{\text{treat}}$ units, which we can then use for estimating effects. I do not discuss this approach here; cf. Chapter 15 of the book for details.

## References

### Book

Imbens, G., & Rubin, D. (2015). Causal Inference for Statistics, Social, and Biomedical Sciences: An Introduction. Cambridge: Cambridge University Press. doi:10.1017/CBO9781139025751

### Articles
Imbens G. Matching Methods in Practice. Journal of Human Resources, 2015;50(2):373-419. 
