---
title: "Weak instrumental variables"
author: "Luis Alvarez"
date: "10/22/2020"
output: 
  html_document:
    keep_md: true
---



# Weak IV is a dividing by zero problem
In this session, we will discuss the problem of weak instrumental variables, and how to work with it. Our presentation will draw heavily on Andrews et al. (2019).

Consider the linear model from Econometrics 1:

$$
y_t  = \beta x_t + \gamma'r_t + \epsilon_t \\
x_t = \pi' z_t + \omega' r_t + u_t \\
\mathbb{E}[z_t \epsilon_t] = 0 \\
\mathbb{E}[z_t u_t] = 0 \\
\mathbb{E}[r_t \epsilon_t] = 0 \\
\mathbb{E}[r_t  u_t] = 0
$$
where $x_t$ is a *scalar* endogenous variables, and $z_t$ is a $k \times 1$ vector of instruments. The $r_t$ are exogenous controls, whereas  $\pi$ is a $k \times 1$ *vector* encompassing the partial association between $x_t$ and the $z_t$.  From Econometrics 1, we know that, if $\pi'$ has rank $1$ (i.e. if at least one entry of $\pi$ is different from 0), then $\beta$ is identified. Indeed we have that:

$$
y_t = \beta\pi' z_t + (\gamma +\beta\omega') r_t +  \beta u_t + \epsilon_t = \kappa'z_t + \psi'r_t + \xi_t \quad \text{(*)} \\
 \mathbb{E}[z_t \xi_t] = \mathbb{E}[z_t (\beta u_t + \epsilon_t)] = 0  \\
 \mathbb{E}[r_t \xi_t ] = \mathbb{E}[r_t (\beta u_t + \epsilon_t)]  = 0
$$ 

which implies that $\kappa := \pi \beta$ is identified from the linear projection of $y_t$ on $z_t$ and $r_t$. Since $\pi$ is identified from the linear projection of $x_t$ on $z_t$ and $r_t$; to show identification, it will be sufficient to express $\beta$ as a function of $\kappa$ and $\pi$. But, by the rank condition, we get:

$$\beta  =  (\pi' W \pi)^{-1}\pi' W \kappa \quad \text{(**)}$$ 

for any positive definite  $k \times k$ matrix $W$, which proves the desired identification result.

The above result shows $\beta$ is identified if $\pi$ is not zero. But what happens if $\pi$ is different than zero, but *close* to zero? Intuitively, we will be dividing by a number, $\pi' W \pi$, which is quite close to zero. If an estimator of $\beta$ has a form similar to (**) (as 2SLS has; more on this below), then we will be dividing by a random quantity "close" to zero. As you might imagine, this will lead to all sorts of problems. This is essentially the weak IV problem.

To get more intuition on the discussion above, let us strengthen the weak IV model above with a fixed regressor and Gaussian assumptions. In order not to overwhelm you with notation, we also consider a case without exogenous controls $r_t$ (you can then generalise to the case of controls by applying the Frisch-Waugh Lovell theorem to each stage of 2SLS).

$$
y_t  = \beta x_t  + \epsilon_t \\
x_t = \pi' z_t  + u_t \\ 
\{z_t\}_{t=1}^T \ \text{fixed} \\
\begin{pmatrix} 
\epsilon_t \\ 
u_t    
\end{pmatrix} \sim N\left(\begin{bmatrix}0 \\ 0 \end{bmatrix} , \begin{bmatrix} \sigma^2_1 & \rho \sigma_1 \sigma_2 \\ 
\rho \sigma_1 \sigma_2  & \sigma^2_2
\end{bmatrix}\right) \\
\{\epsilon_t, u_t\}_{t=1}^T \text{ iid}
$$
In this, case, we note that the 2SLS estimator of $\beta$ may be written as:

$$
\hat{\beta} = (X' P_Z X)^{-1}X' P_Z y =  (\hat{\pi}'Q_Z \hat{\pi})^{-1} \hat{\pi}' Q_Z \hat{\kappa} \quad \text{(***)}
$$
where $Z = \begin{bmatrix} z_1' \\ z_2' \\ \ldots \\ z_T'\end{bmatrix}$; $P_Z = Z(Z'Z)^{-1}Z'$; $Q_Z = (Z'Z) = \sum_{t=1}^T z_t  z_t'$; and $\hat{\pi}$ and $\hat{\kappa}$ are the OLS estimators of, respectively, the first-stage and the reduced-form equation (*). But, under the Gaussian linear model assumptions, we know from Econometrics I that:

$$
\begin{pmatrix}
\hat{\kappa} \\ \hat{\pi}
\end{pmatrix} \sim N\left( \begin{bmatrix} \kappa \\ \pi \end{bmatrix}, \begin{bmatrix} \sigma^{2}_1 (Z'Z)^{-1} & \rho \sigma_1 \sigma_2 (Z'Z)^{-1} \\ \rho \sigma_1 \sigma_2 (Z'Z)^{-1}  & \sigma^2_2  (Z'Z)^{-1}\end{bmatrix}  \right)
$$

This shows that, if $\pi$ is small, **relatively to the variance** of $\hat{\pi}$, $\sigma_2^2(Z'Z)^{-1}$, then the OLS estimator  $\hat{\pi}$ of $\pi$ will be very often close to $0$, meaning that the estimator in $\text{(***)}$ will be quite often divided by something close to 0. This will lead to: (1) serious biases in the estimator $\text{(***)}$ (recall that in Econometrics 1 we discuss *consistency* of 2SLS); (2) confidence intervals based on a Gaussian asymptotic approximation (recall that in Econometrics 1 we discuss *asymptotic normality* of 2SLSL) working quite poorly in finite samples. In other words, the asymptotic approximation of Econometrics 1 will serve as a poor guide to practice. In what follows, we will try to develop tools in order to detect this problem and, when present, circumvent it.

Before concluding this section, it is worth noting that the problem of Weak IV is a problem of $\hat{\pi}'Q_Z\hat{\pi}$ being quite often small. Is there a way to detect whether an instrument is weak? The answer is positive. Consider the test, in the first-stage equation, of $H_0: \pi = 0$ against the alternative of at least one $\pi$ different than zero. **Under homoskedasticity**, the F (Wald) test statistic of this null is (up to a degrees of freedom correction):

$$  \hat{\pi}'\widehat{\mathbb{V}(\hat{\pi})}^{-1} \hat{\pi} = \frac{1}{\hat{\sigma}^2_2}\hat{\pi}' Q_Z \hat{\pi} \quad \text{ (****)}$$
i.e. the F-test under homoskedasticity is a measure of how small is the denominator. This is the statistic behind most approaches for testing if instruments are weak. Importantly, though, since we will be testing a null of weak instruments (to be defined next), and **not** the null $\pi = 0$, critical values will **not** be those of a standard $F$-test of instrument **irrelevance**. In words, an instrument may be relevant and yet weak.

# Weak-IV asymptotics and detecting weak instruments

To produce an asymptotic theory that may account for the fact that, when instruments are weak, things do not go as well as the standard theory predicts -- and to use this theory to work out a solution to the problem --, Staiger and Stock (1997) introduced Weak-IV asymptotics. The idea is to work with the model

$$
y_t  = \beta x_t + \gamma'r_t + \epsilon_t \\
x_t = {\pi_{\color{red}T}}' z_t + \omega' r_t + u_t \\
\mathbb{E}[r_t \epsilon_t] = 0 \\
\mathbb{E}[r_t u_t] = 0 \\
\mathbb{E}[z_t u_t] = 0 \\
\mathbb{E}[z_t \epsilon_t] = 0
$$
where $\pi_T = C/\sqrt{T}$, $C$ being a fixed $s \times 1$ vector. The idea behind this asymptotic approximation is to capture the notion of $\pi$ being of the same order of the variance of $\hat{\pi}$, which was what drove the problem in the previous section.

To see the power of such a theory, let us consider a simpler model, with no covariates and a single instrument, i.e. 

$$
y_t  = \beta x_t + \epsilon_t \\
x_t =  \pi_T z_t  + u_t \\
\mathbb{E}[x_t]=\mathbb{E}[z_t] = 0 \\
\mathbb{E}[z_t u_t] = 0 \\
\mathbb{E}[z_t \epsilon_t] = 0\\
\{y_t, x_t, z_t\}_{t=1}^T \text{ iid }
$$
In this case, the OLS estimator is:

$$\tilde{\beta} =\frac{\sum_{t=1}^T x_t y_t }{\sum_{t=1}^T x_t^2}$$
whereas an IV estimator is:

$$\hat{\beta} =\frac{\sum_{t=1}^T z_t y_t }{\sum_{t=1}^T z_t x_t}$$

What is the behaviour of each estimator under Weak-IV asymptotics? The OLS estimator works similarly as in "strong IV asymptotics"

$$\tilde{\beta} =\frac{\sum_{t=1}^T x_t y_t }{\sum_{t=1}^T x_t^2} = \beta + \frac{T^{-1}\sum_{t=1}^T x_t \epsilon_t }{T^{-1}\sum_{t=1}^T x_t^2}  = \beta + \frac{\pi_T T^{-1}\sum_{t=1}^T \epsilon_t z_t  + T^{-1}\sum_{t=1}^T u_t \epsilon_t}{\pi_T^2 T^{-1}\sum_{t=1}^T z_t^2 + 2 \pi_T T^{-1} \sum_{t=1}^T z_T u_t + T^{-1} \sum_{t=1}^Tu_t^2} \overset{p}{\to} \beta + \frac{\mathbb{E}[u_t \epsilon_t]}{\sigma^2_u}$$
i.e. the OLS estimator is "asymptotically biased" if $$ \mathbb{E} [u_t \epsilon_t]\neq 0 .$$
As for the IV estimator, we get:

$$\hat{\beta} = \beta + \frac{\sum_{t=1}^T z_t \epsilon_t }{ \frac{C}{\sqrt{T}}\sum_{t=1}^T z_t^2 + \sum_{t=1}^T z_t u_t} = \beta + \frac{\frac{1}{\sqrt{T}}\sum_{t=1}^T z_t \epsilon_t}{C  T^{-1 }\sum_{t=1}^T z_t^2 + \frac{1}{\sqrt{T}} \sum_{t=1}^T z_t u_t} \overset{d}{\to} \beta + \frac{A}{C\sigma^2_z + B}$$
where, by the CLT, $A$ and $B$ are zero-mean jointly normal variables:

$$
\begin{pmatrix}
A \\
B
\end{pmatrix} \sim N\left( \begin{bmatrix}
0 \\
0
\end{bmatrix}, \begin{pmatrix}
\mathbb{E}[z_t^2 \epsilon_t^2 ] & \mathbb{E}[z_t^2 \epsilon_t u_t ] \\
\mathbb{E}[z_t^2 \epsilon_t u_t ] & \mathbb{E}[z_t^2  u_t^2 ] 
\end{pmatrix}\right)
$$
i.e. the IV estimator converges to  *a nonstandard random variable*. Notice that this random variable, in general, does not have mean $\beta$, meaning that IV **is asymptotically biased** under weak IV asymptotics. We can further compare the biases of both estimators by noting that, by the properties of the bivariate normal distribution:

$$\text{Bias}_{IV} = \mathbb{E}\left[ \frac{A}{C \sigma^2_z + B}\right] = \frac{\mathbb{E}[z_t^2 \epsilon_t u_t] }{\mathbb{E}[z_t^2 u_t^2]} \times \mathbb{E}\left[\frac{B}{C\sigma^2_z + B}\right] =  \frac{\mathbb{E}[z_t^2 \epsilon_t u_t] }{\mathbb{E}[z_t^2 u_t^2]}  - \frac{\mathbb{E}[z_t^2 \epsilon_t u_t]  }{\mathbb{E}[z_t^2 u_t^2]} \times \mathbb{E}\left[\frac{C\sigma^2_z}{C\sigma_z^2 + B}\right] $$

and, **assuming homoskedasticity** of errors, we further get $\frac{\mathbb{E}[z_t^2 \epsilon_t u_t] }{\mathbb{E}[z_t^2 u_t^2]}  = \frac{\mathbb{E}[\epsilon_tu_t]}{\sigma^2_u}$:

$$\text{Bias}_{IV} = \text{Bias}_{OLS}  - \frac{\mathbb{E}[\epsilon_tu_t]  }{\sigma^2_u} \times \mathbb{E}\left[\frac{C\sigma^2_Z}{C\sigma_z^2 + B}\right] $$

**Under the homoskedastic case**, Stock and Yogo (2005) show how the relation above can lead to a meaningful definition of weak instrument. They propose two alternative definitions of a weak instrument. In the first one, an instrument is deemed weak **if, in a worst-case scenario, the  bias of the 2SLS estimator exceeds the bias of the OLS estimator in more than x\%**. The worst case here is taken as a sup over $\beta$, which is not consistently estimable in our environment.

Yet a second definition of weak-instrument has to do with the size of a t-test on $\beta$. An instrument is deemed weak if, over all possible values of $\beta_0$, a test of the null $\beta = \beta_0$ at level $\alpha$ has worst-case size exceeding the nominal level by some amount.

Stock and Yogo then proceed to show that, in the *serially uncorrelated and homoskedastic case*, an F-test as the one discussed in the previous section can be used to test the null of an instrument being weak. Critical values were tabulated by the authors when different numbers of instruments are included. As discussed in Andrews et al. (2019):

> If we define the instruments as weak when the worst-case bias of two-stage least squares exceeds 10% of the worst-case bias of OLS, then the results of Stock & Yogo show that, for between 3 and 30 instruments, the appropriate critical value for a 5% test of the null of weak instruments ranges from 9 to 11.52 and so is always close to the Staiger & Stock (1997) rule-of-thumb cutoff of 10. By contrast, if we define the instruments as weak when the worst-case size of a nominal 5% t-test based on two-stage least squares exceeds 15%, then the critical value depends strongly on the number of instruments and is equal to 8.96 in cases with a single instrument but rises to 44.78 in cases with 30 instruments.

As argued before, the previous results hinge crucially on the homoskedasticity and serial uncorrelation assumptions. These allowed us to compare the biases of IV and OLS in a "clean" way. In a heteroskedastic and/or serially dependent case, one may be tempted to use a robust version of the Wald/F test, by replacing the homoskedastic variance estimator in $(****)$ by a robust estimate.  **However, there is no formal justification for that.** Intuitively, this approach is not correct, as we are indeed interested in $\hat{\pi}'Q_Z \hat{\pi}$, which is the quantity we want to be large. Replacing the middle matrix won't help.

Building on this intuition, Olea and Pflueger (2013) propose a procedure for detecting weak IVs where we first estimate the homoskedastic (serially uncorrelated) F-statistic, and then multiply the estimate by a factor which captures serial corelation and/or heteroskedasticity and/or clustering etc. They provide critical values based on the first definition of weak instrument of Stock and Yogo (2005). The critical values depend on an "effective critical value", to be computed from the data. As a rule of thumb for a 5\% test that the worst-case relative bias of 2SLS (with respect to OLS) exceeds 10\%, they suggest the critical value of 23.1 for their corrected statistic.

Let us prepare a code that implements the Olea and Pflueger (2013) correction:


```r
#Code to implement the Olea Pflueger correction:
#FS is
#x= Z\pi +  R \gamma + epsilon

#Arguments are 
# edndogenous = character variable with name of dependent variable
# instruments = character vector with name of instruments
# vcov = variance matrix to be used in computing first stage, e.g. vcovCL for cluster robust
# data = data.frame with data
# controls = vector with controls to be included in refression. c() if no controls.
# intercept = shoudl an intercept be included in formulas? Variable named intercept will be included among
#             controls
# weights = vector of weights, if weighted regression is desired
# cluster = if vcov = vcovCL, the name of the cluster variable (defaults to NULL)
# ... = additional arguments, to be passed to vcov function, e.g. degree of freedom correction
olea_pflueger_f <- function(endogenous, instruments, vcov, data, controls = c(), intercept = T, weights = NULL, cluster = NULL, ...)
{
  if(length(controls)>0)
    data.kept = data[,c(endogenous,instruments, controls)] else data.kept = data[,c(endogenous,instruments)]
    
    keep_ind = complete.cases(data.kept)
    
    data.kept = data.kept[keep_ind,]
    
    Nobs = nrow(data.kept)
    
    if(intercept)
    {
      data.kept = cbind(data.kept, "intercept" = 1)
      controls = c(controls, "intercept")
    }
    
    if(length(controls)>0)
    {
      # y = as.vector(residuals(lm(as.formula(paste(endogenous, "~ -1 + ", paste(controls,collapse = "+"),sep="")), data.kept)))
      
      Z = c()
      
      for(instrument in instruments)
      {
        z = as.vector(residuals(lm(as.formula(paste(instrument, "~ -1 + ", paste(controls,collapse = "+"),sep="")), data.kept)))
        
        Z = cbind(Z, z)
      }
      
    } else {
      # y = as.vector(data.kept[,endogenous])
      
      Z = as.matrix(data[,instruments])
      
    }
    
    formula.fs = paste(endogenous,"~ -1 +",paste(c(instruments,controls),collapse = " + "))
    
    if(is.null(weights))
      fs.reg =lm(as.formula(formula.fs), data.kept) else fs.reg =lm(as.formula(formula.fs), data.kept,
                                                                    weights = weights[keep.ind])
    
    # if(is.null(weights))
    #   fs.reg = lm(y~Z-1) else fs.reg = lm(y~Z-1, weights = weights[keep_ind])
    
    coefs = fs.reg$coefficients[names(fs.reg$coefficients)%in%instruments]
    
    if(!is.null(cluster))
     vcov_mat = vcov(fs.reg,cluster = data[keep_ind, cluster], ...) else vcov_mat = vcov(fs.reg, ...)
    
    #Restricting to only instruments
    vcov_mat = vcov_mat[names(fs.reg$coefficients)%in%instruments,names(fs.reg$coefficients)%in%instruments]
    
    if(is.null(weights))
      Q_Z_norm = t(Z)%*%Z/Nobs else Q_Z_norm = t(Z)%*%diag(weights[keep_ind])%*%Z/Nobs
    
    F_eff = t(coefs)%*%Q_Z_norm%*%coefs/sum(diag(vcov_mat%*%Q_Z_norm))
    
    
    return(list("Nobs" = Nobs, "Effective F" = F_eff))

}
```

Let us apply our function to test a simple example. We will use data from 48 US states in 1985 and 1995. Our goal is to estimate the price-elasticity of cigarette demand. As a simple instrument, we will use nonexcise taxes (which excludes sales taxes), which should affect prices through demand but not supply, as an instrument. We will cluster standard errors at the state level. 

The dataset is available on package __AER__, which provides an __ivreg__ function which computes robust F tests (when combined with package __sandwich__). Let us compare this robust F test -- which works as a test for __instrument relevance__ -- with our test for weak instruments constructed above.





```r
#Loading packages
#Package with IVreg function
library(AER)
```

```
## Loading required package: car
```

```
## Warning: package 'car' was built under R version 3.6.2
```

```
## Loading required package: carData
```

```
## Warning: package 'carData' was built under R version 3.6.2
```

```
## Loading required package: lmtest
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

```
## Loading required package: sandwich
```

```
## Loading required package: survival
```

```
## Warning: package 'survival' was built under R version 3.6.2
```

```r
#Robust SEs
library(sandwich)
#Hypothesis testing
library(lmtest)

## Stock and Watson (2007)
## data and transformations -- from AER package
data("CigarettesSW")
CigarettesSW$rprice <- with(CigarettesSW, price/cpi)
CigarettesSW$rincome <- with(CigarettesSW, income/population/cpi)
CigarettesSW$rtax <- with(CigarettesSW, tax/cpi)
CigarettesSW$rtdiff <- with(CigarettesSW, (taxs - tax)/cpi)
CigarettesSW$l_packs = log(CigarettesSW$packs)
CigarettesSW$l_rprice = log(CigarettesSW$rprice)

#IV regression
iv.model =  ivreg(l_packs ~ l_rprice | rtdiff, data = CigarettesSW)

summary(iv.model, diagnostics = T, vcov. = vcovCL, cluster = CigarettesSW$state)
```

```
## 
## Call:
## ivreg(formula = l_packs ~ l_rprice | rtdiff, data = CigarettesSW)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.64918 -0.08151  0.02625  0.08537  0.41205 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   9.9552     0.7458  13.349  < 2e-16 ***
## l_rprice     -1.1322     0.1598  -7.083 2.54e-10 ***
## 
## Diagnostic tests:
##                  df1 df2 statistic  p-value    
## Weak instruments   1  94    90.421 2.03e-15 ***
## Wu-Hausman         1  93     0.024    0.878    
## Sargan             0  NA        NA       NA    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1688 on 94 degrees of freedom
## Multiple R-Squared: 0.5251,	Adjusted R-squared:  0.52 
## Wald test: 50.17 on 1 and 94 DF,  p-value: 2.541e-10
```

```r
#Olea-Pflueger
olea_pflueger_f("l_rprice", "rtdiff", vcovCL, data = CigarettesSW, cluster = "state")
```

```
## $Nobs
## [1] 96
## 
## $`Effective F`
##          [,1]
## [1,] 125.8276
```
So we reject the null of weak instruments by their rule of thumb! Importantly, we once again stress that the p-value in the "Weak instruments" output from the __ivreg__ command is the p-value for the null of an **irrelevant** instrument, not a weak one (the description is wrong). The statistic reported is a robust F (Wald) test and p-values are computed from the corresponding F distribution.

**What about more than one endogenous variable?** The F-test is constructed in the single-endogenous regressor case. When more than one endogenous regressor is present, one can work with the Cragg-Donald test statistic, which was originally suggested for testing the rank of a matrix (the first-stage matrix of coefficients). This statistic specializes to the homoskedastic (serially uncorrelated) F test when a single endogenous regressor is present. Stock and Yogo (2005) tabulate the critical values for the test based on the Cragg-Donald statistic when more than one endogenous regressor is present. However, we note that the Olea and Pflueger correction was only derived for the single regressor case. How to proceed under heteroskedasticity when more than one endogenous regressor is included is still an open question.

# What to do when instruments are (are not) weak

When a researcher rejects the null of weak instruments, she may argue that there is evidence in the data that the Weak IV problem is not "important". Indeed, it is common practice for practitioners to scan the first-stage F (though usually not the correct one), and, if the F is high enough (usually $\geq 10$), they proceed by performing inference in the second stage as usual (t-tests, Normal CIs). This seems rather intuitive. However, it has been recently shown by Lee et al. (2020) that **pre-testing induces size distortions in the second-stage**. In words, if one reports the t-test in the second-stage only if the F-test exceeds some threshold (and reports something else, or nothing, if it is below the threshold), then the t-test from the specifications that survive this pre-test will be rejecting the null more often than the nominal level. They show that, for one to use the conventional critical values from the normal distribution in the second-stage, the F-test should be greater than 104.7!! They then propose a procedure that does not impose such a high bar, by constructing critical values for the second-stage as a function of the first-stage F. We do not pursue this approach here.

More importantly, what one should do if instruments **are weak**? In this case, one should **not** discard the specification!! Andrews et al. (2019) provide simulation evidence that this ONLY reinforces the problem pointed out by Lee et al. (2020). More importantly, there are tools to account for weak IV, so you do not need to discard your instrument if **it is relevant, but not strong**. Even more interestingly, these tools work **both** under Weak and Strong IVs, so you could always report confidence intervals robust to Weak Identification. The most common tool (Anderson-Rubin-based CIs) has even some optimality properties in just-identified models (number of endogenous regressors = number of instruments), being optimal when identification is weak and also when it is strong (no loss of power vis-à-vis t-test in the latter). It may, however, have low power when the model is overidentified (more instruments than endogenous variables). In this case, however, other, well-powered approaches exist (the CLR test of Moreira (2003) being the leading one). Perhaps surprisingly, then, **always** reporting weak-IV-robust CIs does not appear to be common practice. We discuss the weak-IV-robust approach below.

# Anderson-Rubin CIs

We focus on the single endogenous regressor case, though the discussion generalises to more regressors. Recall from Section 1 that the IV model imposes the restriction that:

$$ \kappa = \pi \beta $$
where $\kappa$ is the reduced-form coefficient in $(*)$, $\pi$ is the first-stage relation and $\beta$ is the structural coefficient of interest. The Anderson-Rubin approach is based on the fact that, under the null $H_0: \beta = \beta_0$, the Wald statistic:

\begin{equation}
(\hat{\kappa} - \hat{\pi} \beta_0)'\widehat{\mathbb{V}(\hat{\kappa} - \hat{\pi} \beta_0)}^{-1}(\hat{\kappa} - \hat{\pi} \beta_0)
\end{equation}
has limiting distribution equal to a Chi-squared random variable with 1 degree of freedom, under __both Weak and Strong IV asymptotics__. The variance of $\hat{\kappa} - \hat{\pi} \beta_0$ can be found by "pooling" together  both the reduced-form and first-stage regressions to compute $\mathbb{V}\left[\begin{pmatrix}\hat{\kappa}  \\ 
\hat{\pi} \end{pmatrix}\right]$, and then using the formula for the variance of linear combinations of random variables.

Let us introduce a code that computes the Anderson-Rubin test statistic for the null $H_0: \beta = \beta_0$. We will allow for a single endogenous regressor only, but more than one instrument:



```r
#AR test

#Arguments are 
# outcome = character variable with name of outcome variable of interest
# edndogenous = character variable with name of endogenous variable
# instruments = character vector with name of instruments
# vcov = variance matrix to be used in pooled model.
# data = data.frame with data
# beta_0 = value under the null (defaults to 0)
# controls = vector with controls to be included in refression. c() if no controls.
# intercept = shoudl an intercept be included in formulas? Variable named intercept will be included among
#             controls
# weights = vector of weights, if weighted regression is desired
# cluster = if vcov = vcovCL, the name of the cluster variable (defaults to NULL)
# ... = additional arguments, to be passed to vcov function, e.g. degree of freedom correction

anderson_rubin_test <- function(outcome, endogenous, instruments, vcov, data, beta_0=0, controls = c(), intercept = T, weights = NULL, cluster = NULL, ...)
{
    
  if(length(controls)>0)
    data.kept = data[,c(outcome, endogenous,instruments, controls)] else data.kept = data[,c(outcome, endogenous,instruments)]
    
    keep_ind = complete.cases(data.kept)
    
    data.kept = data.kept[keep_ind,]
    
    Nobs = nrow(data.kept)
    
    if(intercept)
    {
      data.kept = cbind(data.kept, "intercept" = 1)
      controls = c(controls, "intercept")
     }
    
    #We will pool outcome and endogenous now
    data.pooled = rbind(data.kept, data.kept)
    
    data.pooled$pool_variable = c(data.kept[,outcome], data.kept[,endogenous])
    data.pooled$variable_indicator = as.factor(c(rep("reduced_form",    Nobs),rep("first_stage",Nobs)))
    
   
    
    #Constructing the formula for regression
    if(length(controls)>0)
     formula = paste("pool_variable ~ -1 + ", paste(paste("variable_indicator", instruments, sep = ":"),collapse = "+"), "+", paste(paste("variable_indicator", controls, sep = ":"))) else  formula = paste("pool_variable ~ -1 +", paste(paste("variable_indicator", instruments, sep = ":"),collapse = "+")) 
       
       
    if(is.null(weights))
      pool.model = lm(formula, data.pooled) else  pool.model = lm(formula, data.pooled, weights = weights[rep(keep_ind,2)]) 
    
    coefs = pool.model$coefficients
    
     if(!is.null(cluster))
      vcov_model = vcov(pool.model,cluster = rep(data[keep_ind, cluster],2), ...) else vcov_model = vcov(pool.model, ...)

    
    lin_vec = 1*grepl(paste("reduced_form", instruments, sep = ":"), names(coefs)) - 
      beta_0*grepl(paste("first_stage", instruments, sep = ":"), names(coefs))
    
    #constructing test statistic
    val = (coefs%*%lin_vec)
    vcov_lin = t(lin_vec)%*%vcov_model%*%lin_vec
    
    ar = val%*%solve(vcov_lin)%*%val
    
    pvalue =1 - pchisq(ar, 1)
    
    return(list("AR test statistic" = ar, "P-value" = pvalue, "Nobs" = Nobs))
}
```


Let's use the AR test to test the null that the elasticity is zero.


```r
anderson_rubin_test("l_packs", "l_rprice", "rtdiff", vcovCL, CigarettesSW, cluster = "state")
```

```
## $`AR test statistic`
##          [,1]
## [1,] 29.13659
## 
## $`P-value`
##              [,1]
## [1,] 6.745101e-08
## 
## $Nobs
## [1] 96
```

So we reject the null that the elasticity is zero, right?

How can we construct confidence-intervals using the above testing procedure? The idea is to perform __test inversion__ . In particular, a $95\%$ confidence interval can be constructed by finding all values $\beta_0$ for which the AR test does _not_ reject the null at $5\%$. This amounts to performing grid search over a wide range of values of $\beta$. As argued in Andrews et. al, confidence sets based on inverting the AR test can sometimes be of the form $(-\infty, b]$, $[a,b] \cup [c,d]$ etc. They may also be empty. So one should be careful when determining the grids so we get these peculiarities.

The function below computes AR confidence sets by grid search over a grid supplied by the user. 


```r
#AR CI

#Arguments are 
# outcome = character variable with name of outcome variable of interest
# edndogenous = character variable with name of endogenous variable
# instruments = character vector with name of instruments
# vcov = variance matrix to be used in pooled model.
# data = data.frame with data
# grid_beta = grid over which to perform grid search
# confidence = confidence level for CI (defaults to 0.95)
# controls = vector with controls to be included in refression. c() if no controls.
# intercept = shoudl an intercept be included in formulas? Variable named intercept will be included among
#             controls
# weights = vector of weights, if weighted regression is desired
# cluster = if vcov = vcovCL, the name of the cluster variable (defaults to NULL)
# ... = additional arguments, to be passed to vcov function, e.g. degree of freedom correction

anderson_rubin_ci <- function(outcome, endogenous, instruments, vcov, data, grid_beta, confidence = 0.95, controls = c(), intercept = T, weights = NULL, cluster = NULL, ...)
{
    
  if(length(controls)>0)
    data.kept = data[,c(outcome, endogenous,instruments, controls)] else data.kept = data[,c(outcome, endogenous,instruments)]
    
    keep_ind = complete.cases(data.kept)
    
    data.kept = data.kept[keep_ind,]
    
    Nobs = nrow(data.kept)
    
    if(intercept)
    {
      data.kept = cbind(data.kept, "intercept" = 1)
      controls = c(controls, "intercept")
     }
    
    #We will pool outcome and endogenous now
    data.pooled = rbind(data.kept, data.kept)
    
    data.pooled$pool_variable = c(data.kept[,outcome], data.kept[,endogenous])
    data.pooled$variable_indicator = as.factor(c(rep("reduced_form",    Nobs),rep("first_stage",Nobs)))
    
   
    
    #Constructing the formula for regression
    if(length(controls)>0)
     formula = paste("pool_variable ~ -1 + ", paste(paste("variable_indicator", instruments, sep = ":"),collapse = "+"), "+", paste(paste("variable_indicator", controls, sep = ":"))) else  formula = paste("pool_variable ~ -1 +", paste(paste("variable_indicator", instruments, sep = ":"),collapse = "+")) 
       
       
    if(is.null(weights))
      pool.model = lm(formula, data.pooled) else  pool.model = lm(formula, data.pooled, weights = weights[rep(keep_ind,2)]) 
    
    coefs = pool.model$coefficients
    
     if(!is.null(cluster))
      vcov_model = vcov(pool.model,cluster = rep(data[keep_ind, cluster],2), ...) else vcov_model = vcov(pool.model, ...)

    p1 = grepl(paste("reduced_form", instruments, sep = ":"), names(coefs))
    p2 = grepl(paste("first_stage", instruments, sep = ":"), names(coefs))
    
    acc_vec = c() 
    
    #Looping over grid
    for(beta in grid_beta)
    {
      
        lin_vec = p1 - beta*p2
        #constructing test statistic
        val = (coefs%*%lin_vec)
        vcov_lin = t(lin_vec)%*%vcov_model%*%lin_vec
        
        ar = val%*%solve(vcov_lin)%*%val
        
        pvalue = pchisq(ar, 1)
        
        if(pvalue<= confidence)
          acc_vec = c(acc_vec, T) else acc_vec = c(acc_vec, F)
            
    }
    
    if(sum(acc_vec)==0)
      return("Confidence set is empty!") else {
        
        vec_region_start = c()
        vec_region_end = c()
        if(acc_vec[1]==T)
        {
          warning("Lower boundary point was accepted! Perhaps decrease grid lower bound to see what happens?")
          vec_region_start = grid_beta[1]
        }
        
        if(acc_vec[length(acc_vec)]==T)
        {
          warning("Upper boundary point was accepted! Perhaps increase grid upper bound to see what happens?")
          vec_region_end = grid_beta[length(acc_vec)]
        }
        
        vec_region_start = c(vec_region_start, grid_beta[c(F,diff(acc_vec)==1)]  )
        vec_region_end = c(grid_beta[c(diff(acc_vec)==-1, F)],vec_region_end)
        
        CI.text = paste(paste("[",vec_region_start, ",", vec_region_end, "]"),collapse = " U ")
        
        return(CI.text)
      }

    
}
```

Let's compute a 95\% CI in our setting.



```r
grid_beta = seq(-100, 100, 0.01)
anderson_rubin_ci("l_packs", "l_rprice", "rtdiff", vcovCL, CigarettesSW, cluster = "state", grid_beta = grid_beta, confidence = 0.95)
```

```
## [1] "[ -1.53 , -0.75 ]"
```

How do these CIs compare to those constructed using normal critical values?


```r
coefci(iv.model, vcov. = vcovCL, cluster= CigarettesSW$state, level = 0.95)
```

```
##                 2.5 %     97.5 %
## (Intercept)  8.140979 11.7694448
## l_rprice    -1.520507 -0.7439442
```


We note the inversion procedure readily generalises when more than one endogenous regressor is included, though it becomes increasingly computationally demanding to perform test inversion, because we have to perform grid search over all __possible combinations of parameters__.

Finally, as we have argued, the AR approach is optimal in just-identified models, but suffers from power losses when more instruments than endogenous regressors are included. In the latter, the conditional likelihood ratio (CLR) approach of Moreira (2003) and its generalisations to nonhomoskedastic cases ought to be preferred. See Andrews et al. (2019) for a discussion on how to implement these.

# Recommendations

Suppose you have a relevant instrument (i.e. one that rejects the null of $\pi=0$ using a (possibly robust) F-test). You could follow three approaches:

1. Always report Anderson-Rubin-based CIs (when model is just-identified) or CLR-based CIs (when there are more instruments than endogenous variables).

2. **Two-step approach:** perform a weak IV test (Stock and Yogo in homoskedastic non-autocorrelated case; Olea-Montiel more generally). If you reject the null of a weak IV (if the (adjusted) $F$ exceeds the corresponding critical value), proceed by using standard inference in the second stage. If you do not reject the null, use the approach in the previous item. This approach is subject to Lee's critique of size distortions, though with a sample of AER papers, Andrews et al. (2019) argue the distortion is not that large.

3. Perform a (possibly robust) F-test (no need for Olea and Montiel approach here). If the F-test statistic exceeds the 104.7 threshold in Lee et al. (2020), use t-stat-based inference with conventional normal critical values. If the F-test is smaller than that, use Lee's corrected critical values.

Finally, it seems to be always good to compare your normal CIs with those from a robust-to-weak-IV procedure. We also note that, when an instrument is weak, biases in **point estimators will be large**, as no consistent estimator is known to exist -- though there are some alternative approaches that may lead to lower biases than 2SLS. When instruments are weak, you should thus focus mostly on interval estimation (i.e. computing CIs).

# References

Andrews, I., Stock, J. H., & Sun, L. (2019). Weak instruments in instrumental variables regression: Theory and practice. Annual Review of Economics, 11, 727-753.

Lee, D. L., McCrary, J., Moreira, M. J., & Porter, J. (2020). Valid t-ratio Inference for IV. arXiv preprint arXiv:2010.05058.

Moreira, M.J. (2003), A Conditional Likelihood Ratio Test for Structural Models. Econometrica, 71: 1027-1048. doi:10.1111/1468-0262.00438

Olea, J. L. M., & Pflueger, C. (2013). A robust test for weak instruments. Journal of Business & Economic Statistics, 31(3), 358-369.

Staiger, D., & Stock, J. (1997). Instrumental Variables Regression with Weak Instruments. Econometrica, 65(3), 557-586. doi:10.2307/2171753

Stock, J. H., & Yogo, M. (2005). Testing for Weak Instruments in Linear IV Regression. Identification and Inference for Econometric Models: Essays in Honor of Thomas Rothenberg, 80.


