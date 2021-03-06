---
title: "Regression Discontinuity Designs"
author: "Luis Alvarez"
date: "10/29/2020"
output: 
  html_document:
    keep_md: true
---



## RD design

Let $Y_0$ and $Y_1$ denote potential outcomes of a binary treatment $D$ in a population of interest, which is assigned according to a deterministic rule $D = 1\{T\geq 0 \}$, where $T$ is a __continuous running variable__. The observed outcome is $Y = D Y_1 + (1-D) Y_0$. Let  $\mu_0(t)=  \mathbb{E}[Y_0|T =t ]$, $\mu_{1}(t) = \mathbb{E}[Y_1|T=t]$, and $\mu(t) := \mathbb{E}[Y|T = t]$, where expectations are taken wrt the population distribution. You have seen in class that, under the continuity assumption that $\mu_0$ and $\mu_1$ are continuous at $t=0$ (what does this assumption mean?), the conditional average treatment effect at the cutoff, $\tau : = \mathbb{E}[Y_1 = Y_0|T=0]$ is __identified__ as:

$$\tau = \mathbb{E}[Y_1|T=0] - \mathbb{E}[Y_0|T=0] =  \lim_{t \downarrow 0}\mathbb{E}[Y_1|T=t] - \lim_{t \uparrow 0}\mathbb{E}[Y_0|T=t] = \lim_{t \downarrow 0}\mathbb{E}[Y|T=t] - \lim_{t \uparrow 0}\mathbb{E}[Y|T=t] = \\ = \lim_{t \downarrow 0 } \mu(t) - \lim_{t \uparrow 0 } \mu(t)  =: \mu^+(0) - \mu^-(0)$$

Estimation of the treatment effect at the cutoff thus requires estimation of a left- (right-) limit of the conditional expectation of the outcome on the running variable. In principle, this could be done by fitting  flexible polynomials separately to data on the left (right) of the cutoff via (global) linear regression, and then predicting the outcome at the cutoff using each specification. However, since this is operated via linear regression, it is subject to extrapolation biases, as information far away from the cutoff could be very informative in prediction at the cutoff (Imbens and Gelman, 2019). Therefore, the preferred approach resorts to __local regression__, which we briefly discuss below.

## Local nonparametric regression
Consider a random sample $\{(Z_i,X_i)\}_{i=1}^n$ from a distribution $F_{Z,X}$, where $Z$ and $X$ are scalar random variables and, in addition, $X$ is a continuous random variable. Estimation of the condtional expectation $\mu(x) : =\mathbb{E}[Z|X=x]$ could in principle be performed using the __Nadaraya-Watson estimator__, i.e.:

$$\widehat{\mu(x)} = \frac{\sum_{i=1}^n K\left(\frac{X_i - x}{h}\right) Z_i}{\sum_{i=1}^n K\left(\frac{X_i - x}{h}\right)}$$
where $K$ is a __kernel function__ and $h > 0$ is a __bandwidth__ parameter, to be defined by the user. The kernel measures the "distance" -- rescaled by h -- between the realizations of the $X_i$ and the point of interest $x$. Typical choices of kernel include: (1) the uniform kernel $K(s) = \frac{1}{2}\{ - 1 \leq  s \leq 1 \}$; (2) the triangular kernel $K(s) = 1\{|s|\leq 1\} (1 - |s|)$; (3) the Epanechnikov kernel $K(s) =  1\{|s|\leq 1\}  \frac{3}{4}(1-s^2)$; and (4) the Gaussian kernel $K(s) = \frac{1}{\sqrt{2\pi}} \exp\left(-\frac{1}{2}s^2\right)$. The first three are the most common in the RD design, as they effectively discard observations far away from the point of interest. Let us plot these below:


```r
uniform <- function(s){dunif(s, min = -1,max = 1)}
triangular <- function(s){(abs(s)<=1)*(1-abs(s))}
epanechnikov <- function(s){(abs(s)<=1)*(3/4)*(1-s^2)}
gaussian <- function(s){dnorm(s)}
grid_s = seq(-2, 2, 1e-3)

#Plot
plot(grid_s, uniform(grid_s), ylab = "density", xlab = "value", type = "l",ylim = c(0,1))

lines(grid_s, triangular(grid_s), lty =2 )
lines(grid_s, epanechnikov(grid_s), lty =3, col =  "blue" )

lines(grid_s, gaussian(grid_s), lty =4, col = "red" )

legend('topleft', c("Uniform", "Triangular", "Epanechnikov", "Gaussian"), 
       col = c("black", "black", "blue", "red"), lty = 1:4, bty = "n")
```

![](rdd_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

From the plot above, it is clear that the Nadaraya-Watson estimator is a weighted average of the observations $Z_i$, where the weights are given by "closeness" of the corresponding $X_i$ to the point $x$ of interest. The bandwidth allows us to vary the penalization of this closeness: a small bandwidth penalizes deviations from $x$ more strongly than a higher $h$. Clearly, there is a bias-variance tradeoff in the choice of bandwidth: a smaller bandwidth will lead to low bias but high variance, as less observations are "in practice" used in estimation. A higher bandwidth has smaller variance, but higher bias.

How does one choose a bandwidth, then? The usual approach relies on choosing $h$ so as to minimize the mean-squared error:

$$
\mathbb{E}\left[\left(\widehat{\mu(x)} - \mu(x)\right)^2\right] = \text{Bias}\left(\widehat{\mu(x)}\right)^2  + \mathbb{V}\left[\widehat{\mu(x)}\right]
$$
(note that $x$ is __fixed__ here: the expectation is taken wrt to the $Z_i, X_i$, this is the error of estimating $\mu(x)$ with a random sample). The problem is that the mean-squared error is __unknown__. One approach thus relies on asymptotic expansions of the bias and variance, which allows one to estimate an approximate mean-squared error to be minimized by $h$. Yet another approach uses leave-one-out cross-validation to estimate the mean-squared error -- and chooses $h$ so as to minimize the cross-validated approximation. In the Nadaraya-Watson context with iid observations, it can be shown that _asymptotically_ both approximations lead to the same selection rule.

The problem with Nadaraya-Watson regression is that it works poorly (slow convergence rate) when the point $x$ is at the boundary of the support of $X$. This is the case in RDD, where we estimate effects separately to the right (left) of the cutoff. Thankfully, there exists a similar estimation approach which has the __same convergence rate__ irrespective of where we are at the support (interior or boundary). To introduce it, we first note that the Nadaraya-Estimation can be obtained as the solution to the minimisation:

$$
\hat{\mu}(x) \in \operatorname{argmin}_{a \in \mathbb{R}} \sum_{i=1}^n K\left(\frac{X_i - x}{h}\right) (Z_i - a)^2 \quad (*)
$$
for which reason the Nadaraya-Watson estimator is also known as __local constant regression__. Suppose now that the conditional expectation $\mu(s)$ is $k$-times differentiable at $s = x$. Then it admits the Taylor expansion:

\begin{equation}
\mu(s)  = \mu(x) +\sum_{j=1}^k \mu^{(j)}(x)  \frac{(s-x)^j}{j!}  +  o(|s-x|^k)  \quad (**)
\end{equation}
where $\mu^{(j)}(x)$ is the j-th derivative of $\mu$ at $x$. In light of $(*)$ and $(**)$, one is motivated to estimate $\mu(x)$ as the _intercept_ $\hat{b}_0$ of the __k-th order local polynomial__ regression.

$$(\hat{b}_0, \hat{b}_1, \ldots \hat{b}_k) \in \operatorname{argmin}_{b_0, b_1 \ldots b_k} \sum_{i=1}^n K\left(\frac{X_i - x}{h}\right) \left(Z_i - b_0 - \sum_{j=1}^k b_j \frac{(Z_i -x)^j}{j!}\right)^2$$
and the other $\hat{b}_j$ are estimates of derivatives of $\mu$ at $x$.

## Back to RDD
As seen in class, estimation in the RDD context is performed using local linear regression ($k=1$ above) __at each side of the cutoff__. Consider the environment in the first section. In particular let $\{(Y_i, D_i, T_i)\}_{i=1}^n$ be a random sample from the population of interest. One then estimates the quantities $\mu^+$ and $\mu_0^-$ as the intercepts of:

$$\hat{\mu}^+ = \hat{b}_0 \in \operatorname{argmin}_{b_0, b_1 }  \sum_{i=1}^n D_i \cdot  K\left(\frac{T_i - 0}{h^+}\right) \cdot \left(Y_i - b_0 - b_1 \cdot T_i \right)$$

$$\hat{\mu}^- = \hat{b}_0 \in \operatorname{argmin}_{b_0, b_1 }  \sum_{i=1}^n (1-D_i) \cdot  K\left(\frac{T_i - 0}{h^-}\right) \cdot \left(Y_i - b_0 - b_1 \cdot T_i \right)$$
and estimates the treatment effect by the difference $\hat{\tau}(h^+, h^-) = \hat{\mu}^+ - \hat{\mu}^-$. The bandwidth can be chosen differently on the right and on the left, though it is usually constrained to be equal. We assume this throughout ($h^+ = h^- = h_n$ ; the $n$ subscript is to emphasise that the bandwidth is a function of the sample size) and write $\hat{\tau}(h_n)$ for the estimator.

As argued in Calonico et al. (2014), under some technical conditions, it can be shown that, if  <span style="color:red"> $nh_n \to \infty$ and  $nh^5_n \to 0$ </span>,

$$\frac{\left(\hat{\tau}(h_n) - \tau\right)}{\mathbb{V}[\hat{\tau}(h_n)|\mathbf{T}]^{\frac{1}{2}}} \overset{d}{\to} N(0,1)$$
where $\mathbf{T} =(T_1,T_2\ldots T_n)$. This approximation could in principle be used as a basis for an inferential procedure on the treatment effect -- in particular, we note that  $\mathbb{V}[\hat{\tau}(h_n)|\mathbf{T}]^{1/2}$ could be estimated from the standard error output of a weighted regression using the kernel weights.  The problem is that the leading approaches to choose $h_n$ -- which, as we discussed, work by minimizing approximations of the mean squared error -- lead to bandwidth choices (as a function of sample size) "too large" to satisfy $nh^5_n \to 0$. For example, the approach in  Imbens and Kalyanaraman (2012), which consists of minimising an asymptotic expansion of:

$$\mathbb{E}[(\hat{\tau}(h) - \tau)^2] = \text{Bias}(\hat{\tau}(h))^2 + \mathbb{V}[\hat{\tau}(h)]$$
leads to a choice of bandwidth $h_n = C \cdot n^{-1/5}$ where $C > 0$ is a constant, implying that $nh^5_n \nrightarrow 0$. In this case, it can be shown that the RDD estimator __is still consistent__, but that now it satisfies the asymptotic distribution:

$$\frac{\left(\hat{\tau}(h_n) - B(h_n) -  \tau\right)}{\mathbb{V}[\hat{\tau}(h_n)|\mathbf{T}]^{\frac{1}{2}}} \overset{d}{\to} N(0,1)$$
where $B(h_n)$ is a __bias term__. This means the distribution of the estimator is no longer asymptotically centered around $\tau$ -- and so inference based on normal CIs will be misleading.

The usual approach in nonparametric statistics to circumvent the noncentering problem is to __bias-correct__ $\hat{\tau}(h_n)$. In words, one works out an asymptotic approximation to $B(h_n)$; and then estimates this approximation. This will involve estimation of the right- and left- second-order derivative of the conditional expectation functions $\mu$, which will be performed by local regression, and thus choosing another bandwidth -- known as a pilot bandwidth. Let $\hat{B}(h_n; b_n)$ denote the estimated bias, where $b_n$ is the pilot. It is shown in Calonico et al. that if $nh_n\to \infty$, <span style="color:red">  $nh^7_n \to 0$  and $h_n/b_n \to 0$</span>, then the bias corrected estimator satisfies:

$$
\frac{\left(\hat{\tau}(h_n) - \hat{B}(h_n, b_n) -  \tau\right)}{\mathbb{V}[\hat{\tau}(h_n)|\mathbf{T}]^{\frac{1}{2}}} \overset{d}{\to} N(0,1)
$$

as seen above, bias-correction ameliorates the rate concern by only requiring $nh^7 \to 0$. This requirement is satisfied by the MSE optimal rule of Imbens and Kalyanaraman (2011). Moreover, we note that, if $h_n/b_n \to 0$, then bias-correction does not affect the asymptotic variance, as seen in the formula above.  Once again, however, $b_n$ is often chosen via a MSE-optimal criterion, which will lead to bandwidth choices __not__ satisfying $h_n/b_n \to 0$. This will mean that the asymptotic distribution will be affected by the estimation error of the bias, and that inference ignoring this fact will be misleading.

In light of the above, Calonico et al. (2014) develop the asymptotic theory of the bias-corrected RDD estimator dropping the requirement that $h_n/b_n \to 0$. In this case, they have to correct the variance by estimation of the bias. They also propose a MSE-expansion-optimal selection rule for $h_n$ and $b_n$, where $h_n$ is chosen so as to minimize the asymptotic MSE of $\hat{\tau}$; and $b_n$ is chosen so that bias estimation is MSE-optimal. They similarly derive bandwidth rules and inferential methods for the fuzzy RDD design. These methods are implemented by the __rdrobust__ package in R. See https://rdpackages.github.io/rdrobust/ for further reference.

As an example of the methods available. We will estimate the effect of party incumbency in US senate elections on subsequent voting in the same party. The running variable is the margin of victory of the Democratic candidate in an election (negative if margin of loss); the outcome is the fraction of votes of the Democractic candidate in the next election. This dataset is available in the __rdrobust__ package:



```r
#loading package
library(rdrobust)
```

```
## Warning: package 'rdrobust' was built under R version 3.6.2
```

```r
data(rdrobust_RDsenate)

head(rdrobust_RDsenate)
```

```
##        margin     vote
## 1  -7.6885610 36.09757
## 2  -3.9237082 45.46875
## 3  -6.8686604 45.59821
## 4 -27.6680565 48.47606
## 5  -8.2569685 51.74687
## 6   0.7324815 39.80264
```

```r
#Plotting regression discontinuity plot
rdplot(y = rdrobust_RDsenate$vote, x = rdrobust_RDsenate$margin )
```

![](rdd_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

The RD plot above shows the average of the $Y$-axis variable by "bins" of the $X$-axis variable. The number and position of the bins is computed so as to minimise the distance between the conditional variance of the outcome, integrated along the left- and right-supports, and the (integrated) conditional variance of the estimated bins. The idea is that the bins should mimick the variance of the underlying data.  An alternative would be to choose the bins so as to minimise the integrated mean squared error (IMSE) of predicting the conditional expectation along the right- and left- sides of the support (this is implemented by option _binselect_ = "es" in the function). A fourth-order (global) polynomial (the order can be changed by argument p) is also fit separately to the left and to the right of the cutoff. Be careful when interpreting these global polymonials (Imbens and Gelman, 2019)!!!

Next, we estimate the effect of interest:


```r
tri = rdrobust(y = rdrobust_RDsenate$vote, x = rdrobust_RDsenate$margin, kernel =  "tri",all= T)
uni = rdrobust(y = rdrobust_RDsenate$vote, x = rdrobust_RDsenate$margin, kernel =  "uniform", all = T)

summary(tri)
```

```
## Call: rdrobust
## 
## Number of Obs.                 1297
## BW type                       mserd
## Kernel                   Triangular
## VCE method                       NN
## 
## Number of Obs.                 595         702
## Eff. Number of Obs.            360         323
## Order est. (p)                   1           1
## Order bias  (q)                  2           2
## BW est. (h)                 17.754      17.754
## BW bias (b)                 28.028      28.028
## rho (h/b)                    0.633       0.633
## Unique Obs.                    595         665
## 
## =============================================================================
##         Method     Coef. Std. Err.         z     P>|z|      [ 95% C.I. ]       
## =============================================================================
##   Conventional     7.414     1.459     5.083     0.000     [4.555 , 10.273]    
## Bias-Corrected     7.507     1.459     5.146     0.000     [4.647 , 10.366]    
##         Robust     7.507     1.741     4.311     0.000     [4.094 , 10.919]    
## =============================================================================
```

```r
summary(uni)
```

```
## Call: rdrobust
## 
## Number of Obs.                 1297
## BW type                       mserd
## Kernel                      Uniform
## VCE method                       NN
## 
## Number of Obs.                 595         702
## Eff. Number of Obs.            271         235
## Order est. (p)                   1           1
## Order bias  (q)                  2           2
## BW est. (h)                 11.597      11.597
## BW bias (b)                 22.944      22.944
## rho (h/b)                    0.505       0.505
## Unique Obs.                    595         665
## 
## =============================================================================
##         Method     Coef. Std. Err.         z     P>|z|      [ 95% C.I. ]       
## =============================================================================
##   Conventional     7.202     1.613     4.466     0.000     [4.041 , 10.364]    
## Bias-Corrected     7.593     1.613     4.708     0.000     [4.432 , 10.755]    
##         Robust     7.593     1.852     4.100     0.000     [3.963 , 11.224]    
## =============================================================================
```

The three methods are, in order: the standard method (no correction; Conventional); bias-corrected method (Bias-Corrected); and bias-corrected method accounting for variance of bias correction (Robust). The choice of bandwidth(s) was MSE optimal according to the rule in Calonico et. al (2014).

## Manipulation testing

As you have seen in class, an important credibility check for the validity of the RD design in estimating a causal effect is that there is no manipulation of the running variable. A possible test for manipulation would be to check if the density of the running variable to the left and to the right of the cutoff is different.  Intuitively, if there is more mass in either side, this could be evidence that individuals manipulate the cutoff so as to receive (not receive) treatment. Formally, we would like to test the null that:

$$H_0: \lim_{t \uparrow 0} f_T(t) = \lim_{t \downarrow 0} f_T(t)$$

where $f_T$ is the density of $T$; against an alternative that the densities are different. In an influential paper, McCrary (2008) originally suggested testing this null by estimating the density at either side of the cutoff using the _intercept_ of a local linear regression on an histogram estimator of the pdf -- and then performing a t-test on the difference. This approach involved, along with the bandwidth selection of the local linear regression, selecting the binpoints for estimation of the histogram. More recently, Cattaneo et al. (2020) suggested estimating the density at either side by the coefficient $b_1$ of the local polynomial regression:

$$\min_{b_0, b_1, \ldots b_k}  \sum_{i=1}^N D_i \cdot K\left(\frac{T_i}{h}\right)\left(\hat{F}(T_i) - \sum_{j=0}^k \frac{b_j}{j!}T_i^j\right)^2 $$

$$\min_{b_0, b_1, \ldots b_k}  \sum_{i=1}^N (1-D_i) \cdot K\left(\frac{T_i}{h}\right)\left(\hat{F}(T_i) - \sum_{j=0}^k \frac{b_j}{j!}T_i^j\right)^2 $$
where $\hat{F}(t)$ is the empirical cdf $\hat{F}(t) = \frac{1}{n}\sum_{i=1}^n 1\{T_i \leq t\}$. Here, the only choice of tuning parameter is the bandwidth $h$. Cattaneo et al. (2020) provide rules for selecting $h$ in this context.

Package __rddensity__ implements the manipulation test based on Cattaneo et al.'s estimator. Let's use it in the senate data.


```r
library(rddensity)
```

```
## Warning: package 'rddensity' was built under R version 3.6.2
```

```r
rdd = rddensity(rdrobust_RDsenate$margin)
summary(rdd)
```

```
## 
## Manipulation testing using local polynomial density estimation.
## 
## Number of obs =       1390
## Model =               unrestricted
## Kernel =              triangular
## BW method =           estimated
## VCE method =          jackknife
## 
## c = 0                 Left of c           Right of c          
## Number of obs         640                 750                 
## Eff. Number of obs    408                 460                 
## Order est. (p)        2                   2                   
## Order bias (q)        3                   3                   
## BW est. (h)           19.841              27.119              
## 
## Method                T                   P > |T|             
## Robust                -0.8753             0.3814
```

```
## Warning in summary.CJMrddensity(rdd): There are repeated observations. Point
## estimates and standard errors have been adjusted. Use option massPoints=FALSE to
## suppress this feature.
```

```
## 
## P-values of binomial tests (H0: p=0.5).
## 
## Window Length / 2          <c     >=c    P>|T|
## 0.430                       8      12    0.5034
## 0.861                      17      25    0.2800
## 1.291                      25      34    0.2976
## 1.722                      45      47    0.9170
## 2.152                      51      55    0.7709
## 2.583                      66      65    1.0000
## 3.013                      79      71    0.5678
## 3.444                      94      86    0.6020
## 3.874                     105      94    0.4785
## 4.305                     115     107    0.6386
```

```r
plot <- rdplotdensity(rdd, rdrobust_RDsenate$margin)
```

![](rdd_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

The plot above computes a histogram for the running variable. The lines in this case come from the __local polynomial__ density regression of Cattaneo et al. (not a global polynomial!; function defaults to local quadratic regression); and the shaded areas correspond to robust (i.e. accounting for variance due to bias estimation) bias-corrected CIs.

## References

Calonico, S., Cattaneo, M.D. and Titiunik, R. (2014), Robust Nonparametric Confidence Intervals for Regression‐Discontinuity Designs. Econometrica, 82: 2295-2326. doi:10.3982/ECTA11757

Cattaneo, Matias D., Jansson, Michael and Ma, Xinwei (2020) Simple Local Polynomial Density Estimators, Journal of the American Statistical Association, 115:531, 1449-1455, DOI: 10.1080/01621459.2019.1635480

Gelman, Andrew and Imbens, Guido (2019). Why High-Order Polynomials Should Not Be Used in Regression Discontinuity Designs, Journal of Business & Economic Statistics, 37:3, 447-456, DOI: 10.1080/07350015.2017.1366909

Imbens, Guido and Kalyanaraman, Karthik (2012). Optimal Bandwidth Choice for the Regression Discontinuity Estimator, The Review of Economic Studies, Volume 79, Issue 3, July 2012, Pages 933–959, https://doi.org/10.1093/restud/rdr043

McCrary, J. (2008). Manipulation of the running variable in the regression discontinuity design: A density test. Journal of econometrics, 142(2), 698-714.
