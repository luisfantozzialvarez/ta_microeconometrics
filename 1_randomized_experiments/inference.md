---
title: "Some topics on inference"
author: "Luis Alvarez"
date: "8/27/2020"
output: 
  html_document:
    keep_md: true
---



## "Back-of envelope" power calculations

Consider a linear model:

$$Y_i =  \beta X_i + \gamma' Z_i + \epsilon_i$$

and suppose we are interested in parameter $\beta$. Most inferential methods you've seen rely on an estimator $\hat{\beta}$ of $\beta$ such that (at least asymptotically)  $(\hat{\beta}-\beta)/\sqrt{V(\hat{\beta})} \sim N(0,1)$, for all possible true values of $\beta$. Let us assume $V(\hat{\beta})$ is known. In this case, a (asymptotically) size $\alpha$ test of  $H_0: \beta = 0$ against the alternative $H_1: \beta \neq 0$ rejects the null if, and only if:

$$\hat{\beta}/\sqrt{V(\hat{\beta})} > q_{1-\alpha/2}  \text{ or } \hat{\beta}/\sqrt{V(\hat{\beta})} < q_{\alpha/2}  $$

where $q_{s}$ is the $s$ quantile of a standard normal distribution.

What is the (asymptotic) power of such a test against an alternative $\beta = b$, $b \in \mathbb{R}$? Let $\mathbb{P}_{b}$ denote the probability law when the true parameter is $b$. In this case, the probability of rejection is: 

$$\mathbb{P}_{b}[\hat{\beta}/\sqrt{V(\hat{\beta})}> q_{1-\alpha/2}  \text{ or } \hat{\beta}/\sqrt{V(\hat{\beta})} < q_{\alpha/2}] = \\ = 
\mathbb{P}_{b}[(\hat{\beta} - b)/\sqrt{V(\hat{\beta})} > q_{1-\alpha/2} - b/\sqrt{V(\hat{\beta})} ]  + P[(\hat{\beta}-b)/\sqrt{V(\hat{\beta})}< q_{\alpha/2} -b/\sqrt{V(\hat{\beta})}] = \\ = 1 - \Phi[q_{1-\alpha/2} - b/\sqrt{V(\hat{\beta})}] + \Phi[q_{\alpha/2} - b/\sqrt{V(\hat{\beta})}]$$

(here $\Phi$ denotes the standard normal cdf). So we have a quick formula to compute the power under a normal approximation. What if we don't know the variance? We could replace it by a consistent estimator of $V(\hat{\beta})$ (that is how we perform tests in practice anyway). So, if we have data (not even baseline) on a regression, we could use this quick formula to compute power under different alternatives. Importantly, though, there are two crucial things we should be cautious here. First, we are using data for which a particular value of $\beta=b$ is true, so we should be careful so that: (1) $V(\hat{\beta})$ is invariant to the true value of ${\beta}$; and $(2)$ the estimator we are using is consistent for $V(\hat{\beta})$, no matter what the true value of $\beta$ is. Item $(1)$ was used in the derivation above: it basically restricts what kind of alternative we are analysing. For example, in the homoskedastic world, $\mathbb{V}[\hat{\beta}] = \mathbb{E}\left[\begin{pmatrix}X \\ Z\end{pmatrix}\begin{pmatrix} X & Z \end{pmatrix}\right]^{-1} V(\epsilon)$, so we do not allow  $V(\epsilon)$ or the second moment matrix to vary with $\beta$ under alternatives (more on this later). Item $(2)$ works for most estimators we have in practice, though not all: estimators that impose the null when estimating the variance matrix won't fit here, as they are only consistent to the variance if the null is true. We also note that the above formula presumes a Gaussian approximation works, which may not be true. It also precludes more complex forms of inference, such as Wild BS. Notice these problems do not appear in the simulation-based approach to power calculation we discussed in the previous session (we can even compute the ``true'' size from the simulations!). 

Let us implement the above power calculation in the context of *followup* data on the experiment we discussed in the previous session. Before that, we note that, under the potential outcomes framework and random assignment of treatment, i.e:

\begin{equation}
Y_i = D_i Y_i(1) + (1-D_i) Y_i(1) \\
D_i \perp Y_i(0), Y_i(1)
\end{equation}
(where potential outcomes are stochastic because we assume we have a random sample from the population); we have that our model is:

\begin{equation}
Y_i = \alpha + \beta D_i + \epsilon_i \\
\alpha := \mathbb{E}[Y_i(0)] \\
\beta := \mathbb{E}[Y_i(1) - Y_i(0)] \\
\epsilon_i := D_i (Y_i(1)-Y_i(0)- \beta) + (Y_i(0)-\alpha)
\end{equation}

What does this imply? It implies that the power calculation using the formula above will vary the average treatment effect, while *keeping fixed treatment heterogeneity* at what is observed in the data. 


```r
head(data)
```

```
##     escola aluno     logico     verbal   abstrato    espacial mulher    idade
## 89       1    89 -0.2336050  1.1886198 -0.1476038 -0.56702017      1 10.30411
## 151      1   151 -1.1348661  0.5856582 -1.1028113 -1.61972552      0 10.64384
## 152      1   152 -1.1348661 -1.2232269 -0.1476038 -0.56702017      1 11.06027
## 235      1   235 -0.5340254  0.5856582  0.8076036  1.01203810      0 11.26301
## 281      1   281  0.3672357  0.5856582 -0.1476038 -0.04066742      1 11.04658
## 290      1   290 -0.2336050 -0.6202652  0.8076036  1.01203810      1 10.62466
##     tratamento
## 89           0
## 151          0
## 152          0
## 235          0
## 281          0
## 290          0
```

```r
summary(data)
```

```
##      escola          aluno           logico             verbal       
##  Min.   : 1.00   Min.   :  1.0   Min.   :-1.73571   Min.   :-2.4292  
##  1st Qu.: 8.00   1st Qu.:227.2   1st Qu.:-0.83445   1st Qu.:-0.6203  
##  Median :14.00   Median :453.5   Median : 0.06682   Median :-0.0173  
##  Mean   :15.19   Mean   :453.5   Mean   : 0.00000   Mean   : 0.0000  
##  3rd Qu.:23.00   3rd Qu.:679.8   3rd Qu.: 0.66766   3rd Qu.: 0.5857  
##  Max.   :30.00   Max.   :906.0   Max.   : 2.47018   Max.   : 2.9975  
##     abstrato          espacial            mulher           idade       
##  Min.   :-2.0580   Min.   :-1.61973   Min.   :0.0000   Min.   : 8.784  
##  1st Qu.:-0.6252   1st Qu.:-0.56702   1st Qu.:0.0000   1st Qu.:10.661  
##  Median :-0.1476   Median :-0.04067   Median :0.0000   Median :11.060  
##  Mean   : 0.0000   Mean   : 0.00000   Mean   :0.4735   Mean   :11.219  
##  3rd Qu.: 0.8076   3rd Qu.: 0.48569   3rd Qu.:1.0000   3rd Qu.:11.402  
##  Max.   : 2.2404   Max.   : 3.11745   Max.   :1.0000   Max.   :17.005  
##    tratamento    
##  Min.   :0.0000  
##  1st Qu.:0.0000  
##  Median :1.0000  
##  Mean   :0.5342  
##  3rd Qu.:1.0000  
##  Max.   :1.0000
```

```r
#Package to compute CR SE
library(sandwich)

#Package to produce test tables with robust SEs 
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
#Running the regression of interest
modelo = lm(logico~tratamento+idade+mulher, data)

#Inference with CR standard errors
coeftest(modelo, vcov. = vcovCL, cluster = data$escola)
```

```
## 
## t test of coefficients:
## 
##              Estimate Std. Error t value  Pr(>|t|)    
## (Intercept)  2.992809   0.335465  8.9214 < 2.2e-16 ***
## tratamento  -0.113388   0.140396 -0.8076    0.4195    
## idade       -0.273428   0.029664 -9.2175 < 2.2e-16 ***
## mulher       0.285773   0.049216  5.8065 8.836e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
#We'll use the formula to compute power
varCL = vcovCL(modelo, cluster = data$escola)

se = (sqrt(diag(varCL)))["tratamento"]

#Lets vary effects over a  grid
grid.beta = seq(-0.5,0.5, by = 0.05)

power = c()

alpha = 0.05
for(b in grid.beta)
  power = c(power,  1- pnorm(qnorm(1-alpha/2) - b/se ) + pnorm(qnorm(alpha/2) - b/se))

table.power = cbind("Avg. Effect" = grid.beta, "Power" = power)

print(table.power)
```

```
##            Avg. Effect      Power
## tratamento       -0.50 0.94535523
## tratamento       -0.45 0.89348177
## tratamento       -0.40 0.81303220
## tratamento       -0.35 0.70298292
## tratamento       -0.30 0.57020826
## tratamento       -0.25 0.42894883
## tratamento       -0.20 0.29653605
## tratamento       -0.15 0.18754474
## tratamento       -0.10 0.10983927
## tratamento       -0.05 0.06465238
## tratamento        0.00 0.05000000
## tratamento        0.05 0.06465238
## tratamento        0.10 0.10983927
## tratamento        0.15 0.18754474
## tratamento        0.20 0.29653605
## tratamento        0.25 0.42894883
## tratamento        0.30 0.57020826
## tratamento        0.35 0.70298292
## tratamento        0.40 0.81303220
## tratamento        0.45 0.89348177
## tratamento        0.50 0.94535523
```

How do these results compare with the simulation analysis we performed in the previous session?

What if we doubled the sample size? In this case, we know that standard errors are divided by $\sqrt{2}$, in which case we have:


```r
se = se/sqrt(2)


#Lets vary effects over a  grid
grid.beta = seq(-0.5,0.5, by = 0.05)

power = c()

alpha = 0.05
for(b in grid.beta)
  power = c(power,  1- pnorm(qnorm(1-alpha/2) - b/se ) + pnorm(qnorm(alpha/2) - b/se))

table.power = cbind("Avg. Effect" = grid.beta, "Power" = power)

print(table.power)
```

```
##            Avg. Effect      Power
## tratamento       -0.50 0.99895297
## tratamento       -0.45 0.99495756
## tratamento       -0.40 0.98073883
## tratamento       -0.35 0.94127899
## tratamento       -0.30 0.85587086
## tratamento       -0.25 0.71168287
## tratamento       -0.20 0.52182444
## tratamento       -0.15 0.32697244
## tratamento       -0.10 0.17188355
## tratamento       -0.05 0.07953038
## tratamento        0.00 0.05000000
## tratamento        0.05 0.07953038
## tratamento        0.10 0.17188355
## tratamento        0.15 0.32697244
## tratamento        0.20 0.52182444
## tratamento        0.25 0.71168287
## tratamento        0.30 0.85587086
## tratamento        0.35 0.94127899
## tratamento        0.40 0.98073883
## tratamento        0.45 0.99495756
## tratamento        0.50 0.99895297
```

Finally, we provide a function that computes the MDE for a given probability of rejecting the null $\kappa$. We do so by numerically inverting the formula provided above.


```r
#Function that computes MDE for bilateral test.
#alpha: significance level
#se: estimate of coefficient standard error
#kappa: rejection probability
#step.grid: grid step for for computing power
compute.mde.bilateral <- function(alpha, se, kappa, step.grid = 1e-3)
{
  #Grid for evaluating formula
  grid = seq(0, qnorm(1-alpha/2) + 5*se, by = step.grid*se)
  
  power = sapply(grid, function(b){1- pnorm(qnorm(1-alpha/2) - b/se ) + pnorm(qnorm(alpha/2) - b/se)})
  
  mde = grid[which.min(abs(power-kappa))]
  
  print(paste("MDE at ",100*kappa, "% of a ", alpha*100,"% test is absolute value of effect >= ",mde,sep=""))
  return(mde)
}
```


```r
#We'll use the formula to compute power
varCL = vcovCL(modelo, cluster = data$escola)

se = (sqrt(diag(varCL)))["tratamento"]

mde <- compute.mde.bilateral(0.05,se, 0.8)
```

```
## [1] "MDE at 80% of a 5% test is absolute value of effect >= 0.393389216231784"
```



## Wild BS

You have seen in class that Wild BS is a good alternative to Cluster-Robust standard errors when there are not many clusters in the sample. Indeed, it has been recently shown (Canay et al., 2020) that this method works well when: (1) there are many clusters; or (2) there are few clusters, but each cluster has *many* observations. The method also requires a specific ``homogeneity'' assumption on the distribution of covariates across clusters, which is made precise by Canay and coauthors. Recall that, for a model:

$$Y_{ic} = \beta X_{ic} + \gamma' Z_{ic} + \epsilon_{ic} $$
where $c \in \{1,2\ldots C\}$ indexes clusters; a Wild BS level $\alpha$ test for  $H_0: \beta = b$  against $H_1: \beta \neq b$ consists of:

1. Estimate the model imposing the null $\beta = b$. Let $\tilde{\gamma}$ and $\tilde{\epsilon_{ic}}$ denote the estimated coefficients and residuals from this step.
2. For replications $s \in \{1,2 \ldots S\}$:
    + For each $c \in \{1,2\ldots C\}$, draw $e_c^s = \pm 1$ with probability $1/2$.
    + Generate an artificial outcome using the formula $\tilde{y}^s_{ic} = bX_{ic} + \tilde{\gamma}'Z_{ic} + e_c^s \tilde{\epsilon}_{is}$ 
    + Run the unrestricted regression of $\tilde{y}^s_{ic}$ on $X_{ic}$ and $Z_{ic}$. Let $\hat{\beta}_s$ denote the OLS estimator obtained in this step.
    + Store the $s$-th coefficient $\hat{\beta}_s$ or its studentized version, $\hat{t}_s = (\hat{\beta}_s-b)/\sqrt{\hat{V}(\hat{\beta}_s)}$, where $\hat{V}(\hat{\beta}_s)$ is an estimate of the variance of the OLS coefficient (say, CR standard error).
3. Compute the $1-\alpha$ quantile of $\{|\hat{\beta}_s-b|:s \in 1,2\ldots S\}$ ($\{|\hat{t}_s|:s \in 1,2\ldots S\}$).
4. Reject the null if the absolute value of the unrestricted regression coefficient minus the value under the null $\hat{\beta}-b$ (unrestricted $t$-stat) *in the data* strictly exceeds the $(1-\alpha)$ quantile.

If we want to compute a p-value, we just need to find the smallest $\alpha$ for which the null is rejected.

Note that there is an option between studentizing or not the bootstrapped statistic. Canay et al. recommends studentization using cluster-robust standard errors, as it will work better when there are many clusters.

Let's create a function that implements Wild BS:


```r
#Function that computes wild BS. Takes as arguments:
#formula: amodel to test
#coef.to.test: character. name of the variable in formula whose coefficient we will test
#cluster.var: character. name of the cluster indicator variable
#data: dataframe where estimation will be conducted
#b: value of coefficient under the null. Defaults to 0
#S: Number of replications of step 2 in algorithm. Defaults to 1000
#dataset with variables
wild.bs <- function(formula, coef.to.test, cluster.var, data, b = 0, S = 1000)
{
  formula.text = as.character(formula)
  
  #Imposing the null in formula
  formula.null = paste(formula.text[2], "~", gsub(coef.to.test, paste("offset(",b,"*",coef.to.test,")",sep=""),formula.text[3]))
  
  modelo.nulo = lm(formula.null, data = data)
  
  cluster.data = cbind("Cluster"=data[,cluster.var])

  cluster.indexes = unique(cluster.data[,1])
  
  C = length(cluster.indexes)
  
  vec_unstud = c()
  vec_stud = c()
  
  data.artificial  = data
  
  for(s in 1:S)
  {
    e_s = 1 - 2*rbinom(C, 1, 0.5)
    
    vals.cluster = cbind("Cluster"=cluster.indexes, "e_s" = e_s)
    cluster.matched = merge(cluster.data, vals.cluster, by = "Cluster")
    
    #Creating artificial data
    
    data.artificial[,formula.text[2]] = modelo.nulo$fitted.values + cluster.matched$e_s*modelo.nulo$residuals
    
    modelo.s = lm(formula, data = data.artificial) 
    
    coef.s = modelo.s$coefficients[coef.to.test]
    
    vec_unstud = c(vec_unstud, coef.s)
    
    se.s = sqrt(diag(vcovCL(modelo.s, cluster = cluster.data[,1])))[coef.to.test]
    
    vec_stud = c(vec_stud,  (coef.s-b)/se.s)
  }
  
  #Compute estimates from the data now
  modelo.data = lm(formula, data = data)
  
  coef.data = modelo.data$coefficients[coef.to.test]
  
  p.val.unstud = 1 - mean(abs(coef.data-b) > abs(vec_unstud-b))
  se.data =  sqrt(diag(vcovCL(modelo.data, cluster = cluster.data[,1])))[coef.to.test]
  
  p.val.stud = 1 -mean(abs((coef.data-b)/se.data) > abs(vec_stud))
  
  return(list("Unstudentized p-value" = p.val.unstud, "Studentized p-value" = p.val.stud))
}
```


```r
#Setting seed to allow for replication
set.seed(1234)

#Sets seed to allow for replication
formula = logico ~ tratamento + idade + mulher

res.wild.bs = wild.bs(formula, "tratamento", "escola",data)

print(res.wild.bs)
```

```
## $`Unstudentized p-value`
## [1] 0.448
## 
## $`Studentized p-value`
## [1] 0.452
```

## Randomization inference

You have seen in class that randomization inference provides an *exact* method for conducting inference about _in-sample_ individual treatment effects in the context of randomized control trials -- where the treatment assignment mechanism is known. Indeed, we note that, under a *sharp* null of the type:

\begin{equation}
H_0; Y_j(0) + \delta_j = Y_j(1) \quad j = 1,2\ldots N
\end{equation}

where $\delta_j$ are _known constants_, potential outcomes in the sample are fully known. Let $\mathbf{Y}(d) = (Y_1(d)\ldots Y_N(d))'$, $d\in \{0,1\}$;  $\mathbf{D} = (D_1,\ldots D_N)'$; and $\mathbf{Y} = \mathbf{Y}(1)*\mathbf{D} + \mathbf{Y}(0)*(1-\mathbf{D})$ ($*$ denotes pointwise, entry-by-entry, product). We then have that, under the sharp null, for any statistic $\tau(\mathbf{Y},\mathbf{D})$, the distribution of $\tau(\mathbf{Y},\mathbf{D})$ is completely *known*, as the only source of randomness stems from the treatment assignment mechanism; and all potential outcomes are revealed. Indeed, since the treatment assignment mechanism is known, the distribution of $\tau(\mathbf{Y},\mathbf{D})$ under the null can be simulated by:

1. For $s \in \{1,2\ldots S\}$:
    * Draw $\mathbf{D}^s$ from the treatment assignment mechanism.
    * Generate artificial data $\mathbf{Y}^s = \mathbf{Y}(1)^n*\mathbf{D}^s + \mathbf{Y}(0)^n*(1-\mathbf{D}^s)$, where $\mathbf{Y}(d)^n$ denotes potential outcomes under the null.
    * Compute and store $\tau^s = \tau(\mathbf{Y}^s, \mathbf{D}^s)$.

An $\alpha$-level test can then be conducted by computing the test statistic in the observed data, and rejecting the null if this value exceeds the $1-\alpha$ quantile of the simulated distribution $\{\tau^s: s\in 1,2\ldots S\}$. By construction, this test has level $\alpha$.

In the discussion above, any function $\tau$ produces a test with correct level. Which statistic should we use then? Clearly, we also care about _power_, so we should use statistics for which high values indicate evidence against the null. For example, the absolute value of the difference in means between treated and control units; the studentized version of this statistic; coefficients of a treatment effect regression etc. Notice that, if we reject the null, then there is evidence of a violation $Y_j(0) + \delta_j = Y_j(1)$ for at least one individual in the sample. That may not always be useful.

Below, we implement randomization inference to test the null that there was no treatment effect on any individual in the sample: $\delta_j = 0$ for all $j$. We will use a studentized regression coefficient as statistic.


```r
#Running the regression of interest
modelo = lm(logico~tratamento+idade+mulher, data)

coef.data = modelo$coefficients["tratamento"]

#We'll use the formula to compute power
varCL = vcovCL(modelo, cluster = data$escola)

se.data = (sqrt(diag(varCL)))["tratamento"]

test.data = coef.data/se.data

escolas.indexes = unique(data$escola)
  
C = length(escolas.indexes)
  
data.artificial = data

vec_simulation = c()
for(s in 1:1000)
{
  
  #Generates draw from treatment assignment 
  tratamento.artificial  = sample(escolas.indexes, C/2, replace = F)
  data.artificial$tratamento = 1*(data.artificial$escola %in% tratamento.artificial)
  
  modelo.artificial = lm(logico~tratamento+idade+mulher, data.artificial) 
  
  coef.artificial = modelo.artificial$coefficients["tratamento"]
  
  #We'll use the formula to compute power
  varCL = vcovCL(modelo.artificial, cluster = data.artificial$escola)

  se.artificial = (sqrt(diag(varCL)))["tratamento"]

  vec_simulation = c(vec_simulation, coef.artificial/se.artificial)
}

p_value = 1 - mean(abs(test.data)>abs(vec_simulation))

print(p_value)
```

```
## [1] 0.437
```


## Multiple Hypothesis Testing in Randomization Inference

You have seen in class that, in the Randomization Inference framework, MHT can be performed quite easily. Indeed, consider testing the individual effect of a treatment on $K$ outcomes:

$$H_0: Y_j^k(0) + \delta_j^k  = Y_j^k(1) \quad j = 1\ldots N, k = 1 \ldots K$$
Then, as in the framework of the previous section, any statistic $\tau(\mathbf{Y}^1, \mathbf{Y}^2,\ldots \mathbf{Y}^K, \mathbf{D})$ has a known distribution under the null. Let us focus now in the case $\delta_j^k=0$ for all $j$,$k$. A statistic that could be used in testing is a Wald-type statistic

\begin{equation}
\begin{pmatrix} \hat{\beta}_1 & \hat{\beta}_2 & \ldots & \hat{\beta}_K \end{pmatrix} \times W \times \begin{pmatrix} \hat{\beta}_1 \\ \hat{\beta}_2 \\ \vdots \\ \hat{\beta}_K \end{pmatrix}
\end{equation}

where $\hat{\beta}_k$, $k = 1 \ldots K$ are estimates of average treatment effects for outcome $k$; and $W$ is a $K \times K$ weight matrix. For cases where $\delta_j^k \neq 0$, a similar statistic could be constructed by subtracting from each $\hat{\beta}^k$ the quantity $N^{-1}\sum_{j=1}^N \delta_j^k$, i.e. the average treatment effect in the sample under the null. Typical choices of weights $W$ would be:

\begin{equation}
W = \begin{pmatrix}
\hat{V}(\hat{\beta}_1)^{-1} & 0 & \ldots & 0\\
 0 & \hat{V}(\hat{\beta}_2)^{-1} & \ldots & 0\\
0 &  0 & \vdots & 0 \\
0 & 0 & \ldots & \hat{V}(\hat{\beta}_K)^{-1}
\end{pmatrix}
\end{equation}

or, as seen in class

\begin{equation}
W = V\left[\begin{pmatrix} \hat{\beta}_1 \\ \hat{\beta}_2 \\ \vdots \\ \hat{\beta}_K \end{pmatrix}\right]^{-1}
\end{equation}

In the second case, we need to compute the covariance between estimators in order to invert the covariance matrix. The easiest way to do so is, if we are running $K$ regressions:

$$Y_i^k = \beta_k D_i + \gamma_k'Z_{ik} + \epsilon_{ik} $$

is to create an outcome $S_i^k$ that stacks all outcomes and to run a saturated regression:

$$S_{i}^k = \sum_{p = 1}^K \beta_p \mathbb{1}\{p=k\} D_i + \sum_{p=1}^K\mathbb{1}\{p=k\}\gamma_p'Z_{ip} +  {\xi}_{ik}$$

With this regression, we can extract the covariance matrix from the R functions. We should also be cautious to cluster at the appropriate level and _always_ use standard errors that account for heteroskedasticity -- notice that the error in the new regression is heteroskedastic by construction (as it is the sum of errors from different models). Moreover, it is always good to check if you are doing things correctly by comparing standard errors from separate regressions with standard errors from the stacked regression (they should be the same, up to a degrees of freedom correction; we only stack so we can get covariances). We do this below.

Let's test the sharp null that there was no effect in any outcome for any individaul in our dataset.


```r
#We will stack our dataset so we can recover covariances of estimators from one regression
outcome.names = c("logico","verbal", "espacial", "abstrato")

#Lines that will be replicated along the stacking process
data.base = data[,c("escola", "aluno","idade","mulher","tratamento")]

data.stacked = data.frame()
for(outcome in outcome.names)
{
  data.to.stack = cbind(data.base, "outcome" = data[,outcome], "outcome.name" = outcome)
  data.stacked = rbind(data.stacked, data.to.stack )
}

head(data.stacked)
```

```
##     escola aluno    idade mulher tratamento    outcome outcome.name
## 89       1    89 10.30411      1          0 -0.2336050       logico
## 151      1   151 10.64384      0          0 -1.1348661       logico
## 152      1   152 11.06027      1          0 -1.1348661       logico
## 235      1   235 11.26301      0          0 -0.5340254       logico
## 281      1   281 11.04658      1          0  0.3672357       logico
## 290      1   290 10.62466      1          0 -0.2336050       logico
```

```r
#Converts outcome.label to factor (we assign an underlying enumeration to these variables)
data.stacked$outcome.name = as.factor(data.stacked$outcome.name)

#We now run the saturated-by-outcome regression 
modelo.full.data = lm(outcome~-1+outcome.name + outcome.name:tratamento
                      + outcome.name:idade + outcome.name:mulher, data = data.stacked)

#Separate hypothesis testing for avg treatment effect
#Notice inference is the same as running separately each regression, so we are actually doing things correctly
coeftest(modelo.full.data, vcov. = vcovCL, cluster = data.stacked$escola)
```

```
## 
## t test of coefficients:
## 
##                                  Estimate Std. Error t value  Pr(>|t|)    
## outcome.namelogico               2.992809   0.335603  8.9177 < 2.2e-16 ***
## outcome.nameverbal               1.714962   0.318678  5.3815 7.858e-08 ***
## outcome.nameespacial             1.305078   0.352844  3.6987 0.0002199 ***
## outcome.nameabstrato             1.623036   0.373957  4.3402 1.463e-05 ***
## outcome.namelogico:tratamento   -0.113388   0.140454 -0.8073 0.4195475    
## outcome.nameverbal:tratamento   -0.037957   0.097564 -0.3890 0.6972630    
## outcome.nameespacial:tratamento -0.046256   0.089702 -0.5157 0.6061187    
## outcome.nameabstrato:tratamento -0.108489   0.102395 -1.0595 0.2894394    
## outcome.namelogico:idade        -0.273428   0.029676 -9.2137 < 2.2e-16 ***
## outcome.nameverbal:idade        -0.148844   0.027739 -5.3659 8.560e-08 ***
## outcome.nameespacial:idade      -0.113487   0.029192 -3.8876 0.0001031 ***
## outcome.nameabstrato:idade      -0.140890   0.031646 -4.4521 8.761e-06 ***
## outcome.namelogico:mulher        0.285773   0.049237  5.8041 7.031e-09 ***
## outcome.nameverbal:mulher       -0.052430   0.069125 -0.7585 0.4482160    
## outcome.nameespacial:mulher     -0.015152   0.065794 -0.2303 0.8178807    
## outcome.nameabstrato:mulher      0.032835   0.071611  0.4585 0.6466067    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
#We now extract the relevant part of the covariance matrix
vcov.data = vcovCL(modelo.full.data,cluster = data.stacked$escola)

#Extract just subvector/submatrix with treatment coefficients
extract = grepl("tratamento",colnames(vcov.data))

coef.data = modelo.full.data$coefficients[extract]
vcov.data = vcov.data[extract, extract]

#Constructing weight
Weight.data = solve(vcov.data)

#Wald statistic
wald.data = as.numeric(t(coef.data)%*%Weight.data%*%coef.data)

#Randomization inference now
data.artificial = data.stacked
vec_simulation = c()
for(s in 1:1000)
{
  
  #Generates draw from treatment assignment 
  tratamento.artificial  = sample(escolas.indexes, C/2, replace = F)
  data.artificial$tratamento = 1*(data.artificial$escola %in% tratamento.artificial)
  
  modelo.artificial = lm(outcome~-1+outcome.name + outcome.name:tratamento
                      + outcome.name:idade + outcome.name:mulher, data.artificial) 
  
  coef.artificial = modelo.artificial$coefficients[extract]
  
  #We'll use the formula to compute power
  varCL = vcovCL(modelo.artificial, cluster = data.artificial$escola)

  vcov.artificial = varCL[extract, extract]
  
  Weight.artificial =  solve(vcov.artificial)

  vec_simulation = c(vec_simulation, t(coef.artificial)%*%Weight.artificial%*%coef.artificial)
}

#P-value from randomization inference
p_value = 1 - mean(wald.data > vec_simulation)

print(p_value)
```

```
## [1] 0.852
```


## References

Canay, I.A., A. Santos, and A.M. Shaikh (2020): “The Wild Bootstrap with a “Small” Number of “Large” Clusters”,
The Review of Economics and Statistics, forthcoming. 

