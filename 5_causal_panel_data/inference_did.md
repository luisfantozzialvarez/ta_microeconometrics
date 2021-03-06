---
title: "Inference in Differences-in-Differences with few treated units"
author: "Luis Alvarez"
date: "11/22/2020"
output: 
  html_document:
    keep_md: true
---



Consider a panel on $G$ units and $T$ periods. Let $D_{gt}$ denote a binary treatment which begins (and persists thereafter) at period $t^*\leq T$ for a subset $\mathcal{G}_1 \subseteq \{1,\ldots,G\}$ of units. Suppose the following _model_ for an outcome of interest holds:

$$y_{gt} = \alpha_g + \gamma_t + \beta D_{gt} + \epsilon_{gt} \\ \mathbb{E}[\epsilon_{gt}|\mathbf{D}] = 0 \\ g = 1,\ldots G \\
t = 1,\ldots T$$

where $\mathbf{D} := [(D_{11}.\ldots D_{1T})', \ldots, (D_{G1},\ldots D_{GT})']'$. We also define $\mathcal{G}_0:=\{1,\ldots,G\}  \setminus  \mathcal{G}_1$; $G_1:=|\mathcal{G}_1|$ and $G_0 := |\mathcal{G}_0|$ . You have seen in class that, under an asymptotic sequence such that $G_0, G_1 \to \infty$ and the proper technical conditions, the two-way fixed effects (TWFE) estimator of $\beta$ in this context -- which collapses to a pre-post-period-average DID estimator -- is consistent and asymptotically normal. In this context, cluster-robust standard errors provide inference that is robust to (possible) serial correlation within each group.

What happens when $G_1$ is small? In this context, Ferman and Pinto (2019) argue that the cluster-robust variance estimator (CRVE) will tend to underestimate the variance. Similarly, the wild bootstrap we have seen in the first course -- which under some conditions works with both $G_1$ and $G_0$ small __and $T$ large__ -- will tend not to work either, because the symmetry conditions required by it will hardly be satisfied in a setting where there are not,  __simultaneously__ within each cluster, observations from $\mathcal{G}_1$  and $\mathcal{G}_0$ (see Example 2.3. and the Appendix in Canay et al. (2018); in principle, they argue this could be circumvented by clustering at a higher level). One method that __does work__ with small $G_1$ and _large_ $G_0$ is the method of Conley and Taber (2011). They show that, under between-group weak dependence and moment assumptions, the TWFE estimator of the model above converges in probability, as $G_0 \to \infty$:

$$\hat{\beta} \overset{p}{\to} \beta  + W \\
W := \frac{\sum_{g \in \mathcal{G}_1} \sum_{t=1}^T (D_{gt} - \bar{D}_g) (\epsilon_{gt} - \bar{\epsilon}_g) }{\sum_{g \in \mathcal{G}_1} \sum_{t=1}^T(D_{gt} - \bar{D}_g)^2}$$
i.e., the estimator is not consistent, as it converges to a random variable. Conley and Taber (2011) show, nonetheless, that the asymptotic representation above enables valid inference. The idea is to approximate the distribution of $W$ by using the __residuals in the control group__ (why don't we exploit variation in residuals from the treatment group?). For that, they require the assumption:

__CT.Assumption 2__ $(\epsilon_{g1},\ldots \epsilon_{gT})$ is iid across $g$, and __independent__ of $(D_{g1}, D_{g2}\ldots D_{gT})$, with a bounded density.

Under the assumption above, it can be shown that the distribution of the residuals in the control group will approximate the distribution of the _errors_ in the treatment group. The inference procedure then goes as follows. To test the null that $\hat{\beta} = \beta_0$ at the $\alpha$ significance level.

1. Run the TWFE regression. Store the point estimate $\hat{\beta}$ and the residuals $\hat{\epsilon}_{jt}$.
2. For replications $b=1,\ldots B$:

    2.1. Sample (with replacement) $G_1$ units __from the control group__ $\mathcal{G}_0$. Let $\{i_{g,b}: g \in \mathcal{G}_1\}$ denote the sampled indices, which we "pair" to a treated unit.
    
    2.2. Construct and store the statistic $\hat{W}^b := \frac{\sum_{g \in \mathcal{G}_1} \sum_{t=1}^T (D_{gt} - \bar{D}_g) {\color{red}{\hat{\epsilon}}_{{i_{g,b}},t} }}{\sum_{g \in \mathcal{G}_1} \sum_{t=1}^T(D_{gt} - \bar{D}_g)^2}$.
    
3. If $\hat{\beta}-\beta_0 \in (-\infty, c_{\alpha/2}) \cup (c_{1-\alpha/2},\infty)$, where $c_q$ is the $q$-quantile of the empirical distribution of the $\hat{W}^b$, then one rejects the null.

We note that step 3 can be replaced by looking at the distribution of $|\hat{\beta}-\beta_0|$. Moreover, one may want to impose the null when estimating residuals, though in this case a permutation-style test can be constructed by also looking at the residuals in the treatment group (cf. Conley and Taber for details (2011); in their simulations, this approach performs better in terms of size when there are fewer groups, though at the cost of slight power losses). Finally, confidence intervals can be constructed by inverting the test procedure above.

The approach of Conley and Taber easily generalizes to controls entering linearly in the model (we just residualise things) and to nonbinary treatments. In the latter, control groups are those for which $D_{gt}$ does not vary over time. Variation in timing is also accounted for automatically by the results above.

Below, we present a function that implements Conley and Taber's apprach. Once again, we use the __fixest__ package to estimate TWFE models.


```r
library(fixest)

#Setting seed to allow for replication
set.seed(1234)

#Function that implements CT approach to a group-level model:
# y_{gt} = \alpha_g + \gamma_t + \beta'D_{gt} + \omega'X_{gt} \epsilon_{it} 
# Arguments are
# outcome = character vector for the variable
# treatment = character vector for the treatment variable
# timevar = character vector for the time variable. Should be orderable
# groupvar = character vector for the group variable
# data = data.frame
# beta_null = the value being tested under the null
# controls = vector of group-level controls to be included. If c() (default), no controls are included.
# Breps = number of bootstrap replications
conley_taber_test <- function(outcome, treatment, timevar, groupvar, data, beta_null = 0, controls = c(),
                             Breps = 10000){
  
  #Ordering data
  data = data[order(data[,groupvar], data[,timevar]),]
  
  #Model formula now
  if(length(controls)>0)
    formula.model = paste(outcome,"~",treatment,"+",paste(controls, collapse="+"), "|",paste("as.factor(",c(groupvar,timevar),")",sep="",collapse="+"), sep="") else formula.model = paste(outcome,"~",treatment, "|",paste("as.factor(",c(groupvar,timevar),")",sep="",collapse="+"), sep="")
  
  
  model = feols(as.formula(formula.model), data = data)
  
  
  data$group= data[,groupvar]
  data$time = data[,timevar]
  data$treatment = data[,treatment]
  
  data$residuals = resid(model)
  
  #Demeaning treatment
  m.resid_treat = feols(as.formula("treatment~1|as.factor(group)"),data)
  data$treat_demean = resid(m.resid_treat)
  
  
  glist = unique(data$group)
  
  #Computing treatment and control group
  test.G1 = sapply(glist, function(x){
    sum(diff(data[data$group==x,"treatment"])!=0)>0
  })
  
  G1= glist[test.G1]
  G0 = glist[!test.G1]
  
  NG1 = length(G1)
  NG0 = length(G0)
  
  is_G1= data$group%in%G1
  
  
  data.G1 = data[is_G1,]
  
  data.G0 = data[!is_G1,]
  
  W_vec = c()
  
  for(b in 1:Breps){
    #Sampling from G0 with replacement
    draw = sample(G0, NG1, replace=T)
    
    data_artificial = cbind("group_G1" = G1,"group_G0" = draw)
    
    data_artificial =  merge(data_artificial, data.G1[,c("group","time","treat_demean")], by.x = "group_G1", by.y = "group", all.x = T, all.y=F)
    
    data_artificial =  merge(data_artificial, data.G0[,c("group","time","residuals")], by.x = c("group_G0","time"), by.y = c("group","time"), all.x = T, all.y=F)
    
    
    W_b = sum(data_artificial$treat_demean*data_artificial$residuals,na.rm = T)/sum(data_artificial$treat_demean^2, na.rm=T)
    
    W_vec = c(W_vec, W_b)    
  }
  
  #We are ready to compute pvalues
  prob.equitailed = mean((model$coefficients[treatment]-beta_null)<=W_vec)
  
  if(prob.equitailed>=0.5)
    p.value.equitailed =  mean((model$coefficients[treatment]-beta_null)>=W_vec)*2 else
      p.value.equitailed =  mean((model$coefficients[treatment]-beta_null)<=W_vec)*2 

  p.value.absvalue =mean(abs(model$coefficients[treatment]-beta_null)<=abs(W_vec))
  
  
  return(list("G1"=G1, "G0" = G0, "P-value (equitailed)"=p.value.equitailed, "P-value (absolute value)"=p.value.absvalue))
  
  
}
```

In what follows, we apply our function to the data on Conley and Taber (2011), whose empirical application studies the impact of merit-based scholarship programs at the state level on college attendance. We aggregate their data at the state-level -- the authors provide an extension of their method which allows for individual-level covariates which we do not pursue here. Their data can be found at: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/GCBK24. I'm using the "merit" dataset,



```r
dados = read.csv("base_state_final.csv")

merit_data = read.csv("merit.csv",sep = " ")

#Treating data 
merit_data = merit_data[,!is.na(merit_data[1,])]
colnames(merit_data) = c("coll",  "merit", "male", "black", "asian", "year", "state", "chst")

merit_data$nobs = 1

merit_state =  aggregate(merit_data[,!(colnames(merit_data)%in%c("year","state"))], by = list("state"= merit_data$state,
                                                                                              "year" = merit_data$year),
                         FUN = sum)

for(var in colnames(merit_state))
    if(!(var%in%c("year","state","nobs")))
      merit_state[,var] = merit_state[,var]/merit_state$nobs


head(merit_state)
```

```
##   state year      coll merit      male      black      asian chst nobs
## 1    11 1989 0.4358974     0 0.4102564 0.00000000 0.00000000    0   39
## 2    12 1989 0.4615385     0 0.6153846 0.00000000 0.05128205    0   39
## 3    13 1989 0.5555556     0 0.4722222 0.00000000 0.00000000    0   36
## 4    14 1989 0.5229885     0 0.5229885 0.01724138 0.04022989    0  174
## 5    15 1989 0.4324324     0 0.5405405 0.02702703 0.00000000    0   37
## 6    16 1989 0.5000000     0 0.5750000 0.15000000 0.00000000    0   40
```

```r
#Let's take a look at treatment timing
#Ordering data
merit_state = merit_state[order(merit_state$state, merit_state$year),]

#Timing
print(sapply(unique(merit_state$state), function(x){(merit_state$merit[merit_state$state==x])}))
```

```
##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
##  [1,]    0    0    0    0    0    0    0    0    0     0     0     0     0
##  [2,]    0    0    0    0    0    0    0    0    0     0     0     0     0
##  [3,]    0    0    0    0    0    0    0    0    0     0     0     0     0
##  [4,]    0    0    0    0    0    0    0    0    0     0     0     0     0
##  [5,]    0    0    0    0    0    0    0    0    0     0     0     0     0
##  [6,]    0    0    0    0    0    0    0    0    0     0     0     0     0
##  [7,]    0    0    0    0    0    0    0    0    0     0     0     0     0
##  [8,]    0    0    0    0    0    0    0    0    0     0     0     0     0
##  [9,]    0    0    0    0    0    0    0    0    0     0     0     0     0
## [10,]    0    0    0    0    0    0    0    0    0     0     0     0     0
## [11,]    0    0    0    0    0    0    0    0    0     0     0     0     0
## [12,]    0    0    0    0    0    0    0    0    0     0     0     0     1
##       [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25]
##  [1,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [2,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [3,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [4,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [5,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [6,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [7,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [8,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [9,]     0     0     0     0     0     0     0     0     0     0     0     0
## [10,]     0     0     0     0     0     0     0     0     0     0     0     0
## [11,]     0     0     0     0     0     0     0     0     0     0     0     0
## [12,]     0     0     0     0     0     0     0     0     0     0     0     0
##       [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35] [,36] [,37]
##  [1,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [2,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [3,]     0     0     0     0     0     0     0     0     0     1     0     0
##  [4,]     0     0     0     0     0     0     0     0     0     1     0     0
##  [5,]     0     0     0     1     0     0     0     0     0     1     0     0
##  [6,]     0     0     0     1     0     0     0     0     0     1     0     0
##  [7,]     0     0     0     1     0     0     0     0     0     1     0     0
##  [8,]     0     0     0     1     0     0     0     0     1     1     0     0
##  [9,]     0     0     0     1     1     0     0     0     1     1     0     0
## [10,]     0     0     1     1     1     0     0     0     1     1     1     0
## [11,]     0     0     1     1     1     1     0     0     1     1     1     0
## [12,]     0     0     1     1     1     1     0     0     1     1     1     0
##       [,38] [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46] [,47] [,48] [,49]
##  [1,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [2,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [3,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [4,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [5,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [6,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [7,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [8,]     0     0     0     0     0     0     0     0     0     0     0     0
##  [9,]     0     0     0     0     0     1     0     0     0     0     0     0
## [10,]     0     0     0     0     0     1     0     0     0     0     0     0
## [11,]     0     0     0     0     0     1     0     0     0     0     0     0
## [12,]     0     0     0     0     0     1     0     0     1     0     0     0
##       [,50] [,51]
##  [1,]     0     0
##  [2,]     0     0
##  [3,]     0     0
##  [4,]     0     0
##  [5,]     0     0
##  [6,]     0     0
##  [7,]     0     0
##  [8,]     0     0
##  [9,]     0     0
## [10,]     0     0
## [11,]     0     0
## [12,]     0     0
```

```r
#TWFE regression
model = feols(coll~merit|as.factor(year)+as.factor(state),merit_state )
summary(model, se = "cluster", cluster = "state")
```

```
## OLS estimation, Dep. Var.: coll
## Observations: 612 
## Fixed-effects: as.factor(year): 12,  as.factor(state): 51
## Standard-errors: Clustered (state) 
##       Estimate Std. Error t value Pr(>|t|)    
## merit 0.054337   0.010698  5.0792 5.65e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## Log-likelihood: 703.34   Adj. R2: 0.32192 
##                        R2-Within: 0.01473
```

```r
#Conley-Taber
conley_taber_test("coll", "merit", "year", "state", merit_state)
```

```
## $G1
##  [1] 34 57 58 59 61 64 71 72 85 88
## 
## $G0
##  [1] 11 12 13 14 15 16 21 22 23 31 32 33 35 41 42 43 44 45 46 47 51 52 53 54 55
## [26] 56 62 63 73 74 81 82 83 84 86 87 91 92 93 94 95
## 
## $`P-value (equitailed)`
## [1] 0.0094
## 
## $`P-value (absolute value)`
## [1] 0.0116
```

Notice that Assumption 2 of Conley and Taber __prohibits__ heteroskedasticity. Ferman and Pinto (2019) argue this may be restrictive, especially if the group-level model comes from aggregating (taking means of) individual-level observations. For example, in our example above, our state-level panel comes from aggregating individual-level observations. Suppose the individual-level model is:

$$Y_{igt} = \alpha_g + \gamma_t + \beta D_{gt} + v_{gt} +  \xi_{igt} $$
then, once we aggregate the model, $\epsilon_{gt} = v_{gt}+ \frac{1}{N_{gt}} \sum_{i=1}^{N_{gt}} \xi_{igt}$. Ferman and Pinto show that, if, conditional on $\mathbf{D}$, one assumes $\xi_{igt}$ is iid across $i$, $g$ and $t$; $v_{gt}$ is iid across $g$ (though arbitrarily dependent across $t$), and $N_{gt}$ is constant in $t$, then:

$$\mathbb{V}\left[\frac{1}{T-t^*}\sum_{t=t^*}^T \epsilon_{gt} - \frac{1}{t^*-1} \sum_{t=1}^{t^*-1} \epsilon_{gt} |\mathbf{D}\right] = A + \frac{B}{N_{g}}, \quad \forall g \in 1, \ldots G \quad (*)$$

(notice that $\sum_{g \in \mathcal{G}_1} \frac{1}{T-t^*}\sum_{t=t^*}^T \epsilon_{gt} - \frac{1}{t^*} \sum_{t=1}^{t^*} \epsilon_{gt}, g \in \mathcal{G}_1$, is precisely the summand in the $W$ formula previously outlined). Therefore, variation in group sizes induces heteroskedasticity in the group-level model. If the number of observations in treated and control groups is very different, then the approach of Conley and Taber may lead to size distortions. Ferman and Pinto provide simulation evidence showing this is indeed a problem.

The authors thus propose a corrected bootstrap procedure to account for this problem. They suggest estimating the conditional heteroskedasticity model and then rescaling bootstrapped residuals by the conditional variance of the "paired" observation. The main required assumption is that the _same_ model for the conditional variance holds in both groups. Their bootstrap procedure starts from the observation that, in our setting, the TWFE estimator (which collapses to DID):

$$ \hat{\beta} - \beta = \frac{1}{G_1}\left[ \sum_{g \in \mathcal{G}_1} \frac{1}{T-t^*}\sum_{t=t^*}^T \epsilon_{gt} - \frac{1}{t^*-1} \sum_{t=1}^{t^*-1} \epsilon_{gt}\right] - \frac{1}{G_0}\left[ \sum_{g \in \mathcal{G}_0} \frac{1}{T-t^*}\sum_{t=t^*}^T \epsilon_{gt} - \frac{1}{t^*-1} \sum_{t=1}^{t^*-1} \epsilon_{gt}\right] =: \frac{1}{G_1} \sum_{g\in \mathcal{G}_1} W_g - \frac{1}{G_0} \sum_{g\in \mathcal{G}_0} W_g  $$
They then propose the following bootstrap procedure, which uses residuals estimated under the null $\beta = \beta_0$ (so we sample from both groups):

1. Calculate $\hat{\beta}$ via TWFE.
2. Estimate the model __with the NULL imposed__. Store $\tilde{W}_g = \frac{1}{T-t^*}\sum_{t=t^*}^T \tilde{\epsilon}_{gt} - \frac{1}{t^*-1} \sum_{t=1}^{t^*-1} \tilde{\epsilon}_{gt}$ .
3. Estimate the conditional variance model using the $\tilde{W}_g$ (e.g. the variation-in-group-sizes model (*) may be estimated using OLS of $\tilde{W}_g^2$ on an intercept and $1/N_g$). Let $\hat{V}[W_g|\mathbf{D}]$ denote the estimated variance from such a model.
4. For replications $b = 1,2\ldots B$:

    4.1. Sample $G_0+G_1$ units with replacement. Let $i_{g,b}, g \in \mathcal{G}_0 \cup \mathcal{G}_1$ denote the sampled indices, which we "pair" to a unit in the original sample.
  
    4.2. Compute and store $\hat{\beta}^b = \frac{1}{G_1}\sum_{g\in \mathcal{G}_1} \tilde{W}_{i_{g,b}} \sqrt{\frac{\hat{V}[W_{i_{g,b}}|\mathbf{D}]}{\hat{V}[W_g|\mathbf{D}]}} - \frac{1}{G_0}\sum_{g\in \mathcal{G}_1}  \tilde{W}_{i_{g,b}} \sqrt{\frac{\hat{V}[W_{i_{g,b}}|\mathbf{D}]}{\hat{V}[W_g|\mathbf{D}]}}$. 
  
5. Reject the null if $\hat{\beta}$ is at the tails of the distribution $\hat{\beta}^b$.

The authors propose an extension to the above approach which allows for variation in treatment timing by showing that in this case the estimator asymptotically corresponds to a weighted average of DID estimators -- each one with a different treatment timing --, and then, for each simulation draw, computing resampled estimates for each timing and aggregating these. If there are $K$ different treatment starting periods, they estimate $K$ conditional variance models. They also use a consistency result in Conley and Taber (2011) to show that the inclusion of covariates in the outcome model is straightforward: one just includes those in the regressions (steps 1 and 2), and does not have to adapt the $\hat{\beta}^b$ formula.

We implement the approach of Ferman and Pinto in the leading example of variation in group sizes. The function below uses the __average__ number of observations across time within a group in the model (*).


```r
#Function that implements FP approach to a group-level model with variation in group sizes
# y_{gt} = \alpha_g + \gamma_t + \beta'D_{gt} + \alpha'X_{gt} \epsilon_{it} 
# Var(W_g) = A + B/N_g
# Arguments are
# outcome = character vector for the variable
# treatment = character vector for the BINARY treatment variable
# timevar = character vector for the time variable. Should be orderable
# groupvar = character vector for the group variable
# nobsvar = variable indicating the number of observations in each group (group-averages will be taken before regression model is specified)
# data = data.frame
# beta_null = the value being tested under the null
# controls = vector of group-level controls to be included. If c() (default), no controls are included.
# Breps = number of bootstrap replications
ferman_pinto_test <- function(outcome, treatment, timevar, groupvar, nobsvar, data, beta_null = 0, controls = c(),
                             Breps = 10000){
  
  #Ordering data
  data = data[order(data[,groupvar], data[,timevar]),]
  
  #Model formula now
  if(length(controls)>0)
    formula.model = paste(outcome,"~",treatment,"+",paste(controls, collapse="+"), "|",paste("as.factor(",c(groupvar,timevar),")",sep="",collapse="+"), sep="") else formula.model = paste(outcome,"~",treatment, "|",paste("as.factor(",c(groupvar,timevar),")",sep="",collapse="+"), sep="")
 
    
  model = feols(as.formula(formula.model), data = data)


  
  #Model formula under the null
  data$outcome_null = data[,outcome] - beta_null*data[,treatment] 
  
  if(length(controls)>0)
    formula.model.null = paste("outcome_null ~",paste(controls, collapse="+"), "|",paste("as.factor(",c(groupvar,timevar),")",sep="",collapse="+"), sep="") else formula.model.null = paste("outcome_null~1", "|",paste("as.factor(",c(groupvar,timevar),")",sep="",collapse="+"), sep="")
  
  
  model.null = feols(as.formula(formula.model.null), data = data)

  
  data$group= data[,groupvar]
  data$time = data[,timevar]
  data$treatment = data[,treatment]
  data$sizevar = data[,nobsvar]
  
  data$residuals = resid(model.null)
  
  
  #Creating timing variable
  glist = unique(data$group)
  Tt = length(unique(data$time))
  tmax = max(data$time)

  
  mat_timing = data.frame("group" = glist, "tpost" =sapply(glist, function(x){sum(data$treatment[data$group==x],na.rm = T)}),
                          stringsAsFactors = F)
  
  
  mat_to_add = aggregate(data["sizevar"],by = list("group"= data$group), FUN = mean, na.rm = T)
  mat_to_add = mat_to_add[order(mat_timing$group),]
  
  mat_timing = merge(mat_timing, mat_to_add, by = "group")
  
  mat_timing = mat_timing[order(mat_timing$group),]
  
  tpost_list =  unique(mat_timing$tpost)
  
  residual_mat = c()
  sd_mat = c()
  
  for(s in tpost_list)
    if(s!=0){
      
      data_case = data
      data_case$residual_to_aggregate = (data_case$time > (tmax - s))*data_case$residuals/s - (data_case$time <= (tmax - s))*data_case$residuals/(Tt-s) 
      
      agg_mod_case = aggregate(data_case["residual_to_aggregate"], by =list("group" = data_case$group), FUN = sum,na.rm=T)
      
      agg_mod_case = merge(agg_mod_case, mat_timing, by = "group")
      
      agg_mod_case = agg_mod_case[order(agg_mod_case$group),]
      
      agg_mod_case$res_sqd = agg_mod_case$residual_to_aggregate^2
      
      agg_mod_case$inv_size = 1/agg_mod_case$sizevar
      
      cvar_mod = lm(res_sqd~1+inv_size, agg_mod_case)
      
      residual_mat = cbind(residual_mat, agg_mod_case$residual_to_aggregate)
      sd_mat = cbind(sd_mat, sqrt(predict(cvar_mod, agg_mod_case)))
      
    }
  
  beta_vec = c()
  
  
  weight_vec = sapply(tpost_list[tpost_list!=0], function(x){sum(mat_timing$tpost==x)*(Tt-x)*x})
  weight_vec = weight_vec/sum(weight_vec)
  
  for(b in 1:Breps){
    #Sampling from G0 with replacement
    draw = sample(1:length(glist), length(glist), replace=T)
    
    res_reshuffle = residual_mat[draw,]
    sd_correct = sd_mat[draw,]/sd_mat
    res_correct = res_reshuffle*sd_correct
    
    did_pair = sapply(1:sum(tpost_list!=0), function(x){
      tt = (tpost_list[tpost_list!=0])[x]
      
      return(mean(res_correct[mat_timing$tpost==tt,x]) - mean(res_correct[mat_timing$tpost==0,x])) 
    })
    
    
    beta_vec = c(beta_vec, sum(weight_vec*did_pair))
  }
    

  #We are ready to compute pvalues
  prob.equitailed = mean((model$coefficients[treatment])<=beta_vec)
  
  if(prob.equitailed>=0.5)
    p.value.equitailed =  mean((model$coefficients)>=beta_vec)*2 else
      p.value.equitailed =  mean((model$coefficients[treatment])<=beta_vec)*2 

  p.value.absvalue =mean(abs(model$coefficients[treatment])<=abs(beta_vec))
  
  
  return(list("Tpost_distr" = data.frame("Tpost" = tpost_list[tpost_list!=0],"Weight" = weight_vec), "P-value (equitailed)"=p.value.equitailed, "P-value (absolute value)"=p.value.absvalue))

}
```


```r
summary(merit_state$nobs)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   15.00   37.00   48.00   68.89   63.00  383.00
```

```r
ferman_pinto_test("coll", "merit", "year", "state", "nobs", merit_state)
```

```
## $Tpost_distr
##   Tpost     Weight
## 1     1 0.08906883
## 2     3 0.21862348
## 3     8 0.12955466
## 4     4 0.25910931
## 5     2 0.08097166
## 6     5 0.14170040
## 7    10 0.08097166
## 
## $`P-value (equitailed)`
## [1] 0.027
## 
## $`P-value (absolute value)`
## [1] 0.0308
```


## References


Canay, I. A., Santos, A., & Shaikh, A. M. (2018). The wild bootstrap with a "small" number of" large" clusters. Review of Economics and Statistics, 1-45.

Conley, T. G., & Taber, C. R. (2011). Inference with “difference in differences” with a small number of policy changes. The Review of Economics and Statistics, 93(1), 113-125.

Ferman, Bruno and Pinto, Cristine, (2019), <a href="https://EconPapers.repec.org/RePEc:tpr:restat:v:101:y:2019:i:3:p:452-467">Inference in Differences-in-Differences with Few Treated Groups and Heteroskedasticity</a>, <i>The Review of Economics and Statistics</i>, <b>101</b>, issue 3, p. 452-467.

