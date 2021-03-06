---
title: "Two-way fixed effects regressions with heterogenous effects"
author: "Luis Alvarez"
date: "11/16/2020"
output: 
  html_document:
    keep_md: true
---




Consider a setting where we have access to data on $G$ groups over $T$ periods. At each period $t$, there are $N_{g,t}\geq 1$ units in each group, for which we observe outcomes $Y_{i,g,t}$, $i=1,2 \ldots N_{g,t}$. We also have information on a binary __group__-level treatment, $D_{g,t}$.

Suppose one wants to estimate the effect of the treatment on the observed outcomes. A common approach is to estimate the two-way fixed-effect (TWFE) regression:


$$  Y_{i,g,t} = \alpha_{g} + \gamma_t + \beta D_{g,t} + \epsilon_{i,g,t}$$

In a setting where treatment adoption is staggered ($D_{g,t+1}\geq D_{g,t} \ \forall t,g$), and treatment begins __simultaneously__ for all _eventually_ treated groups, it is possible to show that the approach above collapses to a differences-in-differences design comparing pre- and post- treatment averages between treated and control groups. In this case, we know that the OLS estimator $\hat{\beta}$ is unbiased (consistent) for the ATT under a parallel trends assumption. 

However, the staggered simultaneous-start setting is not the only one used in practice. What does the TWFE regression estimate in more general settings? de Chaisemartin and D’Haultfeuille (2020) study this problem in detail. Define $Y_{i,g,t}(1), Y_{i,g,t}(0)$ potential outcomes for treatment ocurring __in $t$__ for unit $i$ in group $g$ at period $t$. Let $Y_{g,t}(d): = \frac{1}{N_{g,t}} \sum_{i=1}^{N_{g,t}} Y_{i,g,t}(d)$, $d \in \{0,1\}$, denote group-level potential outcomes. The authors consider a framework where $(Y_{g,t}(0), Y_{g,t}(1), D_{g,t})$ is random. They show that, under the assumptions:

__Assumption 1 (sampling)__  ${(Y_{g,t}(0), Y_{g,t}(1), D_{g,t})}_{t=1}^T$ is independent across g.


__Assumption 2 (exogeneity)__ For all $g,t$, $\mathbb{E}[Y_{g,t+1}(0) - Y_{g,t}(0)|D_{g,1}\ldots D_{g,T}] = \mathbb{E}[Y_{g,t+1}(0) - Y_{g,t}(0)]$.


__Assumption 3 (common trends) __ $\mathbb{E}[Y_{g,t+1}(0) - Y_{g,t}(0)]$ does not depend on $g$.

the TWFE estimator satisfies:

$$\mathbb{E}[\hat{\beta}_{\text{TWFE}}|\mathbf{D}] = \sum_{(g,t): D_{g,t} = 1} \frac{N_{g,t}}{N_1} \cdot w_{g,t} \cdot \mathbb{E}[Y_{i,t}(1) - Y_{i,t}(0)|\mathbf{D}] $$

where  $\mathbf{D} = ((D_{1,1}\ldots D_{1,T})', \ldots (D_{G,1}, D_{G,2}\ldots D_{G,T})')'$; and $N_1$ is the number of treated observations across all periods $N_1 = \sum_{g,t} N_{g,t} D_{g,t}$. The $w_{g,t}$ are written as:


$$w_{h,s} = \frac{\hat{\xi}_{h,s}}{\sum_{(g,t):D_{g,t}=1} \frac{N_{g,t}}{N_1} \hat{\xi}_{g,t}}$$


where $\hat{\xi_{g,t}}$ is the residual of a weighted TWFE of $D_{g,t}$, $t=1,\ldots T$, $g=1,\ldots G$ on time and group effects, where each $(g,t)$ is assigned weight $N_{g,t}$. Clearly, the weights satisfy $\sum_{(g,t): D_{g,t} = 1} \frac{N_{g,t}}{N_1} \cdot w_{g,t} =1$. However, as shown by the authors, __the weights may also be negative__. In contrast, weights of the (conditional) ATT:

$$\beta_{ATT|\mathbf{D}} = \sum_{(g,t): D_{g,t} = 1} \frac{N_{g,t}}{N_1}  \cdot \mathbb{E}[Y_{g,t}(1) - Y_{g,t}(0)|\mathbf{D}]$$

are always nonnegative. When there is no (conditional) treatment effect heterogeneity, both the conditional ATT and the conditional expectation of the TWFE are the same. Under heterogeneity, the authors provide examples where, even when __all__ individual treatment effects are positive, the expectation of the TWFE estimator is negative due to the behaviour of the weights. They also provide conditions for weights assigned to a treated cell $(g,t)$ to be negative. Weights in a given period are more likely to be negative if there is a higher fraction of treated groups in it. Similarly, weights in a given group are more likely to be negative if there are more treated periods in it. In staggered treatment adoption setups, this means treated groups by the end of the panel are more likely to receive negative weights. Units which have been exposed for longer periods are also more likely to receive negative weights. If there is treatment effect "build-up", this would imply, in a setting with variation in treatment timing, that TWFE is more likely to lead to misleading inference. 

Can we have a measure on how misleading basing analysis on $\hat{\beta}_{TWFE}$ could be? The authors derive minimal values for the conditional variance of the treatment effect heterogeneity, defined as:

$$\sigma^2_{TET} =  \sum_{(g,t): D_{g,t}} \frac{N_{g,t}}{N_1}( \mathbb{E}[Y_{g,t}(1) - Y_{g,t}(0)|\mathbf{D}] -  \beta_{ATT|\mathbf{D}} )^2$$
which are compatible with misleading estimation. In particular, they show that the smallest value $\sigma_{TET}$ compatible with both $\beta_{ATT|\mathbf{D}}$ and $\mathbb{E}[\hat{\beta}_{\text{TWFE}}|\mathbf{D}]$ equal to zero to be:

$$\underline{\sigma} = \frac{|\mathbb{E}[\hat{\beta}_{\text{TWFE}}|\mathbf{D}]|}{\sigma(w)}$$

where

$$\sigma^2(w) =\sum_{(g,t): D_{g,t}} \frac{N_{g,t}}{N_1}( w_{g,t} -  1 )^2$$

is the dispersion of weights. A small value of $\underline{\sigma}$ would indicate that even small amounts of heterogeneity are compatible with misleading results -- in the sense that the ATT and the expectation of the TWFE estimator have opposite signs. Note that we can consistently estimate $\underline{\sigma}$, which the authors recommend one to do.

Let's implement the above approach on data on wages for 51 US states. We will estimate the impact of having a higher state minimum wage, relatively to the federal level, on average real wages. We will use the package __fixest__ which provides a fast and easy way to include dummies in linear models.


```r
library(fixest)

dados = read.csv("base_state_final.csv")

#Restricting to >=1982, when panel is balanced
dados = dados[dados$year>=1982&dados$year<=2019,]

#Aggregating to yearly level -- min wage variation varies yearly
dados = aggregate(dados[,c("lrwage","high_min_wage")], by = list("year"=dados$year,"state" = dados$state), FUN = mean)

head(dados)
```

```
##   year state   lrwage high_min_wage
## 1 1982     1 5.331320             0
## 2 1983     1 5.371751             0
## 3 1984     1 5.340415             0
## 4 1985     1 5.391175             0
## 5 1986     1 5.381586             0
## 6 1987     1 5.356948             0
```

```r
twfe = feols(lrwage ~ high_min_wage|as.factor(year)+as.factor(state), data = dados)
summary(twfe, se = "cluster", cluster = "state")
```

```
## OLS estimation, Dep. Var.: lrwage
## Observations: 1,938 
## Fixed-effects: as.factor(year): 38,  as.factor(state): 51
## Standard-errors: Clustered (state) 
##                Estimate Std. Error   t value Pr(>|t|) 
## high_min_wage -0.007874   0.009237 -0.852471 0.398018 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## Log-likelihood: 3,328.34   Adj. R2: 0.88928 
##                          R2-Within: 0.00275
```

```r
#Computing minimal values for heterogeneity
beta_twfe = twfe$coefficients["high_min_wage"]

weight_reg = feols(high_min_wage~1|as.factor(year)+as.factor(state), data = dados)

dados$res = residuals(weight_reg)
N1 = sum(dados$high_min_wage)
dados$weights = dados$high_min_wage*dados$res/(sum(dados$high_min_wage*dados$res/N1))

summary(dados$weights[dados$high_min_wage!=0])
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -1.1667  0.4155  1.1434  1.0000  1.5423  3.3529
```

```r
minimal_sd = abs(beta_twfe)/sd(dados$weights[dados$high_min_wage!=0])

print(minimal_sd)
```

```
## high_min_wage 
##   0.009287294
```


So a std of treatment effect het. of around 0.93% would be compatible with "misleading inference". 


de Chaisemartin and D’Haultfeuille (2020) then proceed to show that, under some additional assumptions, we are able to estimate an ATE on the subpopulation of "changers" -- those groups which change treatment status -- at the periods of change; by running multiple period-by-period DID estimators and aggregating these. In particular, they consider the _estimand_:

$$\Delta : = \sum_{(g,t): t \geq 2, D_{g,t} \neq D_{g,t-1}} \frac{N_{g,t}}{N_\Delta} \mathbb{E}[Y_{g,t}(1) - Y_{g,t}(0)|\mathbf{D}]$$
where $N_{\Delta} = \sum_{(g,t): t \geq 2, D_{g,t} \neq D_{g,t-1}}  N_{g,t}$. In order to account for "leavers" -- groups for which $D_{g,t+1}<D_{g,t}$, the author imposes the additional exogeneity and parallel trend assumptions:

__Assumption 4 (more exogeneity)__ For all $g,t$, $\mathbb{E}[Y_{g,t+1}(1) - Y_{g,t}(1)|D_{g,1}\ldots D_{g,T}] = \mathbb{E}[Y_{g,t+1}(1) - Y_{g,t}(1)]$.


__Assumption 5 (more common trends) __ $\mathbb{E}[Y_{g,t+1}(1) - Y_{g,t}(1)]$ does not depend on $g$.

Note these assumptions limit dynamic effects (in particular, treatment effects do not depend on the history of past treatments). They are only required when leavers are present, which is __not__ the case in staggered designs.

Next, to be able to define their estimator, the authors replace the independence assumption with an assumption that implies that changers always have a valid "comparison" group. Since independence needs be dropped (it is incompatible with the guarantee of existence of comparison groups), the authors need to further assume that, for all $g,t$ $\mathbf{E}[Y_{g,t}(s)|\mathbf{D}] = \mathbf{E}[Y_{g,t}(s)|D_{g,1},\ldots D_{g,T}]$, $s \in \{0,1\}$. Under these new assumptions and hypotheses 2-5, they show their DID-then-pool estimator is (conditionally) unbiased for $\Delta$.

Below, we provide a function that computes their DID-then-pool estimator by constructing a "stacked" artificial dataset and running a single regression. This allows us to easily compute the covariances between treatment effects required in conducting inference. The user may pass a type of variance estimator, if she prefers to (default is to cluster at group level).


```r
library(sandwich)
#Function that computes the average treatment effect on the subpopulation of changers considered by de Chaisemartin and D’Haultfeuille (2020). Takes as arguments:

# outcome, treatment, time_var, group_var: character variables indicating the respective fields in the datase. Time_var should be date of numeric, allowing for ordering
# dataset: the dataframe to be used in estimation
# vcov.type: should be one of "standard", "hetero", "group", "time". For the latter two, it will cluster on tha specified level.
did_het <- function(outcome, treatment, time_var, group_var, dataset, vcov.type = "group")
{
  #Creating aggregate dataset
  data_agg = aggregate(cbind("Ngt" = 1, "outcome" = dataset[,outcome], "treatment" =  dataset[,treatment]), by = list("time" = dataset[,time_var], "group" = dataset[,group_var]), FUN = mean)
  
  time_vec = unique(data_agg$time)
  time_vec = time_vec[order(time_vec)]
  
  group_vec = unique(data_agg$group)
  
  
  data_reg = c()
  enter_weight = c()
  leave_weight = c()
  
  for(tt in 1:(length(time_vec)-1))
  {
    entering = sapply(group_vec, function(x){
      data_agg$treatment[(data_agg$time==time_vec[tt+1])&(data_agg$group==x)]>data_agg$treatment[(data_agg$time==time_vec[tt])&(data_agg$group==x)]
    })
    
    leaving = sapply(group_vec, function(x){
      data_agg$treatment[(data_agg$time==time_vec[tt+1])&(data_agg$group==x)]<data_agg$treatment[(data_agg$time==time_vec[tt])&(data_agg$group==x)]
    })
    
    never_treated = sapply(group_vec, function(x){
      (data_agg$treatment[(data_agg$time==time_vec[tt+1])&(data_agg$group==x)]==data_agg$treatment[(data_agg$time==time_vec[tt])&(data_agg$group==x)])&(data_agg$treatment[(data_agg$time==time_vec[tt+1])&(data_agg$group==x)]==0)
    })
    
   always_treated = sapply(group_vec, function(x){
      (data_agg$treatment[(data_agg$time==time_vec[tt+1])&(data_agg$group==x)]==data_agg$treatment[(data_agg$time==time_vec[tt])&(data_agg$group==x)])&(data_agg$treatment[(data_agg$time==time_vec[tt+1])&(data_agg$group==x)]==1)
    })
    
   if((sum(entering)>0)&(sum(never_treated)>0))
   {
     
         #Outcome 
        outt = sapply(group_vec[(entering+never_treated)>0], function(x){
      data_agg$outcome[data_agg$group==x&data_agg$time==time_vec[tt+1]]-data_agg$outcome[data_agg$group==x&data_agg$time==time_vec[tt]]
    })
        #Adding weights
        weight_est = sapply(group_vec[(entering+never_treated)>0], function(x){
      data_agg$Ngt[data_agg$group==x&data_agg$time==time_vec[tt+1]]
    })
    
         line_enter_data = data.frame("group" = group_vec[(entering+never_treated)>0], "time" = time_vec[tt+1], "char" = paste("enter_",tt,sep=""), "changer" = 1*(group_vec[(entering+never_treated)>0]%in%group_vec[entering]), "outcome" = outt, "weight_est"= weight_est,stringsAsFactors = F)
    
    #Adding outcome
     data_reg = rbind(data_reg, line_enter_data)
   }
   
   if((sum(leaving)>0)&(sum(always_treated)>0))
  {  
 #Outcome 
        outt = sapply(group_vec[(leaving+always_treated)>0], function(x){
      data_agg$outcome[data_agg$group==x&data_agg$time==time_vec[tt+1]]-data_agg$outcome[data_agg$group==x&data_agg$time==time_vec[tt]]
    })
        #Adding weights
        weight_est = sapply(group_vec[(leaving+always_treated)>0], function(x){
      data_agg$Ngt[data_agg$group==x&data_agg$time==time_vec[tt+1]]
    })
    
         line_enter_data = data.frame("group" = group_vec[(leaving+always_treated)>0], "time" = time_vec[tt+1], "char" = paste("leave_",tt,sep=""), "changer" = -1*(group_vec[(leaving+always_treated)>0]%in%group_vec[(leaving)]), "outcome" = outt, "weight_est"= weight_est, stringsAsFactors = F)
    
    #Adding outcome
     data_reg = rbind(data_reg, line_enter_data)

    
   }
  }
  
  mod = feols(outcome~changer:as.factor(char)|as.factor(char),data_reg, weights = data_reg$weight_est)
  
  if(vcov.type%in%c("standard","hetero"))
    	variance = vcov(mod, vcov.type) else variance = vcov(mod, "cluster", cluster = vcov.type)
  
  weight_mat = aggregate(x = data_reg$weight_est*data_reg$changer/sum(data_reg$weight_est*data_reg$changer), by = list("char" = paste("changer:as.factor(char)", data_reg$char,sep="")), FUN = sum)
  
  reorder = match(weight_mat$char, names(mod$coefficients))
  
  weight_mat = weight_mat[reorder,]
  
  point.est = mod$coefficients%*%weight_mat$x
  sd = sqrt(t(weight_mat$x)%*%variance%*%weight_mat$x)
  
  
  
  return(list("point.est" = point.est, "sd" = sd, "full_res" = mod, "data_2_reg" = data_reg))
  
}
```

Applying it to our data:


```r
resultados = did_het("lrwage","high_min_wage","year", "state", dados)

#Heterogeneous effect
summary(resultados$full_res, se = "cluster", cluster = resultados$data_2_reg[,"group"] )
```

```
## OLS estimation, Dep. Var.: outcome
## Observations: 705 
## Fixed-effects: as.factor(char): 23
## Standard-errors: Clustered 
##                                  Estimate Std. Error   t value    Pr(>|t|)    
## changer:as.factor(char)enter_11  0.003457   0.006562  0.526835 0.600638000    
## changer:as.factor(char)enter_13 -0.008024   0.011993 -0.669046 0.506544000    
## changer:as.factor(char)enter_17 -0.010675   0.013944 -0.765530 0.447555000    
## changer:as.factor(char)enter_19  0.024939   0.004244  5.876300 0.000000341 ***
## changer:as.factor(char)enter_21 -0.008264   0.005105 -1.618600 0.111822000    
## changer:as.factor(char)enter_22  0.008223   0.004514  1.821800 0.074477000 .  
## changer:as.factor(char)enter_23 -0.005402   0.017423 -0.310067 0.757798000    
## changer:as.factor(char)enter_24  0.001499   0.007683  0.195073 0.846126000    
## changer:as.factor(char)enter_25 -0.010488   0.018189 -0.576637 0.566772000    
## changer:as.factor(char)enter_28 -0.010069   0.008284 -1.215600 0.229863000    
## changer:as.factor(char)enter_29  0.005194   0.004224  1.229500 0.224633000    
## changer:as.factor(char)enter_30 -0.014091   0.004993 -2.822400 0.006824000 ** 
## changer:as.factor(char)enter_31 -0.014401   0.012091 -1.191000 0.239258000    
## changer:as.factor(char)enter_32 -0.021396   0.014400 -1.485900 0.143595000    
## changer:as.factor(char)enter_5   0.044766   0.009357  4.784400 0.000016000 ***
## changer:as.factor(char)enter_8  -0.002628   0.014361 -0.182983 0.855551000    
## changer:as.factor(char)leave_11  0.024774   0.009058  2.735100 0.008605000 ** 
## changer:as.factor(char)leave_14  0.005507   0.008661  0.635789 0.527813000    
## changer:as.factor(char)leave_15 -0.041603   0.014647 -2.840400 0.006500000 ** 
## changer:as.factor(char)leave_26  0.023927   0.015673  1.526600 0.133162000    
## changer:as.factor(char)leave_27 -0.002476   0.008248 -0.300166 0.765295000    
## changer:as.factor(char)leave_8   0.009957   0.005971  1.667700 0.101637000    
## changer:as.factor(char)leave_9  -0.005517   0.010405 -0.530287 0.598260000    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## Log-likelihood: 1,600.55   Adj. R2: 0.09715 
##                          R2-Within: 0.05323
```

```r
print(resultados$point.est)
```

```
##             [,1]
## [1,] 0.002888974
```

```r
print(resultados$sd)
```

```
##             [,1]
## [1,] 0.007730266
```




## References

De Chaisemartin, C., & d'Haultfoeuille, X. (2020). Two-way fixed effects estimators with heterogeneous treatment effects. American Economic Review, 110(9), 2964-96.

