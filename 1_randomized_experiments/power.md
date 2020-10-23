---
title: "Power Calculations via simulation"
author: "Luis Alvarez"
date: "8/20/2020"
output: 
  html_document:
    keep_md: true
---



## Power calculations via simulation

Power calculations are an important part of experimental design. In class, you have seen how to perform power analyses in simple settings under random sampling of individuals from the population and Bernoulli assignment of treatments. In this note, we will see how simulation can be used to analyse more complex settings and help design experiments.

We can use baseline data (prior to conducting an experiment) to assess the power of our methods. We will consider a setting where we have baseline student-level data on a sample of several schools in Recife; and would like to randomize a treatment at the school level. We will also analyse the impact of including student-level covariates in power. Since treatment is assigned at the school-level, standard errors should account for clustering. To conduct power analysis, we must be able to:

1) Simulate the sampling process from the population of interest.
2) Simulate a draw from the treatment assignment mechanism.

In the context of RCTs, (2) is known. We assume here that half of the schools will be randomly treated. As for (1), one should ask oneself what is the parameter of interest, as this will change what the sampling is. If the parameter of interest is the average treatment effect *in the sample*, $\tau_s = N^{-1}(\sum_{i=1}^n Y_i(1) - Y_i(0))$ then one shpuld view the sample as fixed and only consider (2) (Abadie et al.,2020). In the case interest lies is the average treatment effect in the population from which the sample was drawn, $\tau_p = \mathbb{E}[Y_i(1)- Y_i(0)]$ (expectations wrt the population distribution), then (1) must be considered. We will consider the second case, i.e. we are interested in the effect on the population of schools from which our sample was drawn. Since we observe all students of interest in the sample, it is natural to assume that the sampling process involved drawing schools from the population of schools in Recife (this indeed was the case). Such a sampling scheme makes a further case for clustering standard errors at the school level (Abadie et al., 2017).


Finally, one must consider different scenarios for the treatment effect. Most analyses consider homogeneous treatment effects, though heterogeneity seems more likely and will affect power. Let us consider the homogeneous case first.


```r
head(data)
```

```
##     escola aluno     logico mulher    idade
## 89       1    89 -0.6110555      1 10.30411
## 151      1   151 -1.0899623      0 10.64384
## 152      1   152 -0.1321487      1 11.06027
## 235      1   235  0.5063938      0 11.26301
## 281      1   281  1.6238429      1 11.04658
## 290      1   290  0.1871226      1 10.62466
```

```r
summary(data)
```

```
##      escola          aluno           logico             mulher      
##  Min.   : 1.00   Min.   :  1.0   Min.   :-1.40923   Min.   :0.0000  
##  1st Qu.: 8.00   1st Qu.:227.2   1st Qu.:-0.77069   1st Qu.:0.0000  
##  Median :14.00   Median :453.5   Median :-0.05233   Median :0.0000  
##  Mean   :15.19   Mean   :453.5   Mean   : 0.00000   Mean   :0.4735  
##  3rd Qu.:23.00   3rd Qu.:679.8   3rd Qu.: 0.66603   3rd Qu.:1.0000  
##  Max.   :30.00   Max.   :906.0   Max.   : 3.06056   Max.   :1.0000  
##      idade       
##  Min.   : 8.784  
##  1st Qu.:10.661  
##  Median :11.060  
##  Mean   :11.219  
##  3rd Qu.:11.402  
##  Max.   :17.005
```



```r
#Sets seed to allow for replication
set.seed(1234)

#We'll consider a fine grid of treatment effects
grid = seq(0, 0.5, 0.05)

#Loading package that calculates cluster robust SEs
library(sandwich)

#We'll consider three specifications.
#Spec 1 Regression of Outcome on treatment indicator
#Sepc 2 Includes sex as control
#Spec 3 Includes both sex and age as controls
spec1.mat = c()
spec2.mat = c()
spec3.mat = c()

#Creating list of schools to sample from
list.schools = unique(data$escola)

for(eff in grid)
{
  rej.col.1 = c()
  rej.col.2 = c()
  rej.col.3 = c()


  for(replications in 1:1000)
  {
    #Drawing a sample, with replacement, of schools from our data. The idea here is that our sample well approximates the population distribution.
    #and that the sampling scheme is close to random sampling from the population of interest. 
    sample.draw = sample(list.schools, replace = T)
    
    #Construct the artificial data from a draw. Note that repeated schools should be treated as distinct so as to mimic random sampling.
    list.new.data = lapply(1:length(sample.draw), function(x){
      extract = data[data$escola==sample.draw[x],]
      extract$escola = x
      return(extract)
    })
    
    #Concatenating data on lists
    data.artificial = do.call(rbind, list.new.data)
    
    #Next, we select that half of the schools will be treated
    treat.draw = sample(unique(data.artificial$escola),size = length(unique(data.artificial$escola))/2, replace = F)
    
    data.artificial$treat = 1*(data.artificial$escola %in% treat.draw)
    
    #Create outcome
    data.artificial$y = data.artificial$logico + data.artificial$treat*eff
    
    #Running models and storing whether we reject the null that effect is 0
    model1 = lm(y~treat, data = data.artificial)
    
    se = sqrt(diag(vcovCL(model1, cluster = data.artificial$escola)))
    
    rej.col.1 = c(rej.col.1, 1*(abs(model1$coefficients["treat"]/se["treat"]) > qnorm(0.975) ))
    
    
    model2 = lm(y~treat+mulher, data = data.artificial)
    
    se = sqrt(diag(vcovCL(model2, cluster = data.artificial$escola)))
    
    rej.col.2 = c(rej.col.2, 1*(abs(model2$coefficients["treat"]/se["treat"]) > qnorm(0.975) ))
    
    
    model3 = lm(y~treat+mulher+idade, data = data.artificial)
    
    se = sqrt(diag(vcovCL(model3, cluster = data.artificial$escola)))
    
    rej.col.3 = c(rej.col.3, 1*(abs(model3$coefficients["treat"]/se["treat"]) > qnorm(0.975) ))

    
  }
  
  spec1.mat = cbind(spec1.mat, rej.col.1)
  spec2.mat = cbind(spec2.mat, rej.col.2)
  spec3.mat = cbind(spec3.mat, rej.col.3)
  
}

  
tabela = cbind("Effect" =  grid, "Spec 1" = colMeans(spec1.mat), "Spec 2" = colMeans(spec2.mat), "Spec 3" = colMeans(spec3.mat))
rownames(tabela) = rep("",nrow(tabela))

print(tabela)
```

```
##  Effect Spec 1 Spec 2 Spec 3
##    0.00  0.080  0.078  0.072
##    0.05  0.086  0.085  0.080
##    0.10  0.127  0.127  0.136
##    0.15  0.189  0.191  0.212
##    0.20  0.293  0.295  0.323
##    0.25  0.337  0.338  0.387
##    0.30  0.489  0.492  0.537
##    0.35  0.598  0.597  0.641
##    0.40  0.690  0.693  0.751
##    0.45  0.806  0.810  0.851
##    0.50  0.853  0.850  0.887
```

What if we want to introduce heterogeneity? For example, individual effects follow $\tau_i \sim N(\tau, 2)$.



```r
#Sets seed to allow for replication
set.seed(1234)

#We'll consider a fine grid of treatment effects
grid = seq(0, 0.5, 0.05)

#Loading package that calculates cluster robust SEs
library(sandwich)

#We'll consider three specifications.
#Spec 1 Regression of Outcome on treatment indicator
#Sepc 2 Includes sex as control
#Spec 3 Includes both sex and age as controls
spec1.mat = c()
spec2.mat = c()
spec3.mat = c()

#Creating list of schools to sample from
list.schools = unique(data$escola)

for(eff in grid)
{
  rej.col.1 = c()
  rej.col.2 = c()
  rej.col.3 = c()


  for(replications in 1:1000)
  {
    #Drawing a sample, with replacement, of schools from our data. The idea here is that our sample well approximates the population distribution.
    #and that the sampling scheme is close to random sampling from the population of interest. 
    sample.draw = sample(list.schools, replace = T)
    
    #Construct the artificial data from a draw. Note that repeated schools should be treated as distinct so as to mimic random sampling.
    list.new.data = lapply(1:length(sample.draw), function(x){
      extract = data[data$escola==sample.draw[x],]
      extract$escola = x
      return(extract)
    })
    
    #Concatenating data on lists
    data.artificial = do.call(rbind, list.new.data)
    
    #Next, we select that half of the schools will be treated
    treat.draw = sample(unique(data.artificial$escola),size = length(unique(data.artificial$escola))/2, replace = F)
    
    data.artificial$treat = 1*(data.artificial$escola %in% treat.draw)
    
    #Create outcome
    data.artificial$y = data.artificial$logico + data.artificial$treat*rnorm(nrow(data.artificial), mean = eff, sd = sqrt(2))
    
    #Running models and storing whether we reject the null that effect is 0
    model1 = lm(y~treat, data = data.artificial)
    
    se = sqrt(diag(vcovCL(model1, cluster = data.artificial$escola)))
    
    rej.col.1 = c(rej.col.1, 1*(abs(model1$coefficients["treat"]/se["treat"]) > qnorm(0.975) ))
    
    
    model2 = lm(y~treat+mulher, data = data.artificial)
    
    se = sqrt(diag(vcovCL(model2, cluster = data.artificial$escola)))
    
    rej.col.2 = c(rej.col.2, 1*(abs(model2$coefficients["treat"]/se["treat"]) > qnorm(0.975) ))
    
    
    model3 = lm(y~treat+mulher+idade, data = data.artificial)
    
    se = sqrt(diag(vcovCL(model3, cluster = data.artificial$escola)))
    
    rej.col.3 = c(rej.col.3, 1*(abs(model3$coefficients["treat"]/se["treat"]) > qnorm(0.975) ))

    
  }
  
  spec1.mat = cbind(spec1.mat, rej.col.1)
  spec2.mat = cbind(spec2.mat, rej.col.2)
  spec3.mat = cbind(spec3.mat, rej.col.3)
  
}

  
tabela = cbind("Effect" =  grid, "Spec 1" = colMeans(spec1.mat), "Spec 2" = colMeans(spec2.mat), "Spec 3" = colMeans(spec3.mat))
rownames(tabela) = rep("",nrow(tabela))

print(tabela)
```

```
##  Effect Spec 1 Spec 2 Spec 3
##    0.00  0.079  0.075  0.063
##    0.05  0.092  0.095  0.087
##    0.10  0.108  0.110  0.114
##    0.15  0.177  0.182  0.198
##    0.20  0.242  0.241  0.262
##    0.25  0.359  0.355  0.384
##    0.30  0.394  0.397  0.444
##    0.35  0.558  0.564  0.619
##    0.40  0.644  0.644  0.672
##    0.45  0.746  0.743  0.780
##    0.50  0.797  0.799  0.843
```

## References

Abadie, A., Athey, S., Imbens, G. W., & Wooldridge, J. (2017). When should you adjust standard errors for clustering? (No. w24003). National Bureau of Economic Research.

Abadie, A., Athey, S., Imbens, G.W. and Wooldridge, J.M. (2020), Sampling‐Based versus Design‐Based Uncertainty in Regression Analysis. Econometrica, 88: 265-296. doi:10.3982/ECTA12675
