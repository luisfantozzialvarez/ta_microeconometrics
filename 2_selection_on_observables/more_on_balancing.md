---
title: "Additional balancing tests"
output:
  html_document:
    df_print: paged
    keep_md: true
---




I'll present an additional function for you to report balancing given the propensity score. This function summarizes the t-stats across propensity score blocks, and includes two additional tests that aggregate information across blocks. The function presumes you have already constructed the blocks (e.g. using the block construction function previously provided) and left a sufficient number of observations within each block to perform analyses.

Let us prepare our data for analysis.


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

#Blocks
lalonde$blocks.stepwise = propensity.score.blocks(lalonde$treat, lalonde$ps.linear, K = length(c("age","educ","black","nodegr","re74","re75","u74","u75","married")))
```


The function below reports t-Tests for equality of means in the $B$ blocks. We also include, in the penultimate column, the t-test suggested by Imbens and Rubin in page 298 of the book, which aggregates t-statistics across blocks (similar in implementation to our subclassification estimator). Since this test may not have much power against the alternative that at least one block is unbalanced (a negative unbalance in a group may compensate a positive discrepancy in another group), we also include the p-value of an F-test of the null that differences in means across all blocks for a given covariate are 0. See page 298-300 of Imbens and Rubin for a discussion. With regards to the t-tests, the function allows you to impose homoskedastic standard errors, which may be a good option if there are only a few treated/control units within each block (and you think variances within each group should not differ much).


```r
#Package for Ftests
library("car")
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

```r
#blocks - vector with block assignment
#treat = vector with treatment indicator label
#covariates = vector with covariate labels for us to test
#data = dataframe with variables
#equal = should equal variances in treatment and control groups be assumed in the t-stests? Defaults to F
balancing.table.blocks <- function(blocks, treat, covariates, data, equal = F)
{
  block.value = unique(blocks[order(blocks)])
  
  #Removing obs discarded in trimming
  block.value = block.value[!is.na(block.value)]
  data = data[!is.na(blocks),]
  
  treat_vec  = data[,treat]
  blocks = blocks[!is.na(blocks)]

  
  #Balancing tests across B blocks
  mat.results = t(sapply(covariates, function(x){
    cov_vec = data[,x]
    
    vec_covs = sapply(block.value, function(j){
      cov_block_treat = cov_vec[treat_vec==1&blocks==j]
      cov_block_control = cov_vec[treat_vec==0&blocks==j]
      
      return(tryCatch(t.test(cov_block_treat, cov_block_control, var.equal = equal)$statistic, error = function(e){
        NA
      }))
    })
    return(vec_covs)
  }))
  
  #Estimates for each block
  est.mat  = t(sapply(covariates, function(x){
    cov_vec = data[,x]
    
    vec_covs = sapply(block.value, function(j){
      cov_block_treat = cov_vec[treat_vec==1&blocks==j]
      cov_block_control = cov_vec[treat_vec==0&blocks==j]
      
      return(mean(cov_block_treat) - mean(cov_block_control))
    })
    return(vec_covs)
  }))
  
  #Variance estimates for each block
  var.mat = t(sapply(covariates, function(x){
    cov_vec = data[,x]
    
    vec_covs = sapply(block.value, function(j){
      cov_block_treat = cov_vec[treat_vec==1&blocks==j]
      cov_block_control = cov_vec[treat_vec==0&blocks==j]
      
      return(tryCatch({(t.test(cov_block_treat, cov_block_control, var.equal = equal)$stderr)^2}, error = function(e){NA}))
    })
    return(vec_covs)
  }))
  
  #Weights for each block
  weight.mat =  t(sapply(covariates, function(x){
    cov_vec = data[,x]
    
    vec_covs = sapply(block.value, function(j){
      cov_block_treat = cov_vec[treat_vec==1&blocks==j]
      cov_block_control = cov_vec[treat_vec==0&blocks==j]
      
      return(sum(blocks==j)/nrow(data))})
    return(vec_covs)
  }))

  #Adding names to table
  rownames(mat.results) = covariates
  colnames(mat.results) = paste("Block", block.value)
  
  #constructing aggregate (across blocks) t-stat
  est.agg = rowSums(weight.mat*est.mat)
  var.agg = rowSums((weight.mat^2)*var.mat)
  
  t.agg = est.agg/sqrt(var.agg)
  
  #Adding aggregaye statistic to results
  mat.results = cbind(mat.results, "Aggregate t-stat" = t.agg)
  
  data$block.variable.ftest = blocks
  #F-tests now
  f.pvalues = sapply(covariates, function(x){
    formula = paste(x, "~ -1 + as.factor(block.variable.ftest) + as.factor(block.variable.ftest):", treat,sep="")
    unres = lm(formula, data = data)
    coef.label = names(unres$coefficients)[grepl(treat, names(unres$coefficients))]
    
    tryCatch({linearHypothesis(unres, paste(coef.label, "= 0"))$`Pr(>F)`[length(coef.label)]}, error = function(e){NA})
  })
  
  mat.results = cbind(mat.results, "P-values F-test" = f.pvalues)
  
  return(mat.results)
}
```

Applying to our data


```r
balancing.table.blocks(lalonde$blocks.stepwise, "treat",c("age","educ","black","nodegr","re74","re75","u74","u75","married"), lalonde)
```

```
##             Block 1    Block 2 Aggregate t-stat P-values F-test
## age      0.42435910  0.3330083       0.53655212       0.8700495
## educ     0.79001308  0.4942682       0.77544246       0.7395194
## black   -0.03827226 -0.1110559      -0.10083693       0.9937161
## nodegr           NA -0.9273706               NA       0.4283150
## re74    -0.03484484  1.1041041       0.73778812       0.5429341
## re75    -0.19810286  1.3240043       0.71462152       0.4239142
## u74      0.42062481 -0.8145360      -0.39806914       0.6124851
## u75      0.34800329 -0.3067021      -0.06340014       0.9048681
## married  1.04064293  0.4340460       1.05061727       0.5176717
```

## References 

### Book

Imbens, G., & Rubin, D. (2015). Causal Inference for Statistics, Social, and Biomedical Sciences: An Introduction. Cambridge: Cambridge University Press. doi:10.1017/CBO9781139025751
