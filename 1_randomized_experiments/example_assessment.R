#Clears previously loaded objects
rm(list=ls())
setwd("~/Dropbox/USP/dout/Monitorias/Microeconometria")

#Loads data
data = read.csv("dados_follow.csv")

#Sets seed to allow for replication

#Package to compute CR SE
library(sandwich)

#Package to produce test tables with robust SEs 
library(lmtest)

#Running the regression of interest
modelo = lm(logico~tratamento+idade+mulher, data)

#Inference with CR standard errors
coeftest(modelo, vcov. = vcovCL, cluster = data$escola)

#We'll perform the assessment on the specification above, to test the size on the coefficient of tratamento

#Estimating model under the null
modelo_nulo = lm(logico~idade+mulher, data)

#Vec p-values
vec_p_values = c()

for(j in 1:1000)
{
  #Generating artificial outcome
  data$artificial_outcome = modelo_nulo$fitted.values + rnorm(nrow(data))
  
  modelo_sim = lm(artificial_outcome~tratamento+idade+mulher, data)
  
  tab = coeftest(modelo_sim, vcov. = vcovCL, cluster = data$escola)
  
  vec_p_values = c(vec_p_values, tab["tratamento",4])
  
}

#Table with assessment. Rejection rates at different significance levels
vec_reject =sapply(c(0.01,0.05,0.1),function(x){mean(vec_p_values<=x)})

table = cbind("Significance level" = paste(c(0.01,0.05,0.1)*100,"%",sep=""), "Rejection rate"= paste(round(100*vec_reject,digits=2),"%",sep=""))

print(table)


#We can plot the cdf of p-values against the 45-degree line. This gives us an idea if we underreject (below the 45-degree line) or overreject 
#(above the 45-degree line) the null
cdf_pval = function(x){mean(vec_p_values<=x)}
grid = seq(0,1,0.001)
plot(grid, sapply(grid,cdf_pval),type = "l", col = "red", ylab = "F(x)", xlab = "x")
lines(grid,grid,lty=2)
