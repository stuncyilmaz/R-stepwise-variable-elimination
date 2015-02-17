# R-stepwise-variable-elimination


The program takes as input an R data frame, and outputs a linear regression model with minimum Akaike information criterion (AIC). The first column in the data frame is used as the response variable. 

The program starts with all the variables in the data frame, and their 2-way interactions. At each step, a term in the model that reduces AIC the most is removed. When a single variable V1 is removed, all the interaction terms that have V1 (for example V1xV2 , V1xV3, …, V1xVn) are removed as well. Variable elimination continues as long as stepwise single variable elimination can reduce AIC.
At each iteration, the code attempts to add a variable to the model as well, in case a single variable or an interaction term can reduce AIC.

“test.r” provides example for how to use the function with “SENIC” data set. 
