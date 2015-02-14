# R-stepwise-variable-elimination


The program takes as input an R dataframe, and outputs a linear regression model with minimum Akaike information criterion (AIC). The first column in the dataframe is used as the response variable. 

The program starts with all the variables in the datframe, and their 2-way interactions. At each step, a term in the model is removed that contributes most to the reduction of AIC. If at a step a single variable v1 is removed, all the intercation terms that have V1 (such as V1xV2 , V1xV3, …, V1xVn) are removed as well. Variable elimination continues until stepwise elimination cannot reduce AIC.

“test.r” provides example for how to use the function with “SENIC” data set. 
