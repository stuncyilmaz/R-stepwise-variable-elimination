this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source(paste( getwd(),"/backward_elimination.R",sep='') )
################## load some  data
filename="SENIC_data.csv"
df <- read.table(filename, sep=",", header=FALSE )
# create categorical binary variables for the column "region"
df$region <- factor(df$V9)
nlevels=nlevels(df$region)
df=cbind(df,model.matrix( ~ region - 1, data=df )[,1:nlevels-1])
df=df[ , -which(names(df) %in% c("V1","V9","region"))]

# find the model with minimum AIC
best_model=bswr(df)
cat("minimum AIC is: ", getAIC(best_model),"\n")
cat("The variables in the proposed model are: \n")
print(names(coef(best_model)))
