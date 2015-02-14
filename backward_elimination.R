# I chose AIC as the criterion because it panalizes overfitting with the p-parameter. It provides a good 
# balance between overfitting and complexity, and the goodness of fit measured by MSE.

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
print(getwd())
print(paste( getwd(),"/backward_elimination.R",sep='') )
rm(list=ls())

############# get the AIC, p-1 # of variables
getAIC = function(fit){
  p=summary(fit)$df[1]
  n=p+summary(fit)$df[2]
  SSE=anova(fit)[["Sum Sq"]][p]
  AIC=n*log(SSE)-n*log(n)+2*p
  AIC
}
##### extract variables
extractVariables= function(df){
  
  n=nrow(df)
  col_names <- colnames(df)
  all_variables_single=col_names[2:ncol(df)]
  all_var_length=length(all_variables_single)
  
  all_interactions <- matrix(data=0,nrow=all_var_length,ncol=all_var_length)
  for (i in seq(1,all_var_length-1, by=1)){
    for (j in seq(i+1,all_var_length, by=1)){
      a=paste(c(all_variables_single[i],all_variables_single[j]),collapse="*")
      all_interactions[i,j]=a
    }
  }
  
  ####### reshape interactions matrix to vector, call it doubles
  doubles=all_interactions
  dim(doubles) <- t(c(1, dim(doubles)[1]**2))
  doubles=doubles[doubles!=0]
  #singe and doubles
  combined_variables=c(all_variables_single,doubles)
  formula_variables <- paste(combined_variables,collapse="+")
  response_str=paste(col_names[1],"~")
  Formula <- formula(paste(response_str, formula_variables))
  all_mlr <- lm(Formula,df)
  minAIC=getAIC(all_mlr) # same as extractAIC(fit)
  
  list(all_interactions=all_interactions,all_variables_single=all_variables_single,col_names=col_names,
       response_str=response_str,all_mlr=all_mlr,minAIC=minAIC)
}
#####################
removeInteraction = function(interactions,i,j){
  interactions[i,j]=0
  interactions
}
## add interaction if corresponding single terms exist
addInteraction = function(interactions,i,j,singles=all_variables_single){
  if ((singles[i]!=0) & (singles[j]!=0) )   {
     a=paste(c(singles[i],singles[j]),collapse="*")
     interactions[i,j]=a  
  }
  interactions
}
# remove single and the corresponding interactions
removeSingle = function(variables_single,i,interactions=all_interactions){
  variables_single[i]=0
  interactions[i,]=0
  interactions[,i]=0
  list(variables_single=variables_single,interactions=interactions)
}
addSingle = function(variables_single,i){
  if (variables_single[i]==0){
    variables_single[i]=all_variables_single[i]
  }
  variables_single
}
#### Step for removing from single terms, finds the best variable to remove
StepSingleRemove = function(minAIC,current_interactions,current_singles,bestMlr,response=response_str){
  indices=which(current_singles!=0,arr.ind = T)
  var_index=0
  
  for (k in seq(1,length(indices), by=1)){
    i=indices[k]
    
    temp=removeSingle(current_singles,i,current_interactions)
    
    doubles=temp$interactions
    dim(doubles) <- t(c(1, dim(doubles)[1]**2))
    doubles=doubles[doubles!=0]
    
    temp_single=temp$variables_single[temp$variables_single!=0]
    
    combined_variables=c(temp_single,doubles)
    formula_variables <- paste(combined_variables,collapse="+")
    Formula <- formula(paste(response, formula_variables))
    mlr <- lm(Formula,df)
    myAIC=getAIC(mlr)
    
    #print(myAIC)
    if (myAIC<minAIC){
      var_index=i
      minAIC=myAIC
      bestMlr=mlr 
    }
    
  }
  list(var_index=var_index,minAIC=minAIC,bestMlr=bestMlr)
}

#### Step for removing from interaction terms, finds the best interaction variable to remove
StepInteractionRemove = function(minAIC,current_interactions,current_singles,bestMlr,response=response_str){
  indices=which(current_interactions!=0,arr.ind = T)
  myrow=0
  mycol=0
  if (!is.numeric(nrow(indices))){
    return(list(myrow=myrow,mycols=mycol,minAIC=minAIC,bestMlr=bestMlr))
  }
  for (k in seq(1,nrow(indices), by=1)){
    i=indices[k,][1]
    j=indices[k,][2]
    temp=removeInteraction(current_interactions,i,j)
    
    doubles=temp
    dim(doubles) <- t(c(1, dim(doubles)[1]**2))
    doubles=doubles[doubles!=0]
    
    current_singles=current_singles[current_singles!=0]
    combined_variables=c(current_singles,doubles)
    formula_variables <- paste(combined_variables,collapse="+")
    Formula <- formula(paste(response, formula_variables))
    mlr <- lm(Formula,df)
    myAIC=getAIC(mlr)
    

    #print(myAIC)
    
    if (myAIC<minAIC){
      myrow=i
      mycol=j
      minAIC=myAIC
      bestMlr=mlr 
    }
    
  }
  list(myrow=myrow,mycols=mycol,minAIC=minAIC,bestMlr=bestMlr)
}

# the step decides whether to remove either a single parameter or an interaction parameter
removeStep = function(minAIC,current_interactions,current_singles,removedSingles,removedDoubles,bestMlr,response_str){
  
  myStop=TRUE
  doubleGood=FALSE
  singleGood=FALSE
  
  doubles=current_interactions
  dim(doubles) <- t(c(1, dim(doubles)[1]**2))
  doubles=doubles[doubles!=0]
  
  if (length(doubles)>1){
    doubleStep=StepInteractionRemove(minAIC,current_interactions,current_singles,bestMlr,response_str)
    if (doubleStep$minAIC<minAIC){
      minAIC=doubleStep$minAIC
      doubleGood=TRUE}
  }
  
  
  if (length(current_singles)>1){
    singleStep=StepSingleRemove(minAIC,current_interactions,current_singles,bestMlr,response_str)
    if (singleStep$minAIC<minAIC){
      minAIC=singleStep$minAIC
      singleGood=TRUE
      doubleGood=FALSE}
  }


  if (doubleGood) {
    i=doubleStep$myrow
    j=doubleStep$mycol
    current_interactions=removeInteraction(current_interactions,i,j)
    bestMlr=doubleStep$bestMlr
    minAIC=doubleStep$minAIC
    removedDoubles[i,j]=TRUE
    myStop=FALSE
  }
  if (singleGood){
    
    i=singleStep$var_index
    removedSingles[i]=TRUE
    removedDoubles[i,]=TRUE
    removedDoubles[,i]=TRUE
    temp=removeSingle(current_singles,i,current_interactions)
    current_interactions=temp$interactions
    current_singles=temp$variables_single
    bestMlr=singleStep$bestMlr
    minAIC=singleStep$minAIC
    myStop=FALSE
    
  }
  
  list(bestMlr=bestMlr,minAIC=minAIC,current_singles=current_singles,current_interactions=current_interactions,
       removedSingles=removedSingles,removedDoubles=removedDoubles,myStop=myStop)
}



#assume there are elements in removedSingles
#step for adding a single parameter
StepSingleAdd = function(minAIC,removedSingles,current_singles,current_interactions,bestMlr,response=response_str){
  indices=which(removedSingles==TRUE,arr.ind = T)
  var_index=0
  
  for (k in seq(1,length(indices), by=1)){
    i=indices[k]

    doubles=current_interactions
    dim(doubles) <- t(c(1, dim(doubles)[1]**2))
    doubles=doubles[doubles!=0]
    
    variables_single=addSingle(current_singles,i)   
    temp_single=variables_single[variables_single!=0]
    
    combined_variables=c(temp_single,doubles)
    formula_variables <- paste(combined_variables,collapse="+")
    Formula <- formula(paste(response, formula_variables))
    mlr <- lm(Formula,df)
    myAIC=getAIC(mlr)
    
    #print(myAIC)
    if (myAIC<minAIC){
      var_index=i
      minAIC=myAIC
      bestMlr=mlr 
    }
    
  }
  list(var_index=var_index,minAIC=minAIC,bestMlr=bestMlr)
}

#assume there are removed interactions
# step for adding an interaction parameter
StepInteractionAdd = function(minAIC,removedDoubles,current_interactions,current_singles,bestMlr,response=response_str){
  indices=which(removedDoubles==TRUE,arr.ind = T)
  indices=indices[indices[,2]>indices[,1],]
  myrow=0
  mycol=0
  if (!is.numeric(nrow(indices))){
    return(list(myrow=myrow,mycols=mycol,minAIC=minAIC,bestMlr=bestMlr))
  }
  for (k in seq(1,nrow(indices), by=1)){
    i=indices[k,][1]
    j=indices[k,][2]
      temp=addInteraction(current_interactions,i,j,current_singles)
      
      doubles=temp
      dim(doubles) <- t(c(1, dim(doubles)[1]**2))
      doubles=doubles[doubles!=0]
      
      singles=current_singles[current_singles!=0]
      combined_variables=c(singles,doubles)
      formula_variables <- paste(combined_variables,collapse="+")
      Formula <- formula(paste(response, formula_variables))
      mlr <- lm(Formula,df)
      myAIC=getAIC(mlr)
      
      
      #print(myAIC)
      
      if (myAIC<minAIC){
        myrow=i
        mycol=j
        minAIC=myAIC
        bestMlr=mlr 
      }
      
    }
  list(myrow=myrow,mycols=mycol,minAIC=minAIC,bestMlr=bestMlr)
}

# the step decides whether to add a single variable or an interaction variable
addStep = function(minAIC,current_interactions,current_singles,removedSingles,removedDoubles,bestMlr,response_str){
  
  myStop=TRUE
  doubleGood=FALSE
  singleGood=FALSE
  
  doubles=current_interactions
  dim(doubles) <- t(c(1, dim(doubles)[1]**2))
  doubles=doubles[doubles!=0]
  
  
  indices=which(removedDoubles==TRUE,arr.ind = T)
  indices=indices[indices[,2]>indices[,1],]
  if (!is.numeric(nrow(indices))){
    return( list(bestMlr=bestMlr,minAIC=minAIC,current_singles=current_singles,current_interactions=current_interactions,
                 removedSingles=removedSingles,removedDoubles=removedDoubles,myStop=myStop))
  }
  
  if (nrow(indices)!=0){

    doubleStep=StepInteractionAdd(minAIC,removedDoubles,current_interactions,current_singles,bestMlr,response_str)
    if (doubleStep$minAIC<minAIC){
      minAIC=doubleStep$minAIC
      doubleGood=TRUE}
  }
  
  indices=which(removedSingles==TRUE,arr.ind = T)
  if (sum(indices)>0){
    singleStep=StepSingleAdd(minAIC,removedSingles,current_singles,current_interactions,bestMlr,response_str)
    if (singleStep$minAIC<minAIC){
      minAIC=singleStep$minAIC
      singleGood=TRUE
      doubleGood=FALSE}
  }

  if (doubleGood){
    i=doubleStep$myrow
    j=doubleStep$mycol
    current_interactions=addInteraction(current_interactions,i,j,current_singles)
    bestMlr=doubleStep$bestMlr
    minAIC=doubleStep$minAIC
    removedDoubles[i,j]=FALSE
    myStop=FALSE
  }
  if (singleGood ){
    i=singleStep$var_index
    removedSingles[i]=FALSE
    current_singles=addSingle(current_singles,i)
    bestMlr=singleStep$bestMlr
    minAIC=singleStep$minAIC
    myStop=FALSE
  }
  
  list(bestMlr=bestMlr,minAIC=minAIC,current_singles=current_singles,current_interactions=current_interactions,
       removedSingles=removedSingles,removedDoubles=removedDoubles,myStop=myStop)
}


bswr=function(df){

    n <<- nrow(df)
    ##################### Extract Variable Names
    
    if (ncol(df)<=2){
      col_names <- colnames(df)
      response_str=paste(col_names[1],"~")
      Formula <- formula(paste(response_str, col_names[2]))
      bestMlr<- lm(Formula,df)
      
      return(bestMlr)
    }
    a=extractVariables(df)
    
    all_interactions <<- a$all_interactions
    all_variables_single <<-a$all_variables_single
    col_names <<- a$col_names
    ####################### Create full model
    response_str <<- a$response_str
    bestMlr <<- a$all_mlr
    minAIC=a$minAIC

    all_var_length=length(all_variables_single)
    removedSingles=rep(FALSE, all_var_length)
    removedDoubles=matrix(data=FALSE,nrow=all_var_length,ncol=all_var_length)
    
    current_singles=all_variables_single
    current_interactions=all_interactions
    
    while (summary(bestMlr)$df[1]>2){
      #remove one variable if it decreases AIC
      a=removeStep(minAIC,current_interactions,current_singles,removedSingles,removedDoubles,bestMlr,response_str)
      #print(minAIC)
      if (a$myStop){
                    break}
    
      minAIC=a$minAIC
      current_singles=a$current_singles
      current_interactions=a$current_interactions
      removedSingles=a$removedSingles
      removedDoubles=a$removedDoubles
      bestMlr=a$bestMlr
      
     
      #print(summary(bestMlr)$df)
      #print(length(current_singles[current_singles!=0]))
      
      # add one variable if it decreases AIC
      b=addStep(minAIC,current_interactions,current_singles,removedSingles,removedDoubles,bestMlr,response_str)
      if (!b$myStop){
      
        minAIC=b$minAIC
        current_singles=b$current_singles
        current_interactions=b$current_interactions
        removedSingles=b$removedSingles
        removedDoubles=b$removedDoubles
        bestMlr=b$bestMlr

      }
 
    }
    return(bestMlr)
}

