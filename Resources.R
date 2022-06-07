###########################################
## Get biomarker for Tabalumab response
## based on transcriptomic response
###########################################
##
## Functions

## Log Normalization
norm.log<-function(data){
  qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) { data[which(data <= 0)] <- NaN
  data <- log2(data) }
  return(data)
}

## Calculate median of expression by gene
calculate_medians = function(gene, ## Gene to annotate
                             genome, ## Table with gene-probe set conections
                             expressionMatrix){ ## Gene-expression data
  
  probe_ids = genome[genome$toGenes==gene,"fromGenes"] ## Select probes for each gene
  if (length(probe_ids)>1){
    res =  colMedians(as.matrix(expressionMatrix[probe_ids,])) ## Median of all probes 
  } else{
    res <- as.numeric(expressionMatrix[probe_ids,])  
  }
  return(res)
}

## Function to annotate genes from a gene-expression matrix
annotateGenes<-function(data,
                        toGenes='external_gene_name',
                        fromGenes='ensembl_gene_id'){
  
  require("biomaRt")
  require("parallel")
  
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  genome = getBM(attributes=c(toGenes,fromGenes),mart = ensembl)
  colnames(genome)<-c("toGenes","fromGenes")
  genome <- genome[genome$fromGenes != "",]
  genome = genome[genome$fromGenes %in% rownames(data),]
  data = data[genome$fromGenes,]
  finalGenes = unique(genome$toGenes)
  
  if(as.character(Sys.info()["sysname"])=="Windows"){
    nCores = 1
  }else{
    nCores = detectCores()
  }
  
  temp = mclapply(finalGenes,calculate_medians,genome = genome,expressionMatrix = data,mc.cores = nCores)
  temp = as.data.frame(do.call("rbind", temp))
  rownames(temp) = finalGenes
  colnames(temp) = colnames(data)
  
  return(temp)
}


#-------
## Chech installed packages
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)>0){
    cat(paste0("\nInstalling ", paste(new.pkg,collapse = ", ")))
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
    }
    bioconductor_packages <- BiocManager::available()
    bioconductor_packages <- bioconductor_packages[bioconductor_packages %in% new.pkg]
    if (length(bioconductor_packages) > 0){
      BiocManager::install(bioconductor_packages, dependencies = TRUE,ask=FALSE)
    }
    new.pkg = new.pkg[!(new.pkg %in% bioconductor_packages)]
    if (length(new.pkg)>0){
      install.packages(new.pkg, dependencies = TRUE)
    }
  }
  res <- lapply(pkg,load.packages)
} ## Character vector contains package names

#-------
## Load packages
load.packages <- function(pkg){
  cat(paste0("\nLoading ",pkg))
  suppressMessages(require(pkg,character.only = T))
} ## Character vector contains package names

##-------
## Fix gramatical symbols 
FixSymbol<-function(vect,
                    symbol=c("/"," "), 
                    toSymbol=c(".","")){
  
  for(i in 1:length(symbol)){
    vect<-gsub(pattern = symbol[i],replacement = toSymbol[i],vect)
  }
  return(vect)
}

##-------
## Split data in a balanced way between multiple classes
splitData<-function(class.vector,  ## Vector with class labels
                    p){            ## %, size of the subsets (in 0-1 scale)
  require("caret")
  res.vect<-as.numeric(createDataPartition(class.vector, p = p, list=F,times=1))
  return(res.vect)
}

##-------
## Select important variables using Lasso
LassoImp<-function(data){
  x_vars <- model.matrix(Group~. , data)[,-1]
  y_var <- data$Group; y_var<-ifelse(y_var=="No",0,1)
  lambda_seq <- 10^seq(2, -2, by = -.1)
  
  cv_output <- cv.glmnet(as.matrix(x_vars), y_var,
                         alpha = 1, lambda = lambda_seq, 
                         nfolds = 5)
  ##-------
  ## Select important variables using Lasso
  LassoImp<-function(data){
    x_vars <- model.matrix(Group~. , data)[,-1]
    y_var <- data$Group; y_var<-ifelse(y_var=="No",0,1)
    lambda_seq <- 10^seq(2, -2, by = -.1)
    
    cv_output <- cv.glmnet(as.matrix(x_vars), y_var,
                           alpha = 1, lambda = lambda_seq, 
                           nfolds = 5)
    
    # identifying best lamda
    best_lam <- cv_output$lambda.min
    
    lasso_best <- glmnet(as.matrix(x_vars), y_var, alpha = 1, lambda = best_lam)
    
    importance<-coef(lasso_best)
    importance<-data.frame(Genes=rownames(importance),Importance=as.numeric(importance))
    
    return(importance)
  }
  # identifying best lamda
  best_lam <- cv_output$lambda.min
  
  lasso_best <- glmnet(as.matrix(x_vars), y_var, alpha = 1, lambda = best_lam)
  
  importance<-coef(lasso_best)
  importance<-data.frame(Genes=rownames(importance),Importance=as.numeric(importance))
  
  return(importance)
}




##-------
## Get classification model
GetModel<-function(data,           ## Data with class to predict in the first columns (named as Group), features in columns and samples in rows
                   method="rf",    ## Method to be applied in classification
                   bootstrap=0.8,  ## Percent to divide training and test (0-1 scale)
                   kfold=10,       ## Either the number of folds or number of resampling iterations
                   repeatedCv=20,  ## For repeated k-fold cross-validation only: the number of complete sets of folds to compute
                   perm.bias=10,   ## Number of permutations in train and test selection
                   seed=1234,
                   accuracy.type="Balanced Accuracy", ## "Balanced Accuracy" or "Accuracy to generate the better model
                   featureSets=c(1,2,3,4,5,6,7,8,9,10,15,20,25)){ ## Subsets of features to calculate minimun number of features
  
  require("caret")
  require("stringr")
  require("doParallel")
  require("glmnet")
  
  set.seed(seed)
  
  ## Allow parallel
  clusters = makeCluster(detectCores()-1)
  registerDoParallel(clusters)
  
  ## Fix colnames from data
  colnames(data)<-FixSymbol(vect=colnames(data),symbol=c("/"," ","-"),toSymbol=c(".","","_"))
  
  fitControl <- trainControl(method = "repeatedcv",
                             number = kfold,
                             repeats = repeatedCv,
                             allowParallel = T,
                             verboseIter = T,
                             savePredictions = TRUE,
                             classProbs = T)
  
  ##--------------------------------------------- STEP 1
  ## Tuning parameters (select best parameters for the model)
  
  ## split TRAIN and TEST
  sel<-splitData(data$Group,p = bootstrap)
  TRAIN<-data[sel,]; TEST<-data[-sel,]
  
  
  cat(paste0("\nTuning parameters for ",method,sep=""," algorithm...\n"))
  
  switch (method,
    rf = { ## RANDOM FOREST

      # mtry: Number of features to consider at every split
      n<-10
      if(n>round(sqrt(ncol(data))+2)){ n<-round(sqrt(ncol(data))+2)}
      mtry<-unique(seq(from=2,to=round(sqrt(ncol(data))+2,digits=0),by=round((round(sqrt(ncol(data))+2))/n,digits=0)))
      # ntree: Number of trees in random forest
      ntree<-seq(from=200,to=2000,by=100)
      
      tunegrid<-expand.grid(mtry=mtry)
      
      res<-list()
      pb = txtProgressBar(min = 0, max = length(ntree), initial = 0)
      for(i in 1:length(ntree)){
        model<- invisible(train(Group ~., data = TRAIN,method="rf",
                                trControl = fitControl,tuneGrid=tunegrid,ntree=ntree[i]))
        res[[i]]<-model
        setTxtProgressBar(pb,i)
      }
      
      ## Search best parameters (high accuracy)
      accuracy<-0
      trees<-NULL
      importance<-NULL
      for(i in 1:length(res)){
        
        tmp<-res[[i]]$results
        tmp<-tmp[order(tmp$Accuracy,decreasing=T),]
        if(tmp$Accuracy[1]>accuracy){
          mtry<-tmp$mtry[1]
          trees<-ntree[i]
          accuracy<-tmp$Accuracy[1]
          
          imp<-varImp(res[[i]])
          features.imp<-as.data.frame(cbind(rownames(imp$importance),imp$importance))
          colnames(features.imp)[1]<-"Features"
          features.imp$Overall<-as.numeric(features.imp$Overall)
          features.imp<-features.imp[order(features.imp$Overall,decreasing=T),]
          features.imp$Features<-factor(x = features.imp$Features,levels = unique(features.imp$Features))
          
          importance<-as.character(features.imp$Features)
        }
      }
      ## Best mtry and ntree selected and feature importance assigned
      
    },
    knn = { ## K-NEAREST NEIGHBORS

     tunegrid.knn<-expand.grid(k = seq(2,50,2))
     
     model<- invisible(train(Group ~., data = TRAIN,method="knn",
                             trControl = fitControl,tuneGrid=tunegrid.knn,tuneLength=25,
                             preProcess = c("center","scale")))
     
     ## Search best parameters (high accuracy)
     tmp<-model$results
     tmp<-tmp[order(tmp$Accuracy,decreasing=T),]
     # k: Number of neigbors
     k=tmp$k[1]
     
     importance<-LassoImp(TRAIN)
     importance$Importance<-abs(importance$Importance)
     importance<-importance[order(importance$Importance,decreasing=T),]
     importance<-as.character(importance$Genes[-1])
     
    },
    svmLinear = {
      
      tunegrid.svm<-expand.grid(C=c(0.01,0.02,0.05,0.1,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3))
      
      model<- invisible(train(Group ~., data = TRAIN,method="svmLinear",
                              trControl = fitControl,tuneGrid=tunegrid.svm))
      
      ## Search best parameters (high accuracy)
      tmp<-model$results
      tmp<-tmp[order(tmp$Accuracy,decreasing=T),]
      # C
      C=tmp$C[1]

      importance<-LassoImp(TRAIN)
      importance$Importance<-abs(importance$Importance)
      importance<-importance[order(importance$Importance,decreasing=T),]
      importance<-as.character(importance$Genes[-1])
      
    },
  )
  cat("\nOptimal parameters selected\n")

  
  ##--------------------------------------------- STEP 2 
  ## Test number of features and TRAN-TEST selection bias
  
  cat("\nGetting model with minimal features...\n")
  
  splitN<-function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  sel<-as.numeric(createDataPartition(data$Group, p = bootstrap, list=F,times=perm.bias))
  perm.bias.vals<-splitN(sel,perm.bias)
  
  resTable<-data.frame(matrix(data=NA,nrow=length(perm.bias.vals),ncol=length(featureSets)))
  colnames(resTable)<-featureSets
  rownames(resTable)<-1:length(perm.bias.vals)
  
  
  for(i in 1:length(featureSets)){ ## feature number
    
    tmp<-data[,c("Group",as.character(importance[featureSets[1:i]]))]
    
    for(b in 1:length(perm.bias.vals)){ ## selection bias
      
      TRAIN<-tmp[as.numeric(perm.bias.vals[b][[1]]),]
      TEST<-tmp[-as.numeric(perm.bias.vals[b][[1]]),]
      
      switch(method,
             rf={
               model<- invisible(train(Group ~., data = TRAIN,method="rf",
                                       trControl = fitControl,tuneGrid=expand.grid(mtry=mtry),
                                       ntree=trees))
             },
             knn={
               
               model<- invisible(train(Group ~., data = TRAIN,method="knn",
                                       trControl = fitControl,tuneGrid=expand.grid(k=k),
                                       preProcess = c("center","scale")))
             },
             svmLinear={
               model<- invisible(train(Group ~., data = TRAIN,method="svmLinear",
                                       trControl = fitControl,tuneGrid=expand.grid(C=C)))
             })
      
      
      res<-predict(model, newdata = TEST,type = "prob")
      
      pred<-as.factor(ifelse(res[,1]>=0.5, colnames(res)[1], colnames(res)[2]))
      real<-as.factor(TEST$Group)
      cm.model<-confusionMatrix(pred,real)              
      if(accuracy.type=="Accuracy"){
        acc<-as.numeric(cm.model$overall[accuracy.type])
      }else if(accuracy.type=="Balanced Accuracy"){
        acc<-as.numeric(cm.model$byClass[accuracy.type])
      }
      resTable[b,i]<-acc
    }
  }
  
  ##--------------------------------------------- STEP 3 
  ## Get final model
  cat("Extracting Final Model...\n")
  
  nFeatures<-apply(resTable,2,mean)
  nFeatures<-as.numeric(names(nFeatures[order(nFeatures,decreasing = T)])[1])
  
  subset<-as.numeric(rownames(resTable)[ order(resTable[,as.character(nFeatures)],decreasing=T)[1]])
  
  tmp<-data[,c("Group",as.character(importance[1:nFeatures]))]
  features.selected<-as.character(importance[1:nFeatures])
  
  TRAIN<-tmp[as.numeric(perm.bias.vals[subset][[1]]),]
  TEST<-tmp[-as.numeric(perm.bias.vals[subset][[1]]),]
  
  switch(method,
         rf={
           fit.model<- invisible(train(Group ~., data = TRAIN,method="rf",
                                   trControl = fitControl,tuneGrid=expand.grid(mtry=mtry),ntree=trees))
         },
         knn={
           
           fit.model<- invisible(train(Group ~., data = TRAIN,method="knn",
                                   trControl = fitControl,tuneGrid=expand.grid(k=k),
                                   preProcess = c("center","scale")))
         },
         svmLinear={
           fit.model<- invisible(train(Group ~., data = TRAIN,method="svmLinear",
                                   trControl = fitControl,tuneGrid=expand.grid(C=C)))
         })
  
  
  res<-predict(fit.model, newdata = TEST,type = "prob")
  cm.model<-confusionMatrix(as.factor(ifelse(res[,1]>=0.5, colnames(res)[1], colnames(res)[2])),
                            as.factor(TEST$Group))   
  results<-list(fit.model,features.selected,cm.model,resTable)
  names(results)<-c("model","features","stats","selection.bias")
  
  stopCluster(clusters)
  registerDoSEQ()
  
  return(results)
}

##-------
## Get Prediction models for categorical variables
getML.cat<-function(training,
                    testing,
                    kfold=10,
                    repeats=10,
                    sel.feature=NULL,
                    set.seed=123456788){
  
  require("caret")
  require("mlbench")
  require("pROC")
  require("rpart")
  require("randomForest")
  require("nnet")
  require("caretEnsemble")
  require("MLeval")
  require("pROC")
  require("ROCR")
  
  ##------
  ## Feature selection on training set
  if(!is.null(sel.feature)){
    
    ## DESARROLLAR UN PAR DE METODOS, NO PARA ESTE CASO  
    
  }
  
  ##-------
  ## TrainControl
  ## Se divide el training en 10 partes aleatorias (9 train, 1 test)
  ## Se repite 10 veces para calcular hiperparametros optimos
  
  my_control <- trainControl(
    method="repeatedcv",
    number=kfold,
    savePredictions="final",
    classProbs=TRUE,
    index=createResample(training$group, kfold),
    repeats=repeats)
  
  ##-------
  ## GetModels ## MEJORA A HACER, INCLUIR MAS MODELOS Y SELECCIONAR DESDE LA FUNCION CUALES SE QUIEREN
  model_list <- caretList(
    group~., data=training,
    trControl=my_control,
    methodList=c("glm","lda"),
    tuneList=list(
      xgbTree=caretModelSpec(method="xgbTree", tuneGrid=expand.grid(max_depth = c(2, 3, 4, 5, 6, 8, 10),nrounds = 50,eta = c(0.01,0.05,0.1),gamma = c(0,1),colsample_bytree=c(0.1, 0.4), min_child_weight=c(1, 10), subsample=c(0.5, 1))),
      rf=caretModelSpec(method="rf", tuneGrid=data.frame(.mtry=seq(2,round(sqrt(ncol(training)-1)),1))),
      knn=caretModelSpec(method="knn", tuneGrid=data.frame(.k=seq(2,50,2))),
      svmLinear=caretModelSpec(method="svmLinear", tuneLength = 15),
      svmRadial=caretModelSpec(method="svmRadial", tuneLength = 15),
      nnet=caretModelSpec(method="nnet", tuneGrid=expand.grid(.size = seq(from = 1, to = 4, by = 1),.decay = seq(from = 0.1, to = 0.5, by = 0.1))),
      nb=caretModelSpec(method="nb", tuneGrid=expand.grid(fL=c(0,0.5,1.0), usekernel = TRUE, adjust=c(0.5,1.0))))
    #
  )
  
  
  ##-------
  ## Extract stats
  res<-as.data.frame(matrix(data=0,ncol=length(model_list),nrow=13))
  colnames(res)<-names(model_list)
  for(cm in 1:length(model_list)){
    res[,cm]<-modelStats(model.i = model_list[[cm]],test = testing)
  }
  rownames(res)<-c("Accuracy","AUC","Balanced Accuracy","Precision","Recall","F1",
                   "Prevalence","Detection Rate","Detection Prevalence",
                   "Sensitivity","Specificity","Pos Pred Value","Neg Pred Value")
  
  res<-res[,order(-res["AUC",], -res["Balanced Accuracy",])]
  
  
  ##-------
  ## Aggregate best two models
  cor.bestModels<-modelCor(resamples(model_list))
  cor.bestModels<-cor.bestModels[colnames(res),colnames(res)]
  tags<-toEnsemble(cor.bestModels)
  
  if(length(tags)>1){
    
    listM<-model_list[ifelse(names(model_list) %in% tags,T,F)]
    aggregatedModel <- caretEnsemble(
      listM,
      trControl=trainControl(
        number=10,
        verboseIter = TRUE
      ))
    
    ensemble<-modelStats(model.i = aggregatedModel,test = testing)
    res<-cbind(ensemble,res)
    res<-res[,order(-res["AUC",], -res["Balanced Accuracy",])]
    
    model_list[["ensemble"]]<-aggregatedModel
  }
  
  model_list<-model_list[colnames(res)]
  
  models.results<-list(res,model_list)
  names(models.results)<-c("stats","models")
  
  return(models.results)
  
}

##-------
## Get Stats for the models
modelStats<-function(model.i,test){
  
  p <- predict(model.i, newdata=testing,type="prob")
  
  if(class(p)=="numeric"){
    cm.model<-confusionMatrix(as.factor(ifelse(p<=0.5,"YES","NO")),as.factor(testing$group))
    pred<-data.frame("NO"=p,"YES"=1-p,"obs"=testing$group)
  }else{ ##data.frame
    cm.model<-confusionMatrix(as.factor(ifelse(p$NO<=0.5,"YES","NO")),as.factor(testing$group))
    pred<-data.frame("NO"=p$NO,"YES"=1-p$NO,"obs"=testing$group)
    
  }
  pred<-na.omit(pred)
  
  stats<-as.numeric(cm.model$overall["Accuracy"])
  
  x<-NA
  tryCatch({
    x <- evalm(pred)
    x<-x[[7]][[1]]["AUC-ROC","Score"]
  }, error=function(e){})
  stats<-c(stats,x)
  names(stats)<-c("Accuracy","AUC")
  
  stats<-c(stats,cm.model$byClass[c("Balanced Accuracy","Precision","Recall","F1",
                                    "Prevalence","Detection Rate","Detection Prevalence",
                                    "Sensitivity","Specificity","Pos Pred Value","Neg Pred Value")])
  return(stats)
}


##-------
## Get Candidates to Ensemble models
toEnsemble<-function(mcor,minCor=0.25){
  candidates<-colnames(mcor)[1]
  forbiden<-rownames(mcor)[ifelse(mcor[,1]<=minCor & mcor[,1]>=0,F,T)]
  
  for(nc in 1:ncol(mcor)){
    if(!colnames(mcor)[nc] %in% forbiden){ 
      candidates<-c(candidates,colnames(mcor)[nc])
      forbiden<-c(forbiden,rownames(mcor)[ifelse(mcor[,nc]<=minCor & mcor[,1]>=0,F,T)])
      forbiden<-unique(forbiden)
    }
  }
  return(unique(candidates))
}

##-------
## Get Stats for the models
modelStats<-function(model.i,test){
  
  p <- predict(model.i, newdata=test,type="prob")
  
  if(class(p)=="numeric"){
    cm.model<-confusionMatrix(as.factor(ifelse(p<=0.5,"YES","NO")),as.factor(test$group))
    pred<-data.frame("NO"=p,"YES"=1-p,"obs"=test$group)
  }else{ ##data.frame
    cm.model<-confusionMatrix(as.factor(ifelse(p$NO<=0.5,"YES","NO")),as.factor(test$group))
    pred<-data.frame("NO"=p$NO,"YES"=1-p$NO,"obs"=test$group)
    
  }
  pred<-na.omit(pred)
  
  stats<-as.numeric(cm.model$overall["Accuracy"])
  
  x<-NA
  tryCatch({
    x <- evalm(pred)
    x<-x[[7]][[1]]["AUC-ROC","Score"]
  }, error=function(e){})
  stats<-c(stats,x)
  names(stats)<-c("Accuracy","AUC")
  
  stats<-c(stats,cm.model$byClass[c("Balanced Accuracy","Precision","Recall","F1",
                                    "Prevalence","Detection Rate","Detection Prevalence",
                                    "Sensitivity","Specificity","Pos Pred Value","Neg Pred Value")])
  return(stats)
}

# https://topepo.github.io/caret/available-models.html
##-------
## Get Prediction models for categorical variables
getML.cat<-function(training,
                    testing,
                    kfold=10,
                    repeats=10,
                    sel.feature=NULL,
                    set.seed=123456788){
  
  require("caret")
  require("mlbench")
  require("pROC")
  require("rpart")
  require("randomForest")
  require("nnet")
  require("caretEnsemble")
  require("MLeval")
  require("pROC")
  require("ROCR")
  
  ##------
  ## Feature selection on training set
  if(!is.null(sel.feature)){
    
    ## DESARROLLAR UN PAR DE METODOS, NO PARA ESTE CASO  
    
  }
  
  ##-------
  ## TrainControl
  ## Se divide el training en 10 partes aleatorias (9 train, 1 test)
  ## Se repite 10 veces para calcular hiperparametros optimos
  
  my_control <- trainControl(
    method="repeatedcv",
    number=kfold,
    savePredictions="final",
    classProbs=TRUE,
    index=createResample(training$group, kfold),
    repeats=repeats)
  
  ##-------
  ## GetModels ## MEJORA A HACER, INCLUIR MAS MODELOS Y SELECCIONAR DESDE LA FUNCION CUALES SE QUIEREN
  model_list <- caretList(
    group~., data=training,
    trControl=my_control,
    methodList=c("glm","lda"),
    tuneList=list(
      xgbTree=caretModelSpec(method="xgbTree", tuneGrid=expand.grid(max_depth = c(2, 3, 4, 5, 6, 8, 10),nrounds = 50,eta = c(0.01,0.05,0.1),gamma = c(0,1),colsample_bytree=c(0.1, 0.4), min_child_weight=c(1, 10), subsample=c(0.5, 1))),
      rf=caretModelSpec(method="rf", tuneGrid=data.frame(.mtry=seq(2,round(sqrt(ncol(training)-1)),1))),
      knn=caretModelSpec(method="knn", tuneGrid=data.frame(.k=seq(2,50,2))),
      svmLinear=caretModelSpec(method="svmLinear", tuneLength = 15),
      svmRadial=caretModelSpec(method="svmRadial", tuneLength = 15),
      nnet=caretModelSpec(method="nnet", tuneGrid=expand.grid(.size = seq(from = 1, to = 4, by = 1),.decay = seq(from = 0.1, to = 0.5, by = 0.1))),
      nb=caretModelSpec(method="nb", tuneGrid=expand.grid(fL=c(0,0.5,1.0), usekernel = TRUE, adjust=c(0.5,1.0))))
    #
  )
  
  
  ##-------
  ## Extract stats
  res<-as.data.frame(matrix(data=0,ncol=length(model_list),nrow=13))
  colnames(res)<-names(model_list)
  for(cm in 1:length(model_list)){
    res[,cm]<-modelStats(model.i = model_list[[cm]],test = testing)
  }
  rownames(res)<-c("Accuracy","AUC","Balanced Accuracy","Precision","Recall","F1",
                   "Prevalence","Detection Rate","Detection Prevalence",
                   "Sensitivity","Specificity","Pos Pred Value","Neg Pred Value")
  
  res<-res[,order(-res["AUC",], -res["Balanced Accuracy",])]
  
  
  ##-------
  ## Aggregate best two models
  cor.bestModels<-modelCor(resamples(model_list))
  cor.bestModels<-cor.bestModels[colnames(res),colnames(res)]
  tags<-toEnsemble(cor.bestModels)
  
  if(length(tags)>1){
    
    listM<-model_list[ifelse(names(model_list) %in% tags,T,F)]
    aggregatedModel <- caretEnsemble(
      listM,
      trControl=trainControl(
        number=10,
        verboseIter = TRUE
      ))
    
    ensemble<-modelStats(model.i = aggregatedModel,test = testing)
    res<-cbind(ensemble,res)
    res<-res[,order(-res["AUC",], -res["Balanced Accuracy",])]
    
    model_list[["ensemble"]]<-aggregatedModel
  }
  
  model_list<-model_list[colnames(res)]
  
  models.results<-list(res,model_list)
  names(models.results)<-c("stats","models")
  
  return(models.results)
  
}

 ##--------

prioML<-function(lista,nameplot=NULL){
  
  require("ggplot2")
  require("reshape")
  
  # Get best models..................
  RES<-matrix(data=0,ncol=9,nrow=13)
  colnames(RES)<-c("rf","nb","xgbTree","svmLinear","svmRadial","nnet","glm","knn","lda")
  rownames(RES)<-rownames(lista[[1]]$stats)
  RES<-as.data.frame(RES)
  
  for(perm in 1:length(lista)){
    for(i in 1:nrow(RES)){
      for(j in 1:ncol(RES)){
        if(is.null(lista[[perm]]$stats[i,colnames(RES)[j]])==FALSE){
          RES[i,j]<-RES[i,j]+lista[[perm]]$stats[i,colnames(RES)[j]]
        }
      }
    }
  }
  
  RES<-RES/length(lista)
  
  
  RES$stats<-rownames(RES)
  #print(RES)
  RES<-melt(RES)
  
  if(is.null(nameplot)==FALSE){
    tiff(filename = nameplot,width = 1500,height = 1200,res = 300)
    p1<-ggplot(RES,aes(x=stats,y=value,group=variable,color=variable))+theme_bw()+
      geom_line()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      scale_color_manual(values=c("rf"="#FF6633","nb"="#990066",
                                  "xgbTree"="#3399CC","svmLinear"="#33FF99",
                                  "svmRadial"="#009900","nnet"="firebrick3",
                                  "glm"="#CC9900","knn"="#FFFF33","lda"="#330000"))+
      geom_hline(yintercept = 0.7,linetype="dashed",color="blue")+ylim(0,1)
    
    plot(p1)
    invisible(dev.off())
  }
  
  tmp<-RES[RES$stats=="AUC",]
  tmp<-tmp[order(tmp$value,decreasing=T),]
  
  #if(tmp$value[1]>=0.6){ ## Filtro de 0.6 AUC medio
  tmp<-as.character(tmp$variable[1])
  
  ## Get most hight model
  m<-data.frame(lista[[1]]$stats[,tmp])
  rownames(m)<-rownames(lista[[1]]$stats)
  
  for(i in 2:length(lista)){
    m<-cbind(m,data.frame(lista[[i]]$stats[,tmp]))
  }
  colnames(m)<-1:length(lista)
  
  m<-m[,order(-m[2,],-m[3,])]
  best<-as.numeric(as.character(colnames(m)[1]))
  
  model<-lista[[best]]$models[[tmp]]
  stats<-as.data.frame(m[,1])
  rownames(stats)<-rownames(m)
  colnames(stats)<-"stats"
  
  res<-list(model,stats)
  names(res)<-c("model","stats")
  
  return(res)
  
  #}else{
  
  #  return("NotSelected")
  #}
  
}













