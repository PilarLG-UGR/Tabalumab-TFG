###########################################
## Get biomarker for Tabalumab response
## based on transcriptomic response
###########################################
##
## Workflow

prioML(lista=RESULTS,nameplot = "modelsTab.tiff")




##------------------------------------------------------------------------- STEP 0
## Set functions and environment 

## Define main work folder
pathMain<-"C:/Users/pilar/Desktop/tabalumab"
pathData<-"C:/Users/pilar/Desktop/tabalumab"

set.seed(1234567)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)

source(paste0(pathMain,sep="","/Resources.R"))

packages<-c("affy","oligo","limma","GEOquery","stringr","glmnet","biomaRt","matrixStats","parallel","car","caret","lars","glmnet")
check.packages(packages)
rm(packages)

setwd(pathMain)

##------------------------------------------------------------------------- STEP 1
## Load data 

load("PilarData.RData")

##------------------------------------------------------------------------- STEP 3
## Identify SLE-important genes
## Traditional gene-expression analysis are not recommended, because we want to select genes that contributing in
## pathogenesis also in small subgroups of patients (High heterogeneity within the disease)ue en controles

#-------
## Select genes more variables in SLE with respect to controls

data<-cbind(H,SLE)
clin<-clin.tab

SLEc<-SLE

## Healthy stats by gene
#H<-data[,ifelse(clin$State=="Normal",T,F)]
meanH<-apply(H,1,mean)
sdH<-apply(H,1,sd)

#clinSLE<-clin[ifelse(clin$State!="Normal",T,F),]
patients<-unique(clin$SubjectID)

selPats<-NULL
for(i in 1:length(patients)){
  
  tmp<-clin[clin$SubjectID==patients[i],]
  
  if("baseline" %in% tmp$Time & "week52" %in% tmp$Time){
    selPats<-c(selPats,patients[i])
  }
  
}

clin<-clin[clin$SubjectID %in% selPats,]
clin<-clin[clin$Time!="week16",]

clin<-clin[clin$Treatment!="None",]
clin$Treatment<-ifelse(clin$Treatment=="Placebo","Placebo","Tabalumab")

SLE<-SLE[,rownames(clin)]

originalSLE<-SLE

data<-cbind(H,SLE)

## Fist, we remove non-variable genes (variance less than a determinate quantile: 0.7, this value is empirically selected based
## on data distribution)
sdg<-apply(data,1,sd)
plot(density(sdg))
sdg<-ifelse(sdg<=quantile(sdg,0.7),F,T) 
table(sdg)
# select 9000 genes

## Healthy samples
#H<-data[,ifelse(clin$State=="Normal",T,F)] 
#SLE<-data[,ifelse(clin$State!="Normal",T,F)] 


#-------
## Parameters for plot significant and no significant gene-selection
Xplot<-apply(H,1,sd)
Yplot<-apply(SLE,1,sd)
Sign<-rep("#999999",length(Xplot))
names(Sign)<-names(Xplot)

SLE<-SLE[sdg,]
H<-H[sdg,]
rm(sdg)

#-------
## Statistical analisis: get most SLE-variant genes
RES<-matrix(data=0,ncol=5,nrow=nrow(H))
colnames(RES)<-c("Healthy_var","All_var","Pvalue","adj.Pvalue","Difference")
rownames(RES)<-rownames(H)
RES<-as.data.frame(RES)
RES[,1]<-apply(H,1,var)
RES[,2]<-apply(SLE,1,var)

norm=FALSE
pb = txtProgressBar(min = 0, max = nrow(data), initial = 0) 
for(i in 1:nrow(SLE)){
  
  # shapiro.test(as.numeric(H[1,])) #test de normalidad
  
  print(i)
  if(norm){
    setTxtProgressBar(pb,i)
    res<-bartlett.test(list(as.numeric(H[i,]),as.numeric(SLE[i,])))
    RES[i,3]<-res$p.value
  }else{
    setTxtProgressBar(pb,i)
    tmp<-matrix(data=0,ncol=2,nrow=ncol(H)+ncol(SLE))
    colnames(tmp)<-c("Values","Group")
    tmp<-as.data.frame(tmp)
    Values<-c(as.numeric(H[i,]),as.numeric(SLE[i,]))
    tmp$Values<-Values
    tmp$Group<-as.factor(c(rep("Healthy",ncol(H)),rep("All",ncol(SLE))))
    lev.res<-leveneTest(y=tmp$Values,group=tmp$Group,center="median")
    RES[i,3]<-lev.res$`Pr(>F)`[1]
  }
}

table(ifelse(RES$Pvalue<=0.05,T,F))

adj.Pvalue<-p.adjust(as.numeric(RES$Pvalue), method = "bonferroni", n = length(as.numeric(RES$Pvalue)))
RES$adj.Pvalue<-adj.Pvalue
table(ifelse(RES$adj.Pvalue<=0.01,T,F)) # 3129

RES$Difference<-RES$All_var-RES$Healthy_var


Sign[rownames(RES)[RES$adj.Pvalue<=0.01 & RES$Difference>0]]<-"#FF9933"
Sign[rownames(RES)[RES$adj.Pvalue<=0.01 & RES$Difference<0]]<-"#66CC00"

Color<-ifelse(RES$adj.Pvalue<=0.01 & RES$Difference>0,"#FF9933",
              ifelse(RES$adj.Pvalue<=0.01,"#66CC00","#999999"))



#plot(x=apply(H,1,sd),y=apply(data,1,sd),col=Color,pch=16,
#     xlab="Healthy variance",ylab="Healthy+SLE variance",main="Levenne test",
#     xlim=c(0,2.5),ylim=c(0,2.5))
#abline(col="#999999",coef = c(0,1))

#-------
## Plot of selected genes (orange)


tiff("GenesVariance.tiff",res = 300,width = 4.5, height = 4,units="in")
plot(x=Xplot,y=Yplot,col=Sign,pch=16,cex=0.6,
         xlab="Healthy variance",ylab="SLE variance",
         xlim=c(0,2.5),ylim=c(0,2.5))
abline(col="#999999",coef = c(0,1))
dev.off()

table(ifelse(RES$adj.Pvalue<=0.01 & RES$Difference>0,T,F)) 

## Select genes most variable in SLE (orange)
sel<-ifelse(RES$adj.Pvalue<=0.01 & RES$Difference>0,T,F) ## 3118 Genes

save.image("tmp.RData")

SLE<-SLE[sel,]
H<-H[sel,]

RES<-RES[sel,]

sel<-ifelse(RES$Difference>=0.5,T,F)
RES<-RES[sel,]
RES<-RES[order(RES$Difference,decreasing=T),]
length(intersect(rownames(SLE),rownames(RES)))

SLE<-SLE[rownames(RES),]
H<-H[rownames(RES),]

#rm(Sign,Xplot,Yplot,lev.res,RES,tmp,adj.Pvalue,Color,i,norm,sel,Values,pb)

rm(list=setdiff(ls(),c("SLE","H","clin","SLEc")))

save.image("geneFiltering.RData")


##------------------------------------------------------------------------- STEP 4
## Define responder and non-responsders patients to the treatment

## Option 1: Transcriptional response, distance to healthy controls

#load("geneFiltering.RData")

#-------## Option 1: Transcriptional response, distance to healthy controls

meanH<-apply(H,1,mean)
sdH<-apply(H,1,sd)
#-------
## Get distances
patients<-unique(clin$SubjectID)

distances<-list()
for(i in 1:length(patients)){
  
  tmp<-clin[clin$SubjectID==patients[i],]
  baseline<-rownames(tmp[tmp$Time=="baseline",])
  w52<-rownames(tmp[tmp$Time=="week52",])
  
  category<-as.character(tmp[w52,"Treatment"])
  x<-as.vector(SLE[,baseline]); names(x)<-rownames(SLE)
  baseline<-(x-meanH)  #/sdH
  x<-as.vector(SLE[,w52]); names(x)<-rownames(SLE)
  w52<-(x-meanH) #/sdH
  
  res<-list("category"=category,"baseline"=baseline,"w52"=w52)
  distances[[i]]<-res
}

rm(res,i,category,meanH,baseline,sdH,w52,x,tmp)
## "Placebo", "LY 120 mg Q4W", "LY 120 mg Q2W"

#-------
## Melt list
distances.matrix<-data.frame(matrix(ncol=3,nrow=0))
for(i in 1:length(distances)){
  
  tmp<-distances[[i]]
  category<-as.character(tmp$category)
  category<-gsub(pattern = " ",replacement = "",category)
  res<-NULL
  
  res<-c(res,category,sum(abs(as.numeric(tmp$baseline))),
                         sum(abs(as.numeric(tmp$w52))))
  distances.matrix<-rbind(distances.matrix,res)         
}
colnames(distances.matrix)<-c("Category","Baseline","w52")

distances.matrix<-cbind(patients,distances.matrix)
distances.matrix$Baseline<-as.numeric(distances.matrix$Baseline)
distances.matrix$w52<-as.numeric(distances.matrix$w52)

rm(i,tmp,category,patients)

#-------
## Delta between time and baseline
Delta<-distances.matrix$w52-distances.matrix$Baseline
## CUanto mas negativo, mayor respuesta

distances.matrix<-cbind(distances.matrix,Delta)

#-------
## Set threshold to consider non-random response (from placebo), less than 0.1 percentile of Delta
Placebo.threshold<-(-55)

plot(density(distances.matrix[ifelse(distances.matrix$Category=="Placebo",T,F),"Delta"]),col="black")
lines(density(distances.matrix[ifelse(distances.matrix$Category!="Placebo",T,F),"Delta"]),col="orange")
abline(v=-55,lty=2,col="grey")
#-------
## Identify responder and non-responder patients
Response<-rep("Not_selected",nrow(distances.matrix))

for(i in 1:nrow(distances.matrix)){
  if(distances.matrix[i,"Category"]!="Placebo"){
    if(distances.matrix[i,"Delta"]<=Placebo.threshold){
      Response[i]<-"Yes"
    }else{
      Response[i]<-"No"
    }
  }
}
distances.matrix<-cbind(distances.matrix,Response)

#-------
## Get samples id from responders and non-responders to construct the gene-expression matrix of baseline samples for caret
pats.responders<-NULL
pats.nonResponders<-NULL
for(i in 1:nrow(distances.matrix)){
  
  if(distances.matrix[i,"Response"]=="No"){
    pats.nonResponders<-c(pats.nonResponders,distances.matrix[i,"patients"])
  }
  if(distances.matrix[i,"Response"]=="Yes"){
    pats.responders<-c(pats.responders,distances.matrix[i,"patients"])
  }
}

response<-c(rep("YES",length(pats.responders)),rep("NO",length(pats.nonResponders)))
pats<-c(pats.responders,pats.nonResponders)

samplesID<-NULL
for(i in 1:length(pats)){
  tmp<-clin[ifelse(clin$SubjectID==pats[i],T,F),]
  tmp<-rownames(tmp[ifelse(tmp$Time=="baseline",T,F),])
  samplesID<-c(samplesID,tmp)
}

SLEc<-SLEc[,samplesID]

#rm(list=setdiff(ls(),c("SLEc","response")))
save.image("ForCells.RData")

SLEc<-data.frame("Genes"=rownames(originalSLE),originalSLE)

write.table(x = SLEc,file="matrix_forcells.txt",sep="\t",row.names = F,quote = F)
  
  
#Treatment<-ifelse(distances.matrix$Category!="Placebo","Tabalumab","Placebo")
#distances.matrix<-cbind(distances.matrix,Treatment)

#dev.off()
#tiff("DensityPlot.tiff",res = 300,width = 4.5, height = 3,units="in")
#p<-ggplot(distances.matrix,aes(x=Delta,color=Treatment))+geom_density()+
#  geom_vline(xintercept=Placebo.threshold, linetype="dashed", 
#             color = "darkgrey", size=0.65)+theme_bw()
#plot(p)
#dev.off()


#-------
## Gene-expression matrix
data.baseline<-SLE[,samplesID]
data.baseline<-t(data.baseline)
colnames(data.baseline)<-gsub(pattern = "-","_",colnames(data.baseline))
data.baseline<-as.data.frame(data.baseline)
data.baseline<-cbind(response,data.baseline)

rm(distances,Response,Delta,i,pats.nonResponders,pats.responders,pats,Placebo.threshold,res,response,samplesID,tmp)

colnames(data.baseline)[1]<-"group"

save.image("preModels.RData")


##------------------------------------------------------------------------- STEP 5
## Get models for Tabalumab response

n.perm.bias<-10
RESULTS<-list()
data<-data.baseline

list.perm.bias<-list()
for(pb in 1:n.perm.bias){
  list.perm.bias[[pb]]<-createDataPartition(y = data$group, p = 0.8, list = FALSE)
}

for(mdl in 1:length(list.perm.bias)){
  
  training <- data[ list.perm.bias[[mdl]],]
  testing <- data[-list.perm.bias[[mdl]],]
  
  res.i<-getML.cat(training = training,testing = testing,
                   kfold = 10,repeats = 30,sel.feature = NULL)
  
  RESULTS[[mdl]]<-res.i
  
}

save.image("TabaModels.RData")

prioML(lista=RESULTS,nameplot = "selModel.tiff")

tabModel<-RESULTS[[2]]$models$rf

#saveRDS(tabModel,"ModeloFinalRF.rds")

train<-RESULTS[[2]]$models$rf$trainingData
colnames(train)[1]<-"group"
sel<-ifelse(rownames(data) %in% rownames(train),F,T)
test<-data[sel,]

modelStats(model.i = tabModel,test = test)

imp<-varImp(tabModel)
imp<-imp$importance
sel<-order(imp$Overall,decreasing=T)

imp<-data.frame(rownames(imp)[sel],
                imp[sel,"Overall"])
colnames(imp)<-c("genes","importance")

imp$genes<-factor(x = imp$genes,levels = unique(imp$genes))

tiff("imp.tiff",width=7,height = 2.5,res=300,units="in")
p1<-ggplot(imp,aes(x=genes,y=importance))+geom_jitter(col="red")+theme_classic()+
  theme(text = element_text(size=5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(p1)
dev.off()

##



#######################################
## PCA
library(Rtsne)

data.b.f<-data.baseline[,c("Group",models$svmLinear$features)]

#x<-prcomp(data.b.f[,-1])  ## data aqui es tu matriz de metilación,
#D<-as.data.frame(x$x[,1:2])  ## Aqui vamos a sacar las 2 dimensiones que hemos calculado con prcomp
#D<-cbind(data.b.f$Group,D); colnames(D)[1]<-"Response"  ## Ahora lo que hacemos es añadirle una columna mas a D, donde tenemos los clusters a los que pertenece cada muestra,
##

ggplot(D,aes(x=PC1,y=PC2,color=Response)) + theme_classic() + geom_point(shape=16,size=2.5)+
  scale_color_manual(values=c("#88af88","#cc3300")) #+
  #stat_ellipse(geom="polygon", aes(fill = Response),
  #             alpha = 0.0, linetype=10,
  #             show.legend = FALSE,
  #             level = 0.9)
## V1 y V2 serian los nombres de columnas en D (no recuerdo como se llaman por defecto), en scale_color_manual puedes meterle tantos colores como clusters tengas, y así pones la figura como mas te guste

set.seed(123456)
Mx<-t(data.b.f[,-1])
x = Rtsne(t(Mx), check_duplicates=FALSE, pca=TRUE, perplexity=40, theta=0.5, dims=3) ## Aqui puedes jugar con los parámetros, en el paquete dicen cuales son los rangos optimos para cada uno
D<-x$Y; D<-as.data.frame(D)
D<-cbind(data.b.f$Group,D); ## Igual que antes, le añadimos los clusters a los que pertenece cada muestra
colnames(D)[1]<-"Response"

dev.off()
tiff("PCA_Responsne.tiff",res = 300,width = 4, height = 3,units="in")
p<-ggplot(D,aes(x=V1,y=V2,color=Response)) + theme_classic() + geom_point(shape=16,size=2.2)+
  scale_color_manual(values=c("#80ad8a","#e69f00")) #+
#stat_ellipse(geom="polygon", aes(fill = Response),
#             alpha = 0.0, linetype=10,
#             show.legend = FALSE,
#            level = 0.9)
plot(p)
dev.off()

rownames(distances.matrix)<-distances.matrix$patients

D<-cbind(as.numeric(as.character(distances.matrix[clinSLE[rownames(data.b.f),"SubjectID"],"Delta"])),D)
colnames(D)[1]<-"Delta"
library(viridis)

dev.off()
tiff("PCA_Responsne_Delta.tiff",res = 300,width = 4, height = 3,units="in")
p<-ggplot(D,aes(x=V1,y=V2,color=Delta)) + theme_classic() + geom_point(shape=16,size=2.2) +
  scale_color_viridis(discrete = FALSE, option = "D")
plot(p)
dev.off()

## Boxplots

ggplot(data.baseline,aes(x=Group,y=UTS2)) + theme_classic() + geom_boxplot()


#+
  #scale_color_manual(values=c("#88af88","#cc3300"))

## Normalized method


## Cell percentage methods













