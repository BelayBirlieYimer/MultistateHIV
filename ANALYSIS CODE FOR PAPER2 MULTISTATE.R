data1=read.csv("D:\\PHDfile\\PHDproject1\\Multistatemodeling\\Data\\prepartionforD4TAZTTDFanalysis.csv")

head(data1)
kk1<-as.data.frame(tapply(data1$New.end.time, data1$card.number,max))
names(kk1)<-c("stoptime")

kk1=read.csv("D:\\PHDfile\\PHDproject1\\Multistatemodeling\\Data\\stoptime.csv" )

hh3=merge(data1, kk1, by = "card.number", sort = FALSE)
head(hh3)


write.csv(hh3, file="D:\\PHDfile\\PHDproject1\\Multistatemodeling\\Data\\actuald4tazttdfdata.csv" )


#----------------changing data to multistate format 
data1=read.csv("D:\\PHDfile\\PHDproject1\\Multistatemodeling\\Data\\dataforpaper2\\d4tazttdfunchangedformultistate.csv" )


head(data1)
tail(data1)
detach(data1)
attach(data1)
pos.ch <- c("d4T","AZT","TDF")
uni<-unique(Card.number)
frq <- table(Card.number)
rep.card.no <- rep(uni,times=frq*2)
length(rep.card.no)
rep.start <- rep(Month.on.ART,each=2) 

endind1 <-(1:length(Endtime))[Month.on.ART!=0]
endtime1 <- Endtime
endtime1[endind1-1] <- Month.on.ART[endind1]
mix.start.end <- endtime1

rep.for.end <- rep(mix.start.end,each=2)

rep.for.from <- rep( Current.ARV.regimen, each=2)

exclude.current <- function(x=Current.ARV.regimen,y=pos.ch){
  y[!(x==y)]
}

for.to <- typeof(as.vector((sapply(Current.ARV.regimen,exclude.current,simplify=TRUE))))
for.to <- as.vector((sapply(Current.ARV.regimen,exclude.current,simplify=TRUE)))

rm.app <- function(x){
  no <- c(rep("no",2))
  ap <- append(x,no)
  ap[-(1:2)]
}
list.rff <- tapply(as.character(rep.for.from),rep.card.no,c)

for.sta <-sapply(list.rff,rm.app,simplify=TRUE)

for.sta.vec <- for.sta[[1]]
for(i in 2:length(for.sta)){
  for.sta.vec <- append(for.sta.vec,for.sta[[i]])
}

status <- as.numeric(for.sta.vec==as.character( for.to))

df.belaya<-data.frame(card.number=rep.card.no,Month.on.ART=rep.start,
                      Endtime=rep.for.end,From.regimen=rep.for.from,To.regimen=for.to,Status=status)
View(df.belaya)

df.belaya$trans <- 0
df.belaya$trans[df.belaya$From.regimen=="d4T" &df.belaya$To.regimen== "AZT"] <- 1
df.belaya$trans[df.belaya$From.regimen=="d4T" &df.belaya$To.regimen== "TDF"] <- 2
df.belaya$trans[df.belaya$From.regimen=="AZT" &df.belaya$To.regimen== "d4T"] <- 3
df.belaya$trans[df.belaya$From.regimen=="AZT" &df.belaya$To.regimen== "TDF"] <- 4
df.belaya$trans[df.belaya$From.regimen=="TDF" &df.belaya$To.regimen== "d4T"] <- 5
df.belaya$trans[df.belaya$From.regimen=="TDF" &df.belaya$To.regimen== "AZT"] <- 6

write.csv(df.belaya,"D:\\PHDfile\\PHDproject1\\Multistatemodeling\\Data\\dataforpaper2\\d4tazttdfchangedformultistate.csv")

# ---------- first do the non-parameteric analysis -------------------
library(survival)
library(mstate)

head(df.belaya)
names(df.belaya)<-c("id","Tstart","Tstop", "from","to","Status", "trans")

tmat1 <- transMat(x = list(c(2,3), c(1,3), c()),
                  names=c("d4T", "AZT", "TDF"))
tmat1

attr(df.belaya, "trans") <- tmat1
class(df.belaya) <- c("msdata","data.frame")

# ++++++++++++++++++++++ Observed transition +++++++++++++++++++++++++++++++


hh<-df.belaya
trans <- attr(hh, "trans")
K <- nrow(trans)
if (!is.null(dimnames(trans))) 
  states <- dimnames(trans)[[1]]
from <- factor(hh$from,labels = states)
to <- factor(hh$to, labels = states)

tbl <- table(hh$from[hh$Status == 1], hh$to[hh$Status == 
                                              1], dnn = c("from", "to"))
counts <- tbl
tbl <- table(hh$from, hh$to, dnn = c("from", "to"))
total <- apply(tbl, 1, max)
noevent <- total - apply(counts, 1, sum)
counts <- cbind(counts, noevent, total)
dn <- dimnames(counts)
dn[[2]][(K + 1):(K + 2)] <- c("no event", "total entering")
names(dn) <- c("from", "to")
dimnames(counts) <- dn
class(counts) <- "table"
freqs <- (counts/total)[, -(K + 2)]
class(freqs) <- "table"
list(Frequencies = counts, Proportions = freqs)


#-------------------- Non-Parameteric.........
fit1=coxph(Surv(Tstart,Tstop,Status)~strata(trans), data=df.belaya,method = "breslow")
summary(fit1)

msf0 <- msfit(object=fit1, trans=tmat1)

pt0 <- probtrans(msf0, predt = 0,method = "greenwood")

pt0 <- probtrans(msf0, predt = 10,method = "greenwood")

x11()
par(mfrow=c(1,2))
plot(pt0, from=1, type=c("single"),legend=c("d4T","AZT","TDF"),lwd=2,col=c(1:3),main="From d4T",xlab="Month on ART")
plot(pt0, from=2, type=c("single"),legend=c("d4T","AZT","TDF"),lwd=2,col=c(1:3),main="From AZT",xlab="Month on ART")
#plot(pt0, from=3, type=c("single"),legend=c(1:3),lwd=2,col=c(1:3),main="From 3",xlab="Month on ART")

#----------- merging the multistate data with covariates...............
hhdesc=read.csv("D:\\PHDfile\\PHDproject1\\Multistatemodeling\\Data\\FulldataMergedfordescriptive.csv")
head(hhdesc)

hhdescselected<-hhdesc[,c(1,11:19,23)]
head(hhdescselected)
names(hhdescselected)<-c("id","Age", "Sex", "Wei", "MarStat","EducLev", "Emp", "ClinStag", "FunStat", "EntryDate", "CD4" )

df.belaya.covariates=merge(df.belaya, hhdescselected, by = "id", sort = FALSE)
head(df.belaya.covariates)


df.belaya.covariates$Sex<-as.factor(df.belaya.covariates$Sex)
df.belaya.covariates$MarStatgrouped<-as.factor(df.belaya.covariates$MarStat>=2)
df.belaya.covariates$Agegrouped<-as.factor(df.belaya.covariates$Age>=41)
df.belaya.covariates$weightgrouped<-as.factor(df.belaya.covariates$Wei>=50)
df.belaya.covariates$EduLevelgrouped<-as.factor(df.belaya.covariates$EducLev>=2)
df.belaya.covariates$Empgrouped<-as.factor(df.belaya.covariates$Emp>=2)
df.belaya.covariates$ClinStaggrouped<-as.factor(df.belaya.covariates$ClinStag>=3)
df.belaya.covariates$FunStatgrouped<-as.factor(df.belaya.covariates$FunStat>=1)
df.belaya.covariates$CD4Bsegrouped<-  as.factor(df.belaya.covariates$CD4>=200)

attr(df.belaya.covariates, "trans") <- tmat1
class(df.belaya.covariates) <- c("msdata","data.frame")

### Model 1: Type = Transition  +  Covariates have identical effect ................. 

fit1=coxph(Surv(Tstart,Tstop,Status)~Sex+MarStatgrouped +Agegrouped +weightgrouped
           +EduLevelgrouped+Empgrouped+ClinStaggrouped+FunStatgrouped+CD4Bsegrouped+strata(trans), data=df.belaya.covariates,method = "breslow")
summary(fit1)

rr1 <- redrank(Surv(Tstart,Tstop,Status)~Sex+MarStatgrouped +Agegrouped +weightgrouped
               +EduLevelgrouped+Empgrouped+ClinStaggrouped+FunStatgrouped+CD4Bsegrouped+strata(trans), data=df.belaya.covariates, R = 1, print.level = 0)


### Model 2: Type = Transition  +  Covariates are type specific effect ................. 

expcovs<-expand.covs (df.belaya.covariates, covs=c("Sex", "MarStatgrouped","Agegrouped","weightgrouped","EduLevelgrouped","Empgrouped","ClinStaggrouped","CD4","FunStatgrouped"),  append = TRUE)
head(expcovs)


### Model 2: Type < Transition  +  Covariates are type specific effect ................. 

fitseparteSex=coxph(Surv(Tstart,Tstop,Status)~Sex1.1+Sex1.2+Sex1.3+Sex1.4+strata(trans), data=expcovs,method = "breslow")
summary(fitseparteSex)


fitseparteCD4B=coxph(Surv(Tstart,Tstop,Status)~CD4BsegroupedTRUE.1+CD4BsegroupedTRUE.2+CD4BsegroupedTRUE.3+CD4BsegroupedTRUE.4
                       +strata(trans), data=expcovs,method = "breslow")
summary(fitseparteCD4B)

fitseparteCD4B2=coxph(Surv(Tstart,Tstop,Status)~CD4.1+CD4.2+CD4.3+CD4.4
                     +strata(trans), data=expcovs,method = "breslow")
summary(fitseparteCD4B2)


fitseparteMarStat=coxph(Surv(Tstart,Tstop,Status)~MarStatgroupedTRUE.1+MarStatgroupedTRUE.2+MarStatgroupedTRUE.3+MarStatgroupedTRUE.4+strata(trans), data=expcovs,method = "breslow")
summary(fitseparteMarStat)


fitseparteAge=coxph(Surv(Tstart,Tstop,Status)~AgegroupedTRUE.1+AgegroupedTRUE.2+AgegroupedTRUE.3+AgegroupedTRUE.4+strata(trans), data=expcovs,method = "breslow")
summary(fitseparteAge)

fitseparteweight=coxph(Surv(Tstart,Tstop,Status)~weightgroupedTRUE.1+weightgroupedTRUE.2+weightgroupedTRUE.3+weightgroupedTRUE.4+strata(trans), data=expcovs,method = "breslow")
summary(fitseparteweight)



fitseparteEduLevel=coxph(Surv(Tstart,Tstop,Status)~EduLevelgroupedTRUE.1+EduLevelgroupedTRUE.2+EduLevelgroupedTRUE.3+EduLevelgroupedTRUE.4+strata(trans), data=expcovs,method = "breslow")
summary(fitseparteEduLevel)


fitseparteEmp=coxph(Surv(Tstart,Tstop,Status)~EmpgroupedTRUE.1+EmpgroupedTRUE.2+EmpgroupedTRUE.3+EmpgroupedTRUE.4+strata(trans), data=expcovs,method = "breslow")
summary(fitseparteEmp)

fitseparteClinStag=coxph(Surv(Tstart,Tstop,Status)~ClinStaggroupedTRUE.1+ClinStaggroupedTRUE.2+ClinStaggroupedTRUE.3+ClinStaggroupedTRUE.4+strata(trans), data=expcovs,method = "breslow")
summary(fitseparteClinStag)


fitseparteFunStat=coxph(Surv(Tstart,Tstop,Status)~FunStatgroupedTRUE.1+FunStatgroupedTRUE.2+FunStatgroupedTRUE.3+FunStatgroupedTRUE.4+strata(trans), data=expcovs,method = "breslow")
summary(fitseparteFunStat)

cbind(round(summary(fitseparteSex)$coefficients[,c(2,5)],digits=3), round(summary(fitseparteCD4B2)$coefficients[,c(2,5)],digits=3), round(summary(fitseparteMarStat)$coefficients[,c(2,5)],digits=3),
      round(summary(fitseparteAge)$coefficients[,c(2,5)],digits=3),round(summary(fitseparteweight)$coefficients[,c(2,5)],digits=3), round(summary(fitseparteEduLevel)$coefficients[,c(2,5)],digits=3),
      round(summary(fitseparteEmp)$coefficients[,c(2,5)],digits=3),round(summary(fitseparteClinStag)$coefficients[,c(2,5)],digits=3), round(summary(fitseparteFunStat)$coefficients[,c(2,5)],digits=3))



fitseparteall=coxph(Surv(Tstart,Tstop,Status)~Sex1.1+Sex1.2+Sex1.3+Sex1.4+CD4BsegroupedTRUE.1+CD4BsegroupedTRUE.2+CD4BsegroupedTRUE.3+CD4BsegroupedTRUE.4+MarStatgroupedTRUE.1+MarStatgroupedTRUE.2+MarStatgroupedTRUE.3+MarStatgroupedTRUE.4+AgegroupedTRUE.1+AgegroupedTRUE.2+AgegroupedTRUE.3+AgegroupedTRUE.4+EduLevelgroupedTRUE.1+EduLevelgroupedTRUE.2+EduLevelgroupedTRUE.3+EduLevelgroupedTRUE.4+EmpgroupedTRUE.1+EmpgroupedTRUE.2+EmpgroupedTRUE.3+EmpgroupedTRUE.4+
                                              ClinStaggroupedTRUE.1+ClinStaggroupedTRUE.2+ClinStaggroupedTRUE.3+ClinStaggroupedTRUE.4+FunStatgroupedTRUE.1+FunStatgroupedTRUE.2+FunStatgroupedTRUE.3+FunStatgroupedTRUE.4+strata(trans), data=expcovs,method = "breslow")
summary(fitseparteall)
