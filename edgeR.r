library("edgeR")
#Change to working directory with R Studio interface

#Load data
biol_reps_data<-read.table("kallisto_master.tsv",stringsAsFactors=FALSE)

#Change column and row names to samples and genes, so that data frame cells are just count values
colnames(biol_reps_data)<-0:16
rownames(biol_reps_data)<-biol_reps_data[,1]
biol_reps_data<-biol_reps_data[,-1]
biol_reps_data<-biol_reps_data[-1,]

#Convert to matrix and convert characters to numeric type
biol_reps_numeric<-as.matrix(biol_reps_data)
class(biol_reps_numeric)<-"numeric"

#Load data matrix into DGEList format to begin edgeR pipepline
#Groups are as follows: 4 morning non-infected, 4 afternoon non-infected, 4 morning infected, 4 afternoon infected
group<-c(rep("non-inf-am",4),rep("non-inf-pm",4),rep("inf-am",4),rep("inf-pm",4))
y<-DGEList(biol_reps_numeric,group=group)

#Filter variants: Recommended defaults here are 10 counts in smallest library, translated to counts/million (here ~11 million = 1 cpm), 
#and in at least 4 libraries (or expressed in all members of one group).
keep<-rowSums(cpm(y)>0.5)>=4
y<-y[keep, , keep.lib.sizes=FALSE] #Recalculate library sizes afterwards

#Calculate normalization factors
y<-calcNormFactors(y)

#Construct multi-dimensional graph for clustering of samples
colors<-rep(c(rep(c("blue"),4),rep(c("red"),4)),2) #Color by infected/non-infected
points<-c(rep(c(16),8),rep(c(17),8)) #Point shape by morning/evening
plotMDS(y, col=colors, pch=points)
abline(h=0, col="green", lty=2, lwd=2) #add line for sperating AM and PM
#Legend will need some work
legend("center", legend=c("Non-infected morning","Infected morning","Non-infected evening","Infected evening"),col=c("blue","red","blue","red"),pch=c(16,16,17,17))
title(main="Multi-Dimensional analysis of Turfgrass transcriptome")

