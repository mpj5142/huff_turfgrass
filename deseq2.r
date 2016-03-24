library("DESeq2")
#Change to working directory with R Studio interface

#Load data
biol_reps_data<-read.table("kallisto_master.tsv",stringsAsFactors=FALSE)

#Change column and row names to samples and genes, so that data frame cells are just count values
colnames(biol_reps_data)<-0:16
rownames(biol_reps_data)<-biol_reps_data[,1]
biol_reps_data<-biol_reps_data[,-1]
biol_reps_data<-biol_reps_data[-1,]

#Convert to matrix and convert characters to numeric type/round to integers
biol_reps_numeric<-as.matrix(biol_reps_data)
class(biol_reps_numeric)<-"numeric"
biol_reps_numeric<-round(biol_reps_numeric)

#Design data frame with condition data for the experiment
#Two separate runs will be one: One (dds_sep) that separates AM and PM samples; one (dds_all) with combined AM/PM samples
#These cannot be done in the same run due to limitations with how DESeq2 factors the desgin matrix; it cannot differentiate between the conditions
condition<-c(rep("non-inf-am",4),rep("inf-am",4),rep("non-inf-pm",4),rep("inf-pm",4)) 
infection<-c(rep("non-inf",4),rep("inf",4),rep("non-inf",4),rep("inf",4))
#time<-c(rep("am",8),rep("pm",8)) #Include this for later analysis if necessary
sample<-c(1:16)
sample.data.sep<-data.frame(sample,condition)
sample.data.all<-data.frame(sample, infection)

#Load data matrix into DESeq format to begin pipepline. Design will account for infection and time
dds_sep<-DESeqDataSetFromMatrix(countData=biol_reps_numeric, colData=sample.data.sep, design= ~ condition)
dds_all<-DESeqDataSetFromMatrix(countData=biol_reps_numeric, colData=sample.data.all, design= ~ infection)

#Pre-filter to remove rows with <2 total counts; additional filtering will be applied later
dds_sep<-dds_sep[rowSums(counts(dds_sep))>1,]
dds_all<-dds_all[rowSums(counts(dds_all))>1,]

#Perform differential expression analysis
#These steps take a while...
dds_sep<-DESeq(dds_sep)
dds_all<-DESeq(dds_all)

#Generate pairwise comparisions
results_all<-results(dds_all, contrast=c("infection", "inf", "non-inf")) 	#All infected or non-infected samples
results_am<-results(dds_sep, contrast=c("condition", "inf-am", "non-inf-am")) 	#Infected/non-infected in AM
results_pm<-results(dds_sep, contrast=c("condition", "inf-pm", "non-inf-pm"))	#Infected/non-infected in PM


#######Working on visualization/output later...


#Construct multi-dimensional graph for clustering of samples
colors<-rep(c(rep(c("blue"),4),rep(c("red"),4)),2) #Color by infected/non-infected
points<-c(rep(c(16),8),rep(c(17),8)) #Point shape by morning/evening
plotMDS(y, col=colors, pch=points)
abline(h=0, col="green", lty=2, lwd=2) #add line for sperating AM and PM
#Legend will need some work
legend("center", legend=c("Non-infected morning","Infected morning","Non-infected evening","Infected evening"),col=c("blue","red","blue","red"),pch=c(16,16,17,17))
title(main="Multi-Dimensional analysis of Turfgrass transcriptome")

#Write as output all (unfiltered) differentialy expressed genes, sorted by FDR
write.csv(topTags(diffexp_all, n=nrow(diffexp_all$table))$table, file="diffexp_all.csv")
write.csv(topTags(diffexp_am, n=nrow(diffexp_am$table))$table, file="diffexp_am.csv")
write.csv(topTags(diffexp_pm, n=nrow(diffexp_pm$table))$table, file="diffexp_pm.csv")
write.csv(topTags(diffexp_inf, n=nrow(diffexp_inf$table))$table, file="diffexp_inf.csv")
write.csv(topTags(diffexp_noninf, n=nrow(diffexp_noninf$table))$table, file="diffexp_noninf.csv")
write.csv(topTags(diffexp_time, n=nrow(diffexp_time$table))$table, file="diffexp_time.csv")
