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

#Generate pairwise comparisions. P-value is set to 0.05, default in DESeq is 0.1
results_all<-results(dds_all, contrast=c("infection", "inf", "non-inf"), alpha=0.05) 	#All infected or non-infected samples
results_am<-results(dds_sep, contrast=c("condition", "inf-am", "non-inf-am"), alpha=0.05) 	#Infected/non-infected in AM
results_pm<-results(dds_sep, contrast=c("condition", "inf-pm", "non-inf-pm"), alpha=0.05)	#Infected/non-infected in PM

summary(results_all)
summary(results_am)
summary(results_pm)

#MA Plot comparing mean expression to log-fold change (quality check)
plotMA(results_all, main="Infected versus Non-Infected--All")
plotMA(results_am, main="Infected versus Non-Infected--AM")
plotMA(results_pm, main="Infected versus Non-Infected--PM")

##Can also get plot counts for individual genes--once annotation occurs, this may be useful. Check documentation

#Export results to text files
write.csv(as.data.frame(results_all), file="results_all_deseq2.csv")
write.csv(as.data.frame(results_am), file="results_am_deseq2.csv")
write.csv(as.data.frame(results_pm), file="results_pm_deseq2.csv")

#######Working on visualization later...
