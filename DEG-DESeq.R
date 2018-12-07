#!/usr/bin/env Rscript

source("http://www.Bioconductor.org/biocLite.R")
library( "DESeq2" )
library("RColorBrewer")
library("ggplot2") #Best plots
library("ggrepel") #Avoid overlapping labels
library("ggpubr")

#PATH TO FILES
inputcondition="/Volumes/AMINABE/tt/wellsonmRNA-lt.csv"
inputMcount="/Volumes/AMINABE/tt/results/mRNA/MergedCountsmRNA-somegroups.csv"

outputfile3="/PATH/TO/mRNAcountsResults.txt"
outputfile2="/PATH/TO/MergedCountsmRNAgenename.csv"
outputfile1="/PATH/TO/CountsmRNAgenename.csv"


###READ FILES
#Lenth file should contain the length of each gene. in this case it is downloaded from ensembl
Length=read.table("/Volumes/AMINABE/tt/gene/W-addapter_gene_name_length.txt", header=T,sep="\t",row.names = 1)

sri = read.csv(inputcondition, stringsAsFactors=FALSE, sep=";")
MergedCount=read.table(inputMcount, header=T,sep=";",row.names = 1)
head(MergedCount)
head(sri)

#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
dds = DESeqDataSetFromMatrix(countData = MergedCount , colData = sri , design =  ~groups)
dds = estimateSizeFactors(dds)
sizeFactors(dds)
idx <- rowSums(counts(dds) >= 5 ) >= 12 #, normalized=TRUE
#idx <- rowSums( counts(dds, normalized=TRUE) >= 5) >=12
dds <- dds[idx,]
#_________________________________________________
##Differential expression analysis
#_________________________________________________
dds <- DESeq(dds)

#_________________________________________________
#Results
#_________________________________________________
res <- results(dds,  contrast=c("groups","B","A"))
mcols(res, use.names=TRUE)
respadjfc <- res[ order(res$padj, decreasing = FALSE), ]
up=sum(res$padj < 0.05 & res$log2FoldChange >= 0.75 , na.rm=TRUE)
down=sum(res$padj < 0.05 & res$log2FoldChange <= -0.75, na.rm=TRUE)
up
down
respadjfc <- res[ order(res$padj, decreasing = FALSE), ]
sum(res$padj < 0.05, na.rm=TRUE)
head(respadjfc)

################################################################################
################################################################################
#write results and add the name of the genes
################################################################################
################################################################################

respadjfc <- res[ order(res$padj, decreasing = FALSE), ]
data2<-counts(dds, normalized=TRUE)

write.table(respadjfc,file="/PATH/TO/mRNA-results.txt",sep="\t",quote=F)
write.table(counts(dds, normalized=TRUE), file="/PATH/TO/mRNA-normalized.txt", sep="\t")

countsfile=respadjfc

#//////////////////////////// some extrat files only ////////////////////////////////////////


length=data.frame(rownames(Length),Length$name, Length$length)

#table of results:
idx=match(length$rownames.Length,rownames(countsfile))
Lengthmatch = as.data.frame(length[!is.na(idx),])
tedgeRmatch = countsfile[idx[!is.na(idx)],]
tt=cbind(tedgeRmatch, Lengthmatch)
#write.table(tt[order(tt$padj),],file=outputfile2,sep=";",quote=F)
#table of counts
idx=match(length$rownames.Length,rownames(data2))
Lengthmatch = as.data.frame(length[!is.na(idx),])
tedgeRmatch = data2[idx[!is.na(idx)],]
tt2=cbind(tedgeRmatch, Lengthmatch)
write.table(tt2[order(tt2$padj),],file=outputfile1,sep=";",quote=F) #needed later to plot the volcano
#Table counts and results
idx=match(rownames(tt),rownames(tt2))
Lengthmatch = as.data.frame(tt[!is.na(idx),])
tedgeRmatch = tt2[idx[!is.na(idx)],]
tt3=cbind(tedgeRmatch, Lengthmatch)
write.table(tt3[order(tt3$padj),],file=outputfile3,sep=";",quote=F)

################################################################################
################################################################################
#plot
################################################################################
################################################################################
#_________________________________________________
#ploting PCA
#_________________________________________________
#rld = rlogTransformation(dds)
library(ggplot2)
rld <- rlogTransformation(dds, blind=TRUE)
png("/PATH/TO/DAmRNAPCA.png", width = 7*300, height = 5*300, res = 400, pointsize = 12)
#plotPCA(rld, intgroup=c('groups', 'samples"),ntop = 6000, returnData = F)
p <- plotPCA(rld,intgroup=c('groups'),ntop = 6000, returnData = F)+theme_light()
#p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")+theme_light()
print(p)
dev.off()
#_________________________________________________
#Ploting dispertionALL
#_________________________________________________
png("/PATH/TO/mRNADispersions.png", width = 7*300, height = 5*300, res = 400, pointsize = 12)
plotDispEsts(dds)
dev.off()
#_________________________________________________
#plot pval
#_________________________________________________

png("/PATH/TO/mRNA_pval.png", width = 7*300, height = 5*300, res = 400, pointsize = 12)
hist(res$pvalue,breaks=100,col="grey", xlab="p-value",main="p-value distribution")
dev.off()
#_________________________________________________
# MA plots
#_________________________________________________

name=read.table("/PATH/TO/MergedCountsmRNAgenename.csv", header=T,sep=";",row.names = 1)

#another type of maplot
png("/PATH/TO/mRNA_MAggplot.png", width = 7*300, height = 5*300, res = 400, pointsize = 12)
ggmaplot(name, main = expression("B. vs A."), #for mRNA
fdr = 0.05, fc = 0.75, size = 1,
palette = c("#B31B21", "#1465AC", "darkgray"),
genenames = as.vector(name[,"Length.name"]),#for mRNA
legend = "top", top = 30,
font.label = c("bold", 5),
font.legend = "bold",
font.main = "bold",
ggtheme = ggplot2::theme_light())
dev.off()

#_________________________________________________
#VOlcano plot
#data=head(mydata, 10) :> the "10" should be modified depending on the number of your genes
#_________________________________________________

png("PATH/TO/mRNA_volcanoggplot.png", width = 7*400, height = 5*300, res = 400, pointsize = 12)
mydata<-name%>%mutate(threshold = ifelse(log2FoldChange >= 0.75 & padj<0.05 ,"UP", ifelse(log2FoldChange<=-0.75 & padj<0.05, "Down", "NA")))
ggplot(mydata, aes(x=log2FoldChange, y=-log10(padj)),fdr = 0.05, fc = 0.75, size = 100,) + xlim(-5, 5)+ylim(0,10)+
geom_point(aes(colour = threshold), size=1.5) +theme_light()+
theme(axis.text.x=element_text(size=rel(1.5)),axis.text.y=element_text(size=rel(1.5)), text = element_text(size=14))+
scale_colour_manual(values = c("UP"= "#B31B21", "Down"="#1465AC",  "NA"= "darkgray"))+
geom_text_repel(data=head(mydata, 10), aes(label=Length.name) ,size = 4.4)+ #adding text for the top 20 genes
geom_vline(xintercept=0, linetype="dashed",color = "black", size=.5)

dev.off()

#END
