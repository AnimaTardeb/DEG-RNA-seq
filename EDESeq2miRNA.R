#!/usr/bin/env Rscript

source("http://www.Bioconductor.org/biocLite.R")
library( "DESeq2" )
library("RColorBrewer")
library("ggplot2") #Best plots
library("ggrepel") #Avoid overlapping labels
library("ggpubr")



#inputcondition="/Volumes/AMINABE/tt/wellsonvrais2.csv"
#inputMcount="/Volumes/AMINABE/tt/results/miRNA/MergedCountsmiRNA.csv"
#some groups
inputcondition="/Volumes/AMINABE/tt/wellsonvraismir-somegroups.csv"
inputMcount="/Volumes/AMINABE/tt/results/miRNA/MergedCounts3-somegroups.csv"

#inputcondition="/Volumes/AMINABE/tt/wellsonvraismirAC.csv"
#inputMcount="/Volumes/AMINABE/tt/results/miRNA/MergedCountsAC.csv"

#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
###read files
sri = read.csv(inputcondition, stringsAsFactors=FALSE, sep=";")
MergedCount=read.table(inputMcount, header=T,sep=";",row.names = 1)
head(MergedCount)
head(sri)

#////////////////////////////////////////////////////////////////////
#BC groups
sri=data.frame(rbind(sri[4,],sri[5,],sri[6,],sri[7,],sri[8,],sri[9,]))
MergedCount2=data.frame(MergedCount[,"X5_TriAlimaturemicrna"],MergedCount[,"X6_TriAlimaturemicrna"],MergedCount[,"X8_TriAlimaturemicrna"],MergedCount[,"X9_TriAlimaturemicrna"],MergedCount[,"X10_TriAlimaturemicrna"],MergedCount[,"X12_TriAlimaturemicrna"])
colnames(MergedCount2)=c("X5_TriAlimaturemicrna","X6_TriAlimaturemicrna","X8_TriAlimaturemicrna","X9_TriAlimaturemicrna","X10_TriAlimaturemicrna","X12_TriAlimaturemicrna")
rownames(MergedCount2)=rownames(MergedCount)
MergedCount=MergedCount2
#////////////////////////////////////////////////////////////////////
#BD
sri=data.frame(rbind(sri[4,],sri[5,],sri[6,],sri[10,],sri[11,],sri[12,]))
MergedCount2=data.frame(MergedCount[,"X5_TriAlimaturemicrna"],MergedCount[,"X6_TriAlimaturemicrna"],MergedCount[,"X8_TriAlimaturemicrna"],MergedCount[,"X13_TriAlimaturemicrna"],MergedCount[,"X14_TriAlimaturemicrna"],MergedCount[,"X15_TriAlimaturemicrna"])
colnames(MergedCount2)=c("X5_TriAlimaturemicrna","X6_TriAlimaturemicrna","X8_TriAlimaturemicrna","X13_TriAlimaturemicrna","X14_TriAlimaturemicrna","X15_TriAlimaturemicrna")
rownames(MergedCount2)=rownames(MergedCount)
MergedCount=MergedCount2
#////////////////////////////////////////////////////////////////////
#CD
sri=data.frame(rbind(sri[7,],sri[8,],sri[9,],sri[10,],sri[11,],sri[12,]))
MergedCount2=data.frame(MergedCount[,"X9_TriAlimaturemicrna"],MergedCount[,"X10_TriAlimaturemicrna"],MergedCount[,"X12_TriAlimaturemicrna"], MergedCount[,"X13_TriAlimaturemicrna"],MergedCount[,"X14_TriAlimaturemicrna"],MergedCount[,"X15_TriAlimaturemicrna"])
colnames(MergedCount2)=c("X9_TriAlimaturemicrna","X10_TriAlimaturemicrna","X12_TriAlimaturemicrna","X13_TriAlimaturemicrna","X14_TriAlimaturemicrna","X15_TriAlimaturemicrna")
rownames(MergedCount2)=rownames(MergedCount)
MergedCount=MergedCount2
#////////////////////////////////////////////////////////////////////
#BC groups
sri=data.frame(rbind(sri[4,],sri[5,],sri[6,],sri[7,],sri[8,],sri[9,]))
MergedCount2=data.frame(MergedCount[,"X5_mapped"],MergedCount[,"X6_mapped"],MergedCount[,"X8_mapped"],MergedCount[,"X9_mapped"],MergedCount[,"X10_mapped"],MergedCount[,"X12_mapped"])
colnames(MergedCount2)=c("X5_mapped","X6_mapped","X8_mapped","X9_mapped","X10_mapped","X12_mapped")
rownames(MergedCount2)=rownames(MergedCount)
MergedCount=MergedCount2
#////////////////////////////////////////////////////////////////////
#BD
sri=data.frame(rbind(sri[4,],sri[5,],sri[6,],sri[10,],sri[11,],sri[12,]))
MergedCount2=data.frame(MergedCount[,"X5_mapped"],MergedCount[,"X6_mapped"],MergedCount[,"X8_mapped"],MergedCount[,"X13_mapped"],MergedCount[,"X14_mapped"],MergedCount[,"X15_mapped"])
colnames(MergedCount2)=c("X5_mapped","X6_mapped","X8_mapped","X13_mapped","X14_mapped","X15_mapped")
rownames(MergedCount2)=rownames(MergedCount)
MergedCount=MergedCount2

#////////////////////////////////////////////////////////////////////
#CD
sri=data.frame(rbind(sri[7,],sri[8,],sri[9,],sri[10,],sri[11,],sri[12,]))
MergedCount2=data.frame(MergedCount[,"X9_mapped"],MergedCount[,"X10_mapped"],MergedCount[,"X12_mapped"], MergedCount[,"X13_mapped"],MergedCount[,"X14_mapped"],MergedCount[,"X15_mapped"])
colnames(MergedCount2)=c("X9_mapped","X10_mapped","X12_mapped","X13_mapped","X14_mapped","X15_mapped")
rownames(MergedCount2)=rownames(MergedCount)
MergedCount=MergedCount2
#////////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////////
##############################################################################################
##############################################################################################

dds = DESeqDataSetFromMatrix(countData = MergedCount , colData = sri , design =  ~groups)
#sizeFactors(dds)
dds = estimateSizeFactors(dds)

idx <- rowSums(counts(dds) >= 5 ) >= 12 #, normalized=TRUE
#idx <- rowSums( counts(dds, normalized=TRUE) >= 5) >=12
dds <- dds[idx,]
#_________________________________________________
##Differential expression analysis
#_________________________________________________
dds <- DESeq(dds)

#ploting PCA
#_________________________________________________
#rld = rlogTransformation(dds)
library(ggplot2)
rld <- rlogTransformation(dds, blind=TRUE)
png("/Volumes/AMINABE/tt/results/miRNA/Figures/DA/DAmRNAPCA.png", width = 7*300, height = 5*300, res = 400, pointsize = 12)
#plotPCA(rld, intgroup=c('groups', 'samples"),ntop = 6000, returnData = F)
p <- plotPCA(rld,intgroup=c('groups'),ntop = 6000, returnData = F)+theme_light()
#p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")+theme_light()
print(p)
dev.off()

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
write.table(respadjfc,file="/Volumes/AMINABE/tt/results/miRNA/Figures/CD/CDmiRNA-results.txt",sep="\t",quote=F)
write.table(counts(dds, normalized=TRUE), file="/Volumes/AMINABE/tt/results/miRNA/Figures/CD/CDmiRNA-normalized.txt", sep="\t")
countsfile=respadjfc
Length=read.table("/Volumes/AMINABE/tt/gene/W-addapter_gene_name_length.txt", header=T,sep="\t",row.names = 1)

#////////////////////////////mRNA only ////////////////////////////////////////
#respadjfc <- res[ order(res$padj, decreasing = FALSE), ]
#data2<-counts(dds, normalized=TRUE)
outputfile3="/Volumes/AMINABE/tt/results/miRNA/Figures/CD/CDmiRNAcountsResults.txt"
outputfile2="/Volumes/AMINABE/tt/results/miRNA/Figures/CD/CDMergedCountsmRNAgenename.csv"
outputfile1="/Volumes/AMINABE/tt/results/miRNA/Figures/CD/CDCountsmRNAgenename.csv"

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
#sata
#write.table(tt2[order(tt2$padj),],file=outputfile1,sep=";",quote=F)

#miRNA
#
idx=match(rownames(countsfile),rownames(data2))
Lengthmatch = as.data.frame(countsfile[!is.na(idx),])
tedgeRmatch = data2[idx[!is.na(idx)],]
tt3=cbind(tedgeRmatch, Lengthmatch)
#sata
write.table(tt3[order(tt3$padj),],file=outputfile3,sep=";",quote=F)

################################################################################
################################################################################
#plot
################################################################################
################################################################################

#Ploting dispertionALL
#_________________________________________________
png("/Volumes/AMINABE/tt/results/miRNA/Figures/CD/CDmiRNADispersions.png", width = 7*300, height = 5*300, res = 400, pointsize = 12)
plotDispEsts(dds)
dev.off()
#_________________________________________________
#plot pval
png("/Volumes/AMINABE/tt/results/miRNA/Figures/CD/CDmiRNA_pval.png", width = 7*300, height = 5*300, res = 400, pointsize = 12)
hist(res$pvalue,breaks=100,col="grey", xlab="p-value",main="p-value distribution")
dev.off()

########################################
#name=read.table("/Volumes/AMINABE/tt/results/miRNA/Figures/CD/CDMergedCountsmRNAgenename.csv", header=T,sep=";",row.names = 1)

#miRNA
name=read.table("/Volumes/AMINABE/tt/results/miRNA/Figures/CD/CDmiRNA-results.txt", header=T,sep="\t",row.names = 1)
#another type of maplot
png("/Volumes/AMINABE/tt/results/miRNA/Figures/CD/CDmiRNA_MAggplot.png", width = 7*300, height = 5*300, res = 400, pointsize = 12)
ggmaplot(res, main = expression("C. vs D."),
fdr = 0.05, fc = 0.75, size = 1,
palette = c("#B31B21", "#1465AC", "darkgray"),
genenames = as.vector(rownames(res)),
legend = "top", top = 4,
font.label = c("bold", 5),
font.legend = "bold",
font.main = "bold",
ggtheme = ggplot2::theme_light())
dev.off()
########################################
########################################

volcanoplot <- function (name, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topright", labelsig=TRUE, textcx=1, ...) {
    with(name, plot(log2FoldChange, -log10(padj), pch=20, main=main, ...))
    with(subset(name, padj<sigthresh ), points(log2FoldChange, -log10(padj), pch=20, col="red", ...))
    #with(subset(name, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(padj), pch=20, col="orange", ...))
    with(subset(name, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(padj), pch=20, col="green", ...))
    if (labelsig) {
        require(calibrate)
        #with(subset(name, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(padj), labs=Length.name, cex=textcx, ...))
        with(subset(name, padj<sigthresh), textxy(log2FoldChange, -log10(padj), labs=Length.name, cex=.75, pos=2, ...))
    }
    
    legend(legendpos, xjust=1, yjust=1, legend=c(paste("padj<",sigthresh,sep=" "), paste("padj<",sigthresh,sep="  ","|log2FC|>",lfcthresh)), pch=20, col=c("red","green"))
    abline(v=c(0), col="blue") ## plot the X axis
}
png("/Volumes/AMINABE/tt/results/miRNA/Figures/CD/CDVolcanoplot_mRNA.png",  width = 10*300, height = 7*300, res = 400, pointsize = 8)
volcanoplot(name, lfcthresh=0.75, sigthresh=0.05)#, xlim=c(-4, 4))

dev.off()

########################################
########################################
name=read.table("/Volumes/AMINABE/tt/results/miRNA/Figures/NEWNOfilter/BA/BAmiRNA-results.txt", header=T,sep="\t",row.names = 1)

png("/Volumes/AMINABE/tt/results/miRNA/Figures/NEWNOfilter/BA/BAmiRNA_volcanoggplot.png", width = 7*400, height = 5*300, res = 400, pointsize = 12)
mydata<-name%>%mutate(threshold = ifelse(log2FoldChange >= 0.75 & padj<0.05 ,"UP", ifelse(log2FoldChange<=-0.75 & padj<0.05, "Down", "NA")))
mydata<-mydata%>%mutate(Length.name = rownames(name))
ggplot(mydata, aes(x=log2FoldChange, y=-log10(padj)),fdr = 0.05, fc = 0.75, size = 100,) + xlim(-5, 5)+ylim(0,6)+
geom_point(aes(colour = threshold), size=1.5) +theme_light()+
theme(axis.text.x=element_text(size=rel(1.5)),axis.text.y=element_text(size=rel(1.5)), text = element_text(size=14))+
scale_colour_manual(values = c("UP"= "#B31B21", "Down"="#1465AC",  "NA"= "darkgray"))+
geom_text_repel(data=head(mydata, 3), aes(label=Length.name) ,size = 4.4)+ #adding text for the top 20 genes
geom_vline(xintercept=0, linetype="dashed",color = "black", size=.5)

dev.off()


#END
####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
/Volumes/AMINABE/tt/results/miRNA/Figures/NEWNOfilter/BA/BAmiRNA-results.txt

#another type of maplot
png("/Volumes/AMINABE/tt/results/mRNA/Figures/BA/miRNA_MAggplot.png", width = 7*300, height = 5*300, res = 400, pointsize = 12)
ggmaplot(name, main = expression("ctrl. vs treat."),
fdr = 0.05, fc = 0.75, size = 1,
palette = c("#B31B21", "#1465AC", "darkgray"),
#genenames = as.vector(rownames(res)),
genenames = as.vector(name[,"Length.name"]),

legend = "top", top = 10,
font.label = c("bold", 5),
font.legend = "bold",
font.main = "bold",
ggtheme = ggplot2::theme_light())
dev.off()

sum(res$padj < 0.05, na.rm=TRUE)





#Plotting counts for individual genes
#_________________________________________________
topGene <- rownames(res)[whiches(res$padj<0.05)]
png("/Volumes/AMINABE/tt/results/mRNA/Figures/BA/countforgenes4.png", width = 5*300, height = 5*300, res = 400, pointsize = 12)
plotCounts(dds, gene=c(topGene[1]), pch=16,intgroup=c("groups"))
dev.off()

#plot MA
#_________________________________________________
topGene <- c(rownames(res)[which(res$padj<0.05)])
png("/Volumes/AMINABE/tt/results/mRNA/Figures/BA/miRNAMA.png", width = 7*300, height = 5*300, res = 400, pointsize = 12)
plotMA(res, alpha = 0.05)
with(res[topGene, ], {
    points(baseMean, log2FoldChange, col="blue", cex=1.5, lwd=1.2)
    #text(baseMean, log2FoldChange, topGene, col="darkblue",cex=0.5)
})
dev.off()

#_________________________________________________

#Plot Volcanos
png("/Volumes/AMINABE/tt/results/mRNA/Figures/BA/Volcano.png", width = 5*300, height = 5*300, res = 400, pointsize = 9)
plot(res$log2FoldChange,-log10(res$padj), main = "ctrl. vs treat.",xlab="log2(FC)", ylab="-log10(padj)", col="black", pch=20)
with(subset(res, padj<.05 & abs(log2FoldChange) > .5), points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(res, padj<.05 & abs(log2FoldChange) > .75), points(log2FoldChange, -log10(padj), pch=20, col="blue"))

abline(v=c(0), col="blue")
dev.off()

#Plot Volcanos basemean>10
res_mean10 <- subset(res, baseMean > 10)

png("/Volumes/AMINABE/tt/results/mRNA/Figures/BA/Volcano_mean10.png", width = 5*300, height = 5*300, res = 400, pointsize = 8)
plot(res_mean10$log2FoldChange,-log10(res_mean10$padj), xlim=c(-5,5),main = "ctrl. vs treat.",xlab="log2(FC)", ylab="-log10(padj)", col="black", pch=16)
with(subset(res_mean10, padj<.05 & abs(log2FoldChange) > .75), points(log2FoldChange, -log10(padj), pch=20, col="red"))
abline(v=c(0), col="blue")
dev.off()


name=read.table("/Volumes/AMINABE/tt/results/mRNA/Figures/BA/MergedCountsmRNAgenename.csv", header=T,sep=";",row.names = 1)
volcanoData <- cbind(name$log2FoldChange, -log10(name$padj))
colnames(volcanoData) <- c("logFC", "negLogFDR")

#print ("==================head of volcanodata==================")
head(volcanoData)
plot(x=name$log2FoldChange, y=-log10(name$padj), xlab="log2(Fold-Change(FC))", ylab="-log10(adjusted P-value(FDR)", col=ifelse(name$padj <= 0.05,"red","black"), pch=16, main="Volcano plot")
legend("topright", c("-log10FDR<=0.05"), col=c("red"), pch=16,cex=0.75,box.lty=0)
abline(v=c(0), col="blue") ## plot the X axis

