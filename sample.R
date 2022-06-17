#all microarray files are supposed to be placed in the current directory
#files whose names start with "SG" -> mRNA expression
#files whose names start with "SH" -> miRNA expression
#a file "data sheet for microarray.csv" is also supposed to be located in the current directory

#loading mRNA expression
files <-list.files(pattern="SG")

x_all <-NULL
for (i in c(1:length(files)))
{
    x <- read.csv(files[i],sep="\t",skip=9)
    x_all <- cbind(x_all,x$gProcessedSignal)
}
x_all <- data.frame(x$SystematicName,x_all)
x_all <- x_all[x$ControlType==0,]

#class sample labels
class <- y[match(files,y[,6]),5]


#loading miRNAexpression
files <-list.files(pattern="SH")
y <- read.csv("data sheet for microarray.csv",sep="\t")
class_m <- as.character(y[match(substring(files,1,13),y[,4]),5])
x_miRNA_all <-NULL
for (i in c(1:length(files)))
{
    cat(i," ")
    x <- read.csv(files[i],sep="\t",skip=7)
    x_miRNA_all <- cbind(x_miRNA_all,x$X635nm)
}
x_miRNA_all <- data.frame(x$G_Name,x_miRNA_all)
x_miRNA_all <- x_miRNA_all[grep("hsa",x_miRNA_all[,1]),]


#PCA to mRNMA
pca <- prcomp(scale(x_all[,-1]))

#PCA to miRNA and miRNA selection
pcam<-prcomp(scale(x_miRNA_all[,-1]))

P <- pchisq(scale(pcam$x[,2])^2,1,lower.tail=F)
index <- p.adjust(P,"BH")<0.01
table(index)
#index
#FALSE  TRUE 
#2541    24 
miRNA <-  x_miRNA_all[index,1]
miRNA
#[1] hsa-miR-4745-5p hsa-miR-5100    hsa-miR-3648    hsa-miR-204-3p 
#[5] hsa-miR-7975    hsa-miR-6729-5p hsa-miR-762     hsa-miR-3665   
#[9] hsa-miR-5787    hsa-miR-1260b   hsa-miR-4787-5p hsa-miR-6090   
#[13] hsa-miR-7704    hsa-miR-4466    hsa-miR-4286    hsa-miR-4488   
#[17] hsa-miR-7977    hsa-miR-3960    hsa-miR-4454    hsa-miR-6869-5p
#[21] hsa-miR-638     hsa-miR-4497    hsa-miR-4508    hsa-miR-6089 

#correlation of PC2 between mRNA and miRNA
PC2 <- unlist(lapply(split(pcam$rotation[,2],class_m),mean))
cor.test(pca$rotation[,2],PC2[match(class,names(PC2))])

#Pearsons product-moment correlation

#data:  pca$rotation[, 2] and PC2[match(class, names(PC2))]
#t = -3.0434, df = 6, p-value = 0.0227
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#    -0.9578669 -0.1648261
#sample estimates:
#    cor 
# -0.7790174 

#figures
pdf(file="pca.pdf")
par(mfrow=c(1,2))
plot(pca$rotation[,2],col=1:2,type="h",xaxt="n",xlab="",ylab="PC2 loadings",main="mRNA")
points(pca$rotation[,2],col=1:2)
abline(0,0,col=2,lty=4)
mtext(class,at=1:8,las=3,side=1)
plot(-PC2[match(class,names(PC2))],col=1:2,type="h",xaxt="n",xlab="",ylab="- PC2 loadings",main="miRNA")
points(-PC2[match(class,names(PC2))],col=1:2)
abline(0,0,col=2,lty=4)
mtext(class,at=1:8,las=3,side=1)
par(mfrow=c(1,1))
COR <- cor.test(pca$rotation[,2],PC2[match(class,names(PC2))])
plot(pca$rotation[,2],PC2[match(class,names(PC2))],xlab="PC2 loadings: mRNA",ylab="PC2 loadiings: miRNA",cex=2,xlim=c(-0.8,0.55),ylim=c(-0.12,0.01),main=paste("P=",sprintf("%3.2e",COR$p.value),"\n COR=",sprintf("%3.2e",COR$estimate)))
text(pca$rotation[,2],PC2[match(class,names(PC2))],labels=class,pos=3)
plot(pca$x,col=index2+1,main="mRNA PC score")
plot(pcam$x,col=index+1,main="miRNA PC score")
dev.off()

pdf(file="boxplot.pdf")
par(mai=c(1.5,0.5,0.5,0.5))
par(mfrow=c(2,2))
for (i in c(1:sum(index)))
{
    boxplot(scale(x_miRNA_all[,-1])[index,][i,]~class_m,las=3,main=x_miRNA_all[,1][index][i],xlab="")
}
par(mfrow=c(1,1))
dev.off()