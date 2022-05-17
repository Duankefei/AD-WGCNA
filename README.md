# AD-WGCNA
AD-WGCNA program
#diff
install.packages("ggplot2")
library(limma)
library(impute)
library(ggplot2)
setwd("C:\\Users\\APPLE\\Desktop\\WGCNA-AD1297\\3_diff")
mydata<-read.table("input.txt",header = T,sep="\t")
mydata<-as.matrix(mydata)
rownames(mydata)=mydata[,1]


geneExp=mydata[,2:ncol(mydata)]
matrixNa=list(rownames(geneExp),colnames(geneExp))
geneExp=matrix(as.numeric(as.matrix(geneExp)),nrow=nrow(geneExp),dimnames=matrixNa)


geneExp1=impute.knn(geneExp)
geneExp=geneExp1$data
geneExp=avereps(geneExp)
geneExp=log2(geneExp+1)
boxplot(geneExp,col="green")


geneExp<-normalizeBetweenArrays(geneExp,method = "scale")
boxplot(geneExp,col="blue")
write.table(geneExp,"geneExp.txt",sep = "\t",quote = F)

design <- model.matrix(~ 0+factor(c(rep("control",9),rep("AD",22))))

colnames(design) <- c("control", "AD")

fit <- lmFit(geneExp, design)
con_matrix<-contrast.matrix <- makeContrasts(AD-control, levels=design)
myfit <- contrasts.fit(fit, con_matrix)
myfit <- eBayes(myfit)
allgene<-topTable(myfit, coef=1, adjust="BH",number = "all")
write.table(allgene,"allgene.txt",sep = "\t",quote = F)
diffgene<-allgene[abs(allgene$logFC)>=0.5 & allgene$P.Value<0.05,]
write.table(diffgene,"diffgene.txt",sep = "\t",quote = F)
upgene<-allgene[allgene$logFC>=0.5 & allgene$P.Value<0.05,]
write.table(upgene,"upgene.txt",sep = "\t",quote = F)
downgene<-allgene[allgene$logFC<=-0.5 & allgene$P.Value<0.05,]
write.table(downgene,"downgene.txt",sep = "\t",quote = F)

#火山图
xmax<-max(allgene$logFC)
ymax<-max(-log10(allgene$P.Value))

allgene$sig = as.factor(ifelse(allgene$P.Value < 0.05 & abs(allgene$logFC) > 0.5, 
                               ifelse(allgene$logFC > 0.5,"Up", "Down"), "Not"))

a=ggplot(allgene,aes(logFC,-log10(P.Value))) 
b=a+ geom_point(aes(color =sig))+xlim(-xmax,xmax) + ylim(0,ymax)+ 
  labs(title="Volcano plot",x="log2FC", y="-log10(P.Value)")+theme(plot.title=element_text(hjust=0.5))+scale_color_manual(values =c("green","black", "red"))

pdf("Volcano.pdf")
b
dev.off()

#热图
library(pheatmap)
heatmapexp<-geneExp[row.names(diffgene),]
heatmapexp=heatmapexp[1:15,]
pdf("heatmap.pdf")
pheatmap(heatmapexp,
         scale = "row",#对行进行归一化
         clustering_method = "average",#选择聚类方法
         legend_labels = c("1.0","2.0","3.0","4.0","5.0"),#添加图例标签
         border=FALSE,#去掉边框线
         show_colnames=F,#是否显示列名，同理
         show_rownames=T,
         display_numbers = F,fontsize_row=8,fontsize_col=12,
         cluster_cols = F,cluster_rows = T,
         main = "GSE1297:normal vs AD heatmap",)
dev.off()


#WGCNA
#下载安装包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("GO.db")
BiocManager::install("preprocessCore")
BiocManager::install("impute")

install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "reshape", "fastcluster", "dynamicTreeCut", "survival"))
install.packages("WGCNA")
install.packages("scatterplot3d")


################################################
library(WGCNA)
setwd("C:\\Users\\APPLE\\Desktop\\WGCNA-AD1297\\4_WGCNA-R")  #设定自己的工作目录

inputdata1="input.txt" #改为自己的数据
data0=read.table(inputdata1,sep="\t",row.names=1,header=T,check.names=F,quote="!")
datSummary=rownames(data0)
datExpr = t(data0)
no.samples = dim(datExpr)[[1]]
dim(datExpr)

#计算软阈值
powers1=c(seq(1,10,by=1),seq(12,20,by=2))
RpowerTable=pickSoftThreshold(datExpr, powerVector=powers1)[[2]]
cex1=1
par(mfrow=c(1,2))
pdf("beta.pdf")
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],xlab="
Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     main = "Scale independence",
     type="n")
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], labels=powers1,cex=cex1,col="red")
abline(h=0.85,col="red")
plot(RpowerTable[,1], RpowerTable[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", 
     main = "Mean connectivity",
     type="n")
text(RpowerTable[,1], RpowerTable[,5], labels=powers1, cex=cex1,col="red")
dev.off()

#软阈值的检验
beta1=5 #根据自己的检验结果来选择
Connectivity=softConnectivity(datExpr,power=beta1)
pdf("scalefree.pdf",15,10)
par(mfrow=c(1,1))
scaleFreePlot(Connectivity, main=paste("soft threshold, power=",beta1), truncated=T)
dev.off()


#模块的检测
ConnectivityCut = 4440#最多能连接的基因数
ConnectivityRank = rank(-Connectivity)	
restConnectivity = ConnectivityRank <= ConnectivityCut
ADJrest = adjacency(datExpr[,restConnectivity], power=beta1)
dissTOM=TOMdist(ADJrest)
hierTOM = hclust(as.dist(dissTOM),method="average")
colorh1= cutreeStaticColor(hierTOM,cutHeight = 0.8, minSize = 3) #调节模块数量
pdf("module.pdf")
par(mfrow=c(2,1),mar=c(2,4,1,1))
plot(hierTOM, main="Cluster Dendrogram", labels=F, xlab="", sub="")
plotColorUnderTree(hierTOM,colors=data.frame(module=colorh1))
title("Module (branch) color")
dev.off()

#拓扑重叠热图
pdf("TOM.pdf")
TOMplot(dissTOM , hierTOM, colorh1, terrainColors=TRUE) 
dev.off()

#模块进行聚类
datME=moduleEigengenes(datExpr[,restConnectivity],colorh1)[[1]]
dissimME=1-(t(cor(datME, method="p")))/2
hclustdatME=hclust(dist(dissimME), method="average" )
pdf("modul_cluster.pdf")
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based on the module eigengenes of modules")
dev.off()

#模块与模块之间的相关性
pdf("modul_cor.pdf")
pairs(datME)
dev.off()


modul<-signif(cor(datME, use="p"), 2)
write.table(modul,"modul_cor.txt",sep="\t",quote=F)

#输出模块基因
datME=moduleEigengenes(datExpr,colorh1)[[1]]
color1=rep("grey",dim(datExpr)[[2]])
color1=as.character(colorh1)
datKME=signedKME(datExpr, datME)
datout=data.frame(datSummary, colorNEW=color1,datKME )
write.table(datout, "gene_module.xls", sep="\t", row.names=F,quote=F)

#将网络导出
exportNetworkToCytoscape(ADJrest,edgeFile="edge.txt",nodeFile="node.txt",threshold = 0.5) #根据网络大小调节

#模块与性状的相关性
inputclinial="clinical data.txt" #改为自己的临床数据
dataclinial = read.table(inputclinial,sep="\t",row.names=1,header=T,check.names=F,quote="!")
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr,colorh1)$eigengenes
MEsFemale = orderMEs(MEs0)
modul_clinical_cor = cor(MEsFemale, dataclinial, use = "p")
write.table(modul_clinical_cor,"module-clinial-cor.xls",sep="\t",quote=F)
modul_clinical_p = corPvalueStudent(modul_clinical_cor, nSamples)
write.table(modul_clinical_p,"modul-clinical-p.xls",sep="\t",quote=F)
textMatrix = paste(signif(modul_clinical_cor, 2), " (", signif(modul_clinical_p, 1), ")",sep = "")
dim(textMatrix) = dim(modul_clinical_cor)
pdf("modul-clinical.pdf",10,8)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = modul_clinical_cor, xLabels = names(dataclinial), 
  yLabels = names(MEsFemale),ySymbols = names(MEsFemale), colorLabels = FALSE, 
  colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, 
  cex.text = 1, zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()

#基因与性状的相关性
geneExp=datExpr
gene_clinial_cor = cor(geneExp, dataclinial, use = "p")
write.table(gene_clinial_cor,"gene_clinial_cor.xls",sep="\t",quote=F)
gene_clinical_p = corPvalueStudent(gene_clinial_cor, nSamples)
write.table(gene_clinical_p,"gene_clinical_p.xls",sep="\t",quote=F)


#GO-KEGG
GO
#安装R包
install.packages("ggplot2")

setwd("C:\\...") #改为自己的工作目录
inputfile="GO_vis.txt" #改为自己的文件
library(ggplot2)

go=read.table(inputfile,header=T,sep="\t")
pdf("go.pdf",15,10)
ggplot(data=go)+geom_bar(aes(x=Term, y=Count, fill=-log10(PValue)), stat='identity')+
  coord_flip() + scale_fill_gradient(low="blue", high = "red")+ 
  xlab("") + ylab("") + theme(axis.text.x=element_text(color="black", size=12),
  axis.text.y=element_text(color="black", size=12)) + 
  scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(expand=c(0,0))
dev.off()

KEGG
#安装软件
install.packages("RSQLite")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

library("org.Hs.eg.db")

setwd("C:\\...") #改为自己的工作目录

inputfile="module-gene-go.txt" #改为自己的文件名

gene_symbol=read.table(inputfile,sep="\t",check.names=F,header=F)
gene_name=as.vector(gene_symbol[,1])
geneID <- mget(gene_name, org.Hs.egSYMBOL2EG, ifnotfound=NA)
geneID <- as.character(geneID)
data=cbind(gene_symbol,entrezID=geneID)
write.table(data,"geneID.txt",sep="\t",quote=F,row.names=F)

#KEGG分析
install.packages("colorspace")
install.packages("stringi")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DOSE")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")

library(clusterProfiler)
gene<-read.table("geneID.txt",header=T,sep="\t")
kegg <- enrichKEGG(gene = gene$entrezID,organism ="human",pvalueCutoff = 0.05)
write.csv(kegg,"KEGG.csv",row.names =F)
#画图
pdf("kegg.pdf")
barplot(kegg, showCategory=20)
dev.off()

#气泡图
pdf("kegg2.pdf")
dotplot(kegg)
dev.off()


#网络图
pdf("kegg3.pdf",15,8)
library(DOSE)
data(geneList)
de <- names(geneList)[abs(geneList) > 2]
edo <- enrichDGN(de)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
dev.off()


#nomergram
#GSE1297-GSE36980-GSE28146
#1
setwd("D:\\...")
hubgene=read.table("hub.txt",header = T,sep = "\t")
genem=read.table("allgene.txt",header = T,sep = "\t")
resultge=merge(hubgene,genem,by="id")
write.table(resultge,"result.txt",sep = "\t",row.names = F,quote = F)

#2
library(rms)
library(ROCR)
setwd("D:\\...")

inputfile2="1297input.txt"
data1<- read.table(inputfile2,header=T,sep="\t")


ddist <- datadist(data1)
options(datadist="ddist") 

mylog<- glm(status~ATP2A2+ATP6V1D+CAP2+SYNJ1+GHITM,family=binomial(link = "logit"),
            data = data1)  

summary(mylog)
coefficients(mylog)
exp(coefficients(mylog))

#nom
mylog<-lrm(status~ATP2A2+ATP6V1D+CAP2+SYNJ1+GHITM,data=data1,x=T,y=T)


mynom<- nomogram(mylog, fun=plogis,fun.at=c(0.0001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.9999),
                 lp=F, funlabel="Risk of AD")

pdf("Nom3.pdf",15,8)
plot(mynom)
dev.off()
#C-index
Cindex <- rcorrcens(data1$status~predict(mylog))
Cindex
#Calibration
mycal<-calibrate(mylog,method="boot",B=1000)

pdf("Calibration.pdf")
plot(mycal,xlab="Nomogram-predicted probability of AD",ylab="Actual diagnosed AD (proportion)",sub=F)
dev.off()

#AUC
pre_rate<-predict(mylog)
ROC1<- prediction(pre_rate,data1$status) 
ROC2<- performance(ROC1,"tpr","fpr")
AUC <- performance(ROC1,"auc")

AUC
print(AUC,max=10000)
AUC<-0.9192 #要改为自己的AUC

pdf("ROC.pdf")
plot(ROC2,col="blue", xlab="False positive rate",ylab="True positive rate",lty=1,lwd=3,main=paste("AUC=",AUC))
abline(0,1,lty=2,lwd=3)
dev.off() 

#immune
#1
setwd("C:\\...")
mr<- read.table("mr.txt",header=T,sep="\t",check.names=F,row.names=1)
mr=t(mr)

pdf("barplot.pdf",10,8)
par(las=1,mar=c(8,4,4,15))
nub1= barplot(mr,col=rainbow(nrow(mr),s=0.9,v=0.9),
              yaxt="n",ylab="Relative Percent",xaxt="n")
nub2=axis(2,tick=F,labels=F)
axis(2,nub2,paste0(nub2*100,"%"))
axis(1,nub1,labels=F)
par(srt=60,xpd=T)
text(nub1,-0.02,colnames(mr),adj=1,cex=0.7);par(srt=0)
ytick2 = cumsum(mr[,ncol(mr)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],
       legend=rownames(mr),col=rainbow(nrow(mr),s=0.9,v=0.9),
       pch=16,bty="n",cex=1)
dev.off()

#2
library(pheatmap)
setwd("F:\\...")
data1=read.table("clinical.txt",header = T,sep = "\t")
data2=read.table("mr.txt",header = T,sep = "\t",check.names = F)
data3=merge(data1,data2,by="id")
write.table(data3,"data3.txt",sep = "\t",quote = F,row.names = F)

ht=read.table("ht.txt",header = T,sep = "\t",row.names = 1)
Group=c(rep("Con",2),rep("AD",2))    
names(Group)=colnames(ht)
Group=as.data.frame(Group)

pdf("heatmap.pdf",10,8)
pheatmap(ht, annotation=Group, 
         color = colorRampPalette(c("blue", "black", "red"))(50),
         cluster_cols =F)
dev.off()


#3
install.packages("corrplot")
install.packages("vioplot")

library(corrplot)
library(vioplot)
setwd("F:\\...")
mr=read.table("data3.txt",sep="\t",header=T,row.names=1,check.names=F)

mr=mr[,-2]
mr=mr[,-3]
mr=mr[,-9]
mr=mr[,-14]
mr=mr[,-15]
mr=mr[,-17]

pdf("corrplot2.pdf",12,10)              
corrplot(corr=cor(mr),method="color",order = "alphabet",tl.col="black",
         addCoef.col = "black",
         number.cex = 0.8,
         col=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()



