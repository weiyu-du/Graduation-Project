library(SingleCellExperiment)
library(cellassign)
library(scran)
library(sqldf)

#导入数据
dat = read.delim("Heart-counts.csv", sep=",",header=TRUE)
dim(dat)


#删除第一列当列名
rownames(dat)<-dat[,1]
dat<-dat[,-1]
dat[1:5,1:5]

#对数转换
log_dat<-apply(dat,2,function(x) log2(x+1))
log_dat[1:5,1:5]



#创建SingleCellExperiment
sce<-SingleCellExperiment(assays = list(counts = as.matrix(log_dat)),
                          rowData = data.frame(gene = rownames(log_dat)),
                          colData = data.frame(cell = colnames(log_dat)))
#标准化
sce<-normalize(sce)
print(sce)

#标记基因
rho<- list(Fibroblast = c("Col1a2", "Col3a1", "Fbln2", "Fstl1", "Gsn", "Sparc"
                          , "Vim","Thy1","Ckap4"),
           Endothelial_cell = c("Ednrb", "Egfl7", "Emcn", "Epas1", "Fabp4"
                                , "Flt1", "Pecam1")
           )
print(str(rho))
rho <- marker_list_to_mat(rho)
#rho <- rho[,-4]
print(rho)

#筛选
sce <- sce[intersect(rownames(rho), rownames(sce)),]

#大小因子
s1<-computeSumFactors(as.matrix(log_dat))

#具体计算
fit <- cellassign(exprs_obj = sce,
                  marker_gene_info = rho,
                  s = s1, 
                  learning_rate = 1e-2,
                  shrinkage = TRUE,
                  verbose = FALSE)
print(fit)

Cell_Types<-celltypes(fit)
Cell_Probs<-cellprobs(fit)

pheatmap::pheatmap(cellprobs(fit),cluster_row = FALSE)

#导出数据
write.csv(Cell_Types,file="Heart_out.csv")
rownames(Cell_Probs)<-colnames(sce)
write.csv(Cell_Probs,file="/Heart_out_probs.csv")


#再次筛选
Heart_true = read.delim("Heart.csv",sep=",",header=TRUE)

Heart_test = read.delim("Heart_out_probs.csv",sep=",",header=TRUE)

Heart <- sqldf("select * from Heart_true true,Heart_test test where test.Cell = true.Cell")

write.csv(Heart,file="Heart_out_out.csv")

#计算F1分数
f1_fun(Heart[,2],Heart[,4])
