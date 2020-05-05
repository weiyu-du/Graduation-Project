library(SingleCellExperiment)
library(cellassign)
library(scran)
library(sqldf)

#导入数据
dat = read.delim("Brain_Neurons-counts.csv",sep=",",header=TRUE)
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
rho<- list(Oligodendrocyte = t(read.delim("Oligodendrocyte_Brain.csv",sep=",",header=TRUE)),
           Neuron = t(read.delim("Neuron_Brain.csv",sep=",",header=TRUE)),
           Endothelial_cell = t(read.delim("Endothelial cell_Brain.csv",sep=",",header=TRUE)),
           Bergmann_glial_cell = t(read.delim("Bergmann_glial_cell_Brain.csv",sep=",",header=TRUE))
)

rho<- list(Oligodendrocyte = t(read.delim("Oligodendrocyte_Brain.csv",sep=",",header=TRUE)),
           Neuron = t(read.delim("Neuron_Brain.csv",sep=",",header=TRUE)),
           Endothelial_cell = t(read.delim("/Endothelial cell_Brain.csv",sep=",",header=TRUE))
)
print(str(rho))
rho <- marker_list_to_mat(rho)
#rho <- rho[,-4]
print(rho)

#筛选
sce <- sce[intersect(rownames(rho), rownames(sce)),]

#查找不同基因
setdiff(rownames(rho),rownames(sce))

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
write.csv(Cell_Types,file="Brain_out.csv")
rownames(Cell_Probs)<-colnames(sce)
write.csv(Cell_Probs,file="Brain_out_probs.csv")

#再次筛选
Brain_true = read.delim("Brain.csv",sep=",",header=TRUE)

Brain_test = read.delim("Brain_out_probs.csv", sep=",",header=TRUE)

Brain <- sqldf("select * from Brain_true true,Brain_test test where test.Cell = true.Cell")

write.csv(Brain,file="Brain_out_out.csv")

#计算F1分数
f1_fun(Brain[,2],Brain[,4])
