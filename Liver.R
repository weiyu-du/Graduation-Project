library(SingleCellExperiment)
library(cellassign)
library(scran)
library(sqldf)

#导入数据
dat = read.delim("Liver-counts.csv", sep=",",header=TRUE)
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
sce<-normalize(sce)
print(sce)

#标记基因
rho<- list(Hepatocyte = c("Acly", "Alb", "Apoa1", "Asl", "Ass1", "Cyp2e1"
                          , "Cyp2f2", "G6pc", "Glul", "Mup3", "Pck1"),
           Kupffer_cell = c("Clec4f", "Irf7", "Spic","Cd68"),
           B_cell = c("Cd22"),
           Natural_Killer_cell = c("Cxcr6")
           )
#优化后的标记基因
rho<- list(Hepatocyte = c("Acly", "Alb", "Apoa1", "Asl", "Ass1", "Cyp2e1"
                          , "Cyp2f2"),
           Kupffer_cell = c("Clec4f", "Irf7", "Spic","Cd68"),
           B_cell = c("Cd22"),
           Natural_Killer_cell = c("Cxcr6")
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

#热图绘制
pheatmap::pheatmap(cellprobs(fit),cluster_row = FALSE)
pheatmap::pheatmap(cellprobs(fit))


#导出数据
write.csv(Cell_Types,file="Liver_out.csv")
rownames(Cell_Probs)<-colnames(sce)
write.csv(Cell_Probs,file="Liver_out_probs.csv")

#再次筛选
Liver_true = read.delim("Liver.csv",sep=",",header=TRUE)

Liver_test = read.delim("Liver_out_probs.csv",sep=",",header=TRUE)

Liver <- sqldf("select * from Liver_true true,Liver_test test where test.Cell = true.Cell")

write.csv(Liver,file="Liver_out_out.csv")

#计算F1分数
f1_fun(Liver[,2],Liver[,4])

