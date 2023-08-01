a=read.table("C:\\Users\\Monishka Battula\\Documents\\GSE200146_counts.Cell.Lines.txt",header=T,row.names=1,sep='\t')
a
a=a[,-1]
a
sizefactors = colSums(a)
columnSums = apply(a,2,sum)
cpm = apply(a,2,function(x)(x/sum(x))*1000000)
View(cpm)
log_transform = function(cpm){
  cpm = log2(cpm + 1)
  return(cpm)
}

#matrix1=cpm[1:100,]
#matrix1

log2_cpm = log_transform(cpm)
log2_cpm

d= log2_cpm

calculate_zscore = function(d){
  gMean = apply(d, 1, mean)
  gsd = apply(d, 1, sd)
  z_scores <- d
  
  for(i in 1:nrow(d)){
    z_scores[i, 1:ncol(z_scores)] = (d[i, 1:ncol(d)] - gMean[i])/gsd[i]
  }
  return(z_scores)
}
m = calculate_zscore(d)
matrix1 = m[!rowSums(is.na(m)),]
dim(matrix1)

mat1 = matrix1[1:100,]
dim(mat1)

library(ComplexHeatmap)
library(circlize)
Heatmap( mat1,col=colorRamp2(c(-10,0,5),c("orange","white","violet")))

#calculate the variance for all genes 
variance = apply(cpm, 1, var)
variance

#take top 100 genes with most variance
sortmat=sort(variance,decreasing=TRUE)
sortmat
length(sortmat)
ma=sortmat[1:100]
ma

matx1=names(ma)
var_CPM = cpm[matx1, ]

z_score1=calculate_zscore(var_CPM)

library(ComplexHeatmap)
library(circlize)

temp=var_CPM
z_score1=z_score1[,order(colnames(z_score1))]
Heatmap(z_score1)


meta_data=read.csv("C:\\Users\\Monishka Battula\\Downloads\\meta.csv")
meta_data
head(meta_data)
ha=(HeatmapAnnotation(cell_line=meta_data$cell_line,age=meta_data$Age,gender=meta_data$Gender))
Heatmap(z_score1,top_annotation = ha)
