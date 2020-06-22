#This script will generate PCA figure using normalized gene levels
sc_exp=read.table('Normalized_gene_abundance.txt',header=TRUE)[,2:33]
sc_exp_log=log2(sc_exp+1)
sc_exp_mean_CC=rowMeans(sc_exp_log[,1:32])
genes_pass_filter= rownames(sc_exp_log)[sc_exp_mean_CC>0]
sc_exp_filter = sc_exp_log[rownames(sc_exp_log) %in% genes_pass_filter, ]
sc_pca=prcomp(t(sc_exp_filter[]),scale. = FALSE,center = FALSE,tol = 0)
pca_matrix=data.frame(sc_pca$x)
pca_matrix$Cell_type=rep(c('ECd','ECt','CCd','CCt'),each=8)
p=ggplot(pca_matrix,aes(PC1,PC2,color=Cell_type))+geom_point(size=1)
p=p+theme(panel.background = element_rect(fill='white'),axis.line.x=element_line(color='black'),
        axis.line.y=element_line(color='black'),legend.title=element_blank(),legend.text=element_text(size=11,family="sans"), axis.title=element_text(size=11,family="sans"),axis.text=element_text(size=11,color='black'),legend.position= c(0.86,0.78), legend.box =  "vertical")
p
ggsave("pca.png", width = 9, height = 7, units = "cm")

