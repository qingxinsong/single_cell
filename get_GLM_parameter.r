#This script will generate slope and intercept in Poisson generalized linear model using ERCC level data
data=read.table('raw_data/ERCC_raw_reads/ERCC_level.txt',header=TRUE)
slope=matrix(rep(0,64),ncol=2)
row.names(slope)=colnames(data)[2:33]
colnames(slope)=c('intercept','slope')
for (i in 2:33) {
x=data[,1]
y=data[,i]
fit=glm(y~log2(x), family="poisson")
slope[i-1,1] <- fit$coefficients[1]
slope[i-1,2] <- fit$coefficients[2]
write.table(slope,'GLM_param.txt',sep='\t',quote = FALSE)
}
