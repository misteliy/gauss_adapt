data <- read.table("fort.111")
col <- rainbow(length(data))
plot(data$V1,data$V2,type='l')
for (i in 1:length(data)-1)
{
lines(data[,1],data[,i+1],type='l',col=col[i])
}
