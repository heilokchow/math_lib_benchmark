library("ggplot2")
library("ggpubr")

x <- read.csv("result.txt", sep = ",", header = FALSE)
x <- as.data.frame(x)
x <- x[which(x[,1]=="DGEMM"),2:4]
# x <- x[which(x[,1]=="DGESV"),2:4]
# x <- x[which(x[,1]=="DPOTRF"),2:4]
colnames(x) <- c("V1", "V2", "V3")

## Bar plot

X <- list()
Q <- list()
l <- c(1000, 2000, 4000, 8000)

for (i in 1:4) {
  X[[i]] <- x[which(x[,2]==l[i]),]
  X[[i]]$V1 <- paste(X[[i]]$V1)
  X[[i]]$V2 <- paste(X[[i]]$V2)
  Q[[i]] <- ggplot(X[[i]], aes(V2, V3)) + geom_bar(aes(fill = V1), width = 0.4, position = position_dodge(width=0.5), stat="identity")  +  
    theme(legend.position="right", legend.title = element_blank(),axis.title.x=element_blank(), axis.title.y=element_blank()) +
    scale_fill_discrete(labels = c("MKL", "OpenBLAS", "Eigen", "Eigen+MKL", "Eigen+OpenBLAS"))
}

ggarrange(Q[[1]], Q[[2]], Q[[3]], Q[[4]], ncol=4 , nrow=1, common.legend = TRUE, legend="bottom") 
  
## Curve plot

x[, 1] <- sapply(x[, 1], function(x) {ifelse(x==1, "MKL", "OpenBLAS")})
x$V3 <- x$V3^(1/3)
q <- ggplot(x, aes(V2, V3)) + geom_point(aes(colour = V1)) + geom_line(aes(colour = V1)) + xlab("dimension") + ylab("t^(1/3) (s)") +
  theme(legend.title = element_blank()) + scale_fill_discrete(labels = c("MKL", "OpenBLAS"))