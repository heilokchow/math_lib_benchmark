library("ggplot2")
library("ggpubr")

x <- read.csv("result.txt", sep = ",", header = FALSE)
x <- as.data.frame(x)

x <- x[which(x[,1]=="DPOTRF"),2:4]
colnames(x) <- c("V1", "V2", "V3")
x1 <- x[which(x[,2]==1000),]
x1$V1 <- paste(x1$V1)
x1$V2 <- paste(x1$V2)
q1 <- ggplot(x1, aes(V2, V3)) + geom_bar(aes(fill = V1), width = 0.4, position = position_dodge(width=0.5), stat="identity")  +  
  theme(legend.position="right", legend.title = 
          element_blank(),axis.title.x=element_blank(), 
        axis.title.y=element_blank()) + scale_fill_discrete(labels = c("MKL", "OpenBLAS", "Eigen", "Eigen+MKL", "Eigen+OpenBLAS"))

x2 <- x[which(x[,2]==2000),]
x2$V1 <- paste(x2$V1)
x2$V2 <- paste(x2$V2)
q2 <- ggplot(x2, aes(V2, V3)) + geom_bar(aes(fill = V1), width = 0.4, position = position_dodge(width=0.5), stat="identity")  +  
  theme(legend.position="right", legend.title = 
          element_blank(),axis.title.x=element_blank(), 
        axis.title.y=element_blank()) + scale_fill_discrete(labels = c("MKL", "OpenBLAS", "Eigen", "Eigen+MKL", "Eigen+OpenBLAS"))

x3 <- x[which(x[,2]==4000),]
x3$V1 <- paste(x3$V1)
x3$V2 <- paste(x3$V2)
q3 <- ggplot(x3, aes(V2, V3)) + geom_bar(aes(fill = V1), width = 0.4, position = position_dodge(width=0.5), stat="identity")  +  
  theme(legend.position="right", legend.title = 
          element_blank(),axis.title.x=element_blank(), 
        axis.title.y=element_blank()) + scale_fill_discrete(labels = c("MKL", "OpenBLAS", "Eigen", "Eigen+MKL", "Eigen+OpenBLAS"))

x4 <- x[which(x[,2]==8000),]
x4$V1 <- paste(x4$V1)
x4$V2 <- paste(x4$V2)
q4 <- ggplot(x4, aes(V2, V3)) + geom_bar(aes(fill = V1), width = 0.4, position = position_dodge(width=0.5), stat="identity")  +  
  theme(legend.position="right", legend.title = 
          element_blank(),axis.title.x=element_blank(), 
        axis.title.y=element_blank()) + scale_fill_discrete(labels = c("MKL", "OpenBLAS", "Eigen", "Eigen+MKL", "Eigen+OpenBLAS"))

ggarrange(q1, q2, q3, q4, ncol=4 , nrow=1, common.legend = TRUE, legend="bottom") 
 
