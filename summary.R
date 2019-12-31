library("ggplot2")
library("ggpubr")
library("foreach")
library("parallel")
library("iterators")
library("doParallel")
library("readxl")

setwd("C:/Users/micha/Documents/Course/c++/20181208 matrix_multi/matrix_multi/matrix_multi/math_lib_benchmark")
x <- read.csv("result.txt", sep = ",", header = FALSE)
x <- as.data.frame(x)
# x <- x[which(x[,1]=="DGEMM"),2:4]
x <- x[which(x[,1]=="DGESV"),2:4]
x <- x[which(x[,1]=="DPOTRF"),2:4]
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
  
## R benchmark

rep =5
nn <- c(10, 50, 100, 400, 800, 1200, 1600, 2000)

### DGEMM

ss <- rep(0, length(nn))
for (j in 1:length(nn)) {
  n <- nn[j]
  for (i in 1:rep) {
    A <- matrix(rnorm(n^2), nrow = n, ncol = n)
    B <- matrix(rnorm(n^2), nrow = n, ncol = n)
    t1 <- Sys.time()
    C = B %*% A
    t2 <- Sys.time()
    ss[j] = ss[j] + (t2 - t1)
    cat(n, "--------", i, "\n")
  }
  ss[j] = ss[j]/rep
}

for (i in 1:(length(nn)-1)) {
  if (ss[i] > ss[i+1]) break()
}
ss[(i+1):length(nn)] = ss[(i+1):length(nn)] * 60
ss[8] = ss[8]/60
ss <- ss[-1]

### DGESV

ss1 <- rep(0, length(nn))
for (j in 1:length(nn)) {
  n <- nn[j]
  for (i in 1:rep) {
    A <- matrix(rnorm(n^2), nrow = n, ncol = n)
    A <- t(A) %*% A 
    x <- matrix(rnorm(n), nrow = n, ncol = n)
    b <- A %*% x
    t1 <- Sys.time()
    x1 = solve(A, b)
    t2 <- Sys.time()
    ss1[j] = ss1[j] + (t2 - t1)
    cat(n, "--------", i, "\n")
  }
  ss1[j] = ss1[j]/rep
}

for (i in 1:(length(nn)-1)) {
  if (ss1[i] > ss1[i+1]) break()
}
ss1[(i+1):length(nn)] = ss1[(i+1):length(nn)] * 60
ss1[8] = ss1[8]/60
ss1 <- ss1[-1]

### DPORTF

ss2 <- rep(0, length(nn))
for (j in 1:length(nn)) {
  n <- nn[j]
  for (i in 1:rep) {
    A <- matrix(rnorm(n^2), nrow = n, ncol = n)
    A <- t(A) %*% A 
    t1 <- Sys.time()
    B <- chol(A)
    t2 <- Sys.time()
    ss2[j] = ss2[j] + (t2 - t1)
    cat(n, "--------", i, "\n")
  }
  ss2[j] = ss2[j]/rep
}

for (i in 1:(length(nn)-1)) {
  if (ss2[i] > ss2[i+1]) break()
}
ss2 <- ss2[-1]

## Read Matlab benchmark
ss <- rbind(ss, ss1, ss2)
mat <- read_excel("matlab.xlsx")
mat1 <- x
for (i in 1:7) {
  mat1[i, 3] <- mat[2, i]
  mat1[i, 1] <- "Matlab"
}

rp <- x
for (i in 1:7) {
  rp[i, 3] <- ss[2, i]
  rp[i, 1] <- "R"
}

x$V2 = "Eigen"

### Curve plot

x <- rbind(x, rp, mat1)
x$V4 <- x$V4^(1/3)
q1 <- ggplot(x, aes(V3, V4)) + geom_point(aes(colour = V2)) + geom_line(aes(colour = V2)) + ggtitle("DGEMM") +
  theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

q2 <- ggplot(x, aes(V3, V4)) + geom_point(aes(colour = V2)) + geom_line(aes(colour = V2)) + ggtitle("DGESV") +
  theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

q3 <- ggplot(x, aes(V3, V4)) + geom_point(aes(colour = V2)) + geom_line(aes(colour = V2)) + ggtitle("DPOTRF") +
  theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

ggarrange(q1, q2, q3, ncol=3 , nrow=1, common.legend = TRUE, legend="bottom") 