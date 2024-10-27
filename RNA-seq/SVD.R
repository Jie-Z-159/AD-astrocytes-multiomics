#SVD for RNA-seq data

rm(list=ls())
library(ggplot2)
library(dplyr)
#preprocess for data ----
femData <-read.csv("gene_fpkm_A.csv",row.names = 1,header=TRUE,sep=",")
femData<-femData[,-c(31:39)]
selectgenes<-order(exp.mad,decreasing = T)[1:10000]
femData.filtererd<-femData[selectgenes,]
mydata<-femData.filtererd
logmydata<-log10(mydata)
logmydata<-logmydata[is.finite(rowSums(logmydata)),]
z_scoredata<-scale(logmydata,center = TRUE, scale = TRUE)
z1<-z_scoredata[,2]
mean(z1)
sd(z1)
finitedata<-round(z_scoredata,5)
finitedata<-finitedata[is.finite(rowSums(finitedata)),]
#SVD----
svd_res<-svd(finitedata)
df_zscore.svd<-svd_res
df_zscore.svd.modes <- diag(df_zscore.svd$d)%*%t(df_zscore.svd$v)
df_singular_values <- data.frame(n=1:length(df_zscore.svd$d),singular_value=df_zscore.svd$d)

p1<-ggplot(df_singular_values,aes(x = n, y = singular_value)) + 
  geom_point(shape=1,color='#3299ff',size=1.5) +
  geom_line(color='#3299ff')+
  labs(y = 'Singular value', x = 'Singular value rank', title = 'Singular values of the astrocyte expression data')+
  theme_classic()+
  theme(plot.title=element_text(hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        text = element_text(size=12))

#mode1----
Modes_Comb <- df_zscore.svd.modes
nkeep <- 1
Time <- c(2, 4, 6, 9, 12)
# create empty matrix to hold t-test results
PT2Mat <- matrix(nrow=nkeep, ncol=5)
for (i in 1:nkeep) {
  
  XAD <- matrix(Modes_Comb[i,1:15], nrow=3)
  XWT <- matrix(Modes_Comb[i,16:ncol(Modes_Comb)], nrow=3)
  
  for (j in 1:5) {
    PT2Mat[i,j] <- t.test(XAD[,j], XWT[,j])$p.value
  }
  
  df <- data.frame(Time=Time, Group=rep(c('XAD', 'XWT'), each=5), Values=c(colMeans(XAD), colMeans(XWT)))
  
  p <- ggplot(df, aes(Time, Values, color=Group)) + 
    geom_point(size=2, shape=1) +
    geom_line(lwd=1) +
    labs(title=paste("Mode", i), x='month', y='Mean mode value') +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major.y = element_line(linetype = "dashed", color = "gray"),
          panel.grid.minor = element_blank(),
          text = element_text(size=12)) +
    scale_color_manual(values=c("red", "black"), labels = c("AD", "WT")) +
    guides(colour=guide_legend(title="Group")) +
    scale_x_continuous(breaks = c(2, 4, 6, 9, 12))
  
max_value <- max(df$Values)
padding <- 0.003 * max_value  
y <- max_value - padding   

  # Add t-test results to the plot
  for (j in 1:5) {
    max_value <- max(df$Values)
    padding <- 0.001 * max_value  
    y <- max_value - padding    
    p_value <- PT2Mat[i, j]  
    if (p_value < 0.01) {
      label <- "**"  
    } else {
      label <- " "  
    }
    p <- p + annotate("text", x = Time[j], y = y, label = label, size = 5,color="red")
  }
  if (i == nkeep) {
    xlab("month")
  }
  print(p)
}


