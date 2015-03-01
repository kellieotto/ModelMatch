##### Make boxplots of estimation simulation results
##### last edited: March 1, 2014
##### Kellie Ottoboni

library(ggplot2)
library(reshape)
dat <- data.frame()
fin <- c("simulation1.Rdata","simulation2.Rdata","simulation3.Rdata")
design <- c("Weak Separation", "Strong Separation", "Skewed Separation")
outcome <- c(": Linear", ": Moderately Nonlinear", ": Quadratic")
for(f in 1:3){
  load(fin[f])
  for(i in 1:3){
    temp <- est.store[[i]][,-c(3,5,7,9)] # remove MMATT columns
    temp[,2:5] <- -temp[,2:5]
    temp <- melt(temp)
    temp[,1] <- paste("Outcome ", i, outcome[i], sep="")
    temp[,4] <- rep(design[f], nrow(temp))
    dat <- rbind(temp,dat)
  }
}
pdf("sim_estimate_boxplots.pdf", width = 12)
qplot(X2, value, data=dat, geom="boxplot", ylim = c(-5,6)) + facet_grid(V4~X1) + xlab("") + ylab("Avg Treatment Effect") + theme(axis.text.x=element_text(angle=-90))
dev.off()
