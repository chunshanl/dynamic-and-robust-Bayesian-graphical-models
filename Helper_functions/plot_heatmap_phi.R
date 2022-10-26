# Plot the heatmap of the posterior mean of Phi

library(GGally)   
library(mvtnorm)

Sigma <- as.matrix(read.csv('Phi_post.csv', header = F))
p <- dim(Sigma)[1]
Sigma <- Sigma + t(Sigma) + diag(p)
# Generate data
set.seed(1)
sample_df <- data.frame(rmvnorm(p*20000, mean=rep(0,p), sigma=Sigma))
temp = c("State_1", "State_2", "State_3", "State_4", "State_5", "State_6", "State_7")
colnames(sample_df) <- temp[1:p]

# Matrix of plots
p2 <- ggcorr(sample_df, label = TRUE, label_round = 2, low = "#EEEEEE" ,
             high = "#3B9AB2", limits = c(0,1), label_size = 4.5, legend.size = 13)
print(p2)
png(filename="post_phi.png")
dev.off()



