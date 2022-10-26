library(coda)

num_edges =  t(as.matrix(read.csv('Helper_functions/geweke_input.csv', header = F)))

# z_scores = geweke.diag(num_edges, frac1=0.5, frac2=0.5)$z
z_scores = geweke.diag(num_edges)$z
p_values = 2*pnorm(abs(z_scores), lower.tail = FALSE)

write.table(p_values, "Helper_functions/geweke_pvalues.csv", 
            row.names=FALSE, col.names=FALSE, sep=",")
