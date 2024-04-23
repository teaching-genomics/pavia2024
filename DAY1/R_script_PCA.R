rm(list = ls())

#load packages
library(tidyverse)
library(readxl)

#read in data
pca <- read_table2("./ecoli_pop.eigenvec", col_names = FALSE)
eigenval <- scan("./ecoli_pop.eigenval")

#remove duplicated column
pca <- pca[,-1]

# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

##Plot the percentage of variation explained by PCs

#convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()                  

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

#plot pca
b <- ggplot(pca, aes(PC1, PC2)) + geom_point(size = .5)
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))



