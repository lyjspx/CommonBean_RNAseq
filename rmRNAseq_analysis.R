# identify DEGs with rmRNAseq ------------------------------
library(rmRNAseq)
rm(list = ls())
counts <- read.csv('PvulgarisUI111.v1.1.COUNTS.csv',row.names = 1)
counts[1:5,]
counts <- counts[apply(counts, 1,function(x) sum(x > 5))>23,] #quality control
design <- read.csv('design.csv')[,c(1:4,6:9)]
design0 <- as.matrix(design)

meta.data <- read.csv('PvulgarisUI111.Metadata.csv')
Subject <- as.numeric(factor(paste0(meta.data$tissue,"_",meta.data$replicate)))
Time <- design$condition
Nboot <- 100
ncores <- 20
print.progress <- T
C.matrix0 <- list()
C.matrix0[[1]] <- limma::makeContrasts(condition, levels = design0)
C.matrix0[[2]] <- limma::makeContrasts(Basal-Roots, levels = design0)
names(C.matrix0) <- c('condition','Basal_Roots')
counts.test <- counts
date()
bootCAR1 <- rmRNAseq:::TC_CAR1(counts = counts.test, design = design0, Subject = Subject,
                              Time = Time, C.matrix = C.matrix0, Nboot = Nboot,
                              ncores = ncores, print.progress = print.progress)
date()

sum(bootCAR1$pqvalue$qv$condition < 0.01)
sum(bootCAR1$pqvalue$qv$Basal_Roots < 0.01)

output.genes <- cbind.data.frame(gene=row.names(counts),
                                 condition.pvalue=bootCAR1$pqvalue$pv$condition,
                                 Basal_Roots.pvalue=bootCAR1$pqvalue$pv$Basal_Roots,
                                 condition.qvalue=bootCAR1$pqvalue$qv$condition,
                                 Basal_Roots.qvalue=bootCAR1$pqvalue$qv$Basal_Roots)
dim(output.genes)
write.csv(output.genes,file = 'output.p.q.values.csv')

# PCA ----------------------------------------------------------------------
count.pca <- prcomp(t(counts), center = TRUE,scale. = TRUE)
count.pca$x
summary.pca <- summary(count.pca)

PC.summary.plot <- 
  rbind.data.frame(cbind.data.frame(name=1:24,
                                    portion=as.numeric(summary.pca$importance[2,]),
                                    group='prop'),
                   cbind.data.frame(name=1:24,
                                    portion=as.numeric(summary.pca$importance[3,]),
                                    group='cumu prop'))

library('ggplot2')
p<-ggplot(PC.summary.plot, aes(x=name, y=portion, group=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group))
p


PC.plot <- as.data.frame(count.pca$x[,1:2])
PC.plot <- cbind.data.frame(PC.plot,tissue=meta.data$tissue,Time=meta.data$condition)
PC.plot$Time <- factor(PC.plot$Time)

ggplot(PC.plot, aes(x=PC1, y=PC2, color=tissue, shape=Time))+
  geom_point(size=3)+
  geom_text(label=row.names(PC.plot),vjust=-1,size=3)
