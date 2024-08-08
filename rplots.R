

library(DESeq2)
library(tidyverse)
library(forcats)
library(data.table)
library(fgsea)

###START of deseq
counts <-  as.matrix(read.csv('results/verse_concat_filtered.csv', row.names = 'gene'))
head(counts)

as_tibble(counts)
coldata <- data.frame(samples = colnames(counts), case = c(rep('CTL', 3), rep('KO', 3)), row.names='samples')
coldata$case <- as.factor(coldata$case)
coldata$case <- relevel(coldata$case, ref='CTL')
print(coldata)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~case)
dds <- DESeq(dds)
res <- results(dds, contrast=c('case', 'KO', 'CTL'))

res2 <- res[order(res$pvalue),]
res2 <- as_tibble(res2, rownames='geneid')

id2 <- read_delim('results/id2gene.txt', col_names=c('geneid', 'genenames'))

res2 <- filter(res2, padj < 0.05)
write_csv(res2, "DESEQ_results.csv")

#will use later
res3 <- res2 %>%
  left_join(id2, by='geneid')

#### DAVID Stuff
DAVID_res <- filter(res2)
print(res2)
DAVID_res <- res2 %>%
  left_join(id2, by='geneid')
DAVID_res[8]
write_csv(DAVID_res[8], "DAVID_res2.csv")

####Ranked l2fc
make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  
  labeled_results<-arrange(labeled_results,desc(log2FoldChange))
  x <- as_tibble(read.delim(id2gene_path,header = FALSE))
  x<-rename(x,"geneid" = V1)
  joined <- inner_join(labeled_results,x)
  print(joined)
  xvec <- pull(joined, V2 )

  logvec <- pull(joined, log2FoldChange)
  foldvec <- setNames(logvec, xvec)
  foldvec<-foldvec[!is.na(foldvec)]
  return(foldvec)
}
ranked <- make_ranked_log2fc(res2,'results/id2gene.txt')


#write_csv(ranked[1], "DAVID_res.csv")


####FGSEA 
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  x<-gmtPathways(gmt_file_path)
  rnk_list<-rnk_list[!is.na(rnk_list)]
  #x<-x[!is.na(x)]
  fgseaRes <- fgsea(x, 
                    rnk_list,
                    minSize  = min_size,
                    maxSize  = max_size)
  fgseaRes<-as_tibble(fgseaRes)
}


fgseares1 <- run_fgsea('c2.cp.reactome.v2023.2.Hs.symbols.gmt',ranked,15,500)

###FGSEA plot
#from rshiny
fgseares1 <- filter(fgseares1, padj < 0.2)
fgseares1 %>% arrange(padj)
fgseares1 %>%
mutate(pathway = forcats::fct_reorder(pathway, NES)) %>%
  
ggplot() +
geom_bar(aes(x=pathway, y=NES, fill = NES <  0), stat='identity') +
scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
theme_minimal() +
theme(legend.position = "none")+
ggtitle('FGSEA results with padj < 0.05 ') +
ylab('Normalized Enrichment Score (NES)') +
xlab('') +
coord_flip()
ggsave("FGSEA.png")


####PCA just using counts data


plot_pca <- function(data,pcx,pcy) {

  #pca
  pcar <- prcomp(t(data))
  newtib <- as_tibble(pcar$x)
  newtib$name <- c("CTL","CTL","CTL","KO","KO","KO")
  print(newtib)
  
  #variance
  summ <- summary(pcar)
  print(summ)
  xlabel <- paste(pcx, ": ",(summ$importance[2,pcx]*100), "% variance", sep = "" )
  ylabel <- paste(pcy, ": ",(summ$importance[2,pcy]*100), "% variance", sep = "" )
  
  plot1 <- ggplot(newtib, aes(x=!!sym(pcx),y=!!sym(pcy),color=name ) ) +
    geom_point()+
    ggtitle('PCA plot') +
    xlab(xlabel) +
    ylab(ylabel)
  
  return(plot1)
}

plot_pca(counts,"PC1","PC2")
ggsave("pcplot.png")






####Histogram

histoplot <- ggplot(res2, aes(x=log2FoldChange)) + geom_histogram(bins = 100) + labs(title="Histogram of DESEQ results")
ggsave("histoplot.png")





####VOLC

label_res <- function(deseq2_res, padj_threshold) {
  newtib <-tibble(genes = rownames(deseq2_res),deseq2_res)
  newtib <-mutate(newtib, volc_plot_status = ifelse(padj >= padj_threshold, "NS",ifelse(log2FoldChange>0,"UP","DOWN")), .after = "genes")
}


lab_res<-label_res(res3,0.05)
topten <- lab_res[1:10,]
topten


plot_volcano <- function(labeled_results) {
  #labeled_results<- drop_na(labeled_results)
  ggplot(topten, aes(x=log2FoldChange, y=-log10(padj),color = volc_plot_status, label=genenames)) + 
    
    geom_point(data = labeled_results)+
    geom_point(data = topten, size = 1, show.legend = FALSE) +
    labs(title="Volcano plot of DESeq2 differential expression results")
}
plot_volcano(lab_res)
ggsave("plotvolc.png")

