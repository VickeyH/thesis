library(sva)
library(edgeR)
library(limma)
library(data.table)
library(rtracklayer)
library(parallel)

filelist <- list.files('/well/got2d/apayne/GTEx_v7/gene_counts/')
gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode <- gencode[gencode$seqid!='MT',]

phenotype <- read.table('/well/got2d/GTEx_v7/sample_annotations/GTEx_Analysis_2016-01-15_v7_SubjectPhenotypesDS.txt', header=TRUE, stringsAsFactors=FALSE, quote="\"", sep="\t", comment.char="")

filelist <- filelist[filelist!="Cells_Leukemia.reads.txt"]

expr_filter <- function(currtiss) {
  #Data cleaning
  gc()
  expr <- fread(paste('/well/got2d/apayne/GTEx_v7/gene_counts/', currtiss, sep=''), sep='\t', stringsAsFactors=FALSE, data.table=FALSE)
  rownames(expr) <- expr[,1]
  expr <- expr[,-c(1, 2)]
  expr <- expr[rownames(expr)%in%gencode$gene_id,]
  expr <- as.matrix(expr)

  sampids <- c()

  for (j in 1:ncol(expr)) {
    sampids <- c(sampids, paste(strsplit(colnames(expr), '-', fixed=TRUE)[[j]][1], strsplit(colnames(expr), '-', fixed=TRUE)[[j]][2], sep='-'))
  }

  sex <- phenotype[match(sampids, phenotype$SUBJID),'SEX']

  if (ncol(expr)<50) return(NULL)
  if (length(unique(sex))==1) return(NULL)

  genereads <- expr[rowSums(expr>6)>=10,]
  gc()
  dge <- DGEList(counts=genereads)
  dge <- dge[rowSums(cpm(dge)>1, na.rm=TRUE)>=10,]

  cpms <- as.data.frame(cpm(dge))
  cpms <- cbind(data.frame(gene_id=rownames(cpms), stringsAsFactors=FALSE), cpms)
  write.table(cpms, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t', file=paste('/well/got2d/apayne/GTEx_v7/filtered_final_expr_data_2sex/cpms/', substring(currtiss, 1, nchar(currtiss)-10), '.cpms', sep=''))

  dge$counts <- dge$counts + 1

  design <- model.matrix(~ as.factor(sex))
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot=F)$E

  voomed <- cbind(data.frame(gene_id=rownames(v), stringsAsFactors=FALSE), as.data.frame(v))
  write.table(voomed, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/GTEx_v7/filtered_final_expr_data_2sex/v_transformed/', substring(currtiss, 1, nchar(currtiss)-10), '.vnorm', sep=''))

  svars <- sva(dat=v, mod=design)$sv
  colnames(svars) <- paste('SVA', 1:ncol(svars), sep='')

  covars <- cbind(data.frame(sampids=sampids, stringsAsFactors=FALSE), svars, data.frame(sex=sex, stringsAsFactors=FALSE))
  write.table(covars, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, file=paste('/well/got2d/apayne/GTEx_v7/filtered_final_expr_data_2sex/final_covariates/', substring(currtiss, 1, nchar(currtiss)-10), '.covars', sep=''))
}

mclapply(filelist, FUN=expr_filter, mc.cores=27)


######################
# Tissue specificity #
######################
library(data.table)
library(rtracklayer)
library(parallel)
library(gtools)

######################################
# Tissues with DE per unique DE gene #
######################################
multiplot <- function(..., plotlist=NULL, file, cols=1, colprops=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow=1, ncol=cols, widths=colprops)))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

filelist <- list.files('/well/got2d/apayne/GTEx_v7/filtered_final_expr_data_2sex/v_transformed/')
filelist <- substring(filelist, 1, nchar(filelist)-6)
gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode <- gencode[gencode$seqid!='MT',]
gencode_genes <- subset(gencode, type=='gene')

vooms <- list()
cpms <- list()

for (i in filelist) {
  vooms[[i]] <- fread(paste('/well/got2d/apayne/GTEx_v7/filtered_final_expr_data_2sex/v_transformed/', i, '.vnorm', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  cpms[[i]] <- fread(paste('/well/got2d/apayne/GTEx_v7/filtered_final_expr_data_2sex/cpms/', i, '.cpms', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
}

tvoom <- list()

for (i in filelist) {
  currvoom <- vooms[[i]]
  rownames(currvoom) <- currvoom[,1]
  currvoom <- currvoom[,-1]
  tvoom[[i]] <- t(currvoom)
}

fullvoommatrix <- tvoom[[1]]

for (i in filelist[2:length(filelist)]) {
  gc()
  fullvoommatrix <- smartbind(fullvoommatrix, tvoom[[i]])
  print(i)
}

write.table(fullvoommatrix, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t', file='/well/got2d/apayne/thesis/chapter3/all_genes_voom_values_all_samples.txt')





tisscounts <- matrix(0, nrow=nrow(gencode_genes), ncol=length(filelist)+1)
rownames(tisscounts) <- gencode_genes$gene_id
colnames(tisscounts) <- c(filelist, 'number_of_tissues')

#exprmeans <- matrix(0, nrow

#expr_metrics <- matrix(NA, nrow=nrow(gencode_genes), ncol=5)
#expr_metrics <- as.numeric(expr_metrics)
#rownames(expr_metrics) <- gencode_genes$gene_id
#colnames(expr_metrics) <- c('mean', 'median', 'min', 'max', 'var')


for (i in filelist) {
  tisscounts[match(vooms[[i]][,1], rownames(tisscounts)),i] <- 1
  print(i)
}

tisscounts[,'number_of_tissues'] <- rowSums(tisscounts)
tisscounts <- tisscounts[tisscounts[,44]>0,]

nrow(tisscounts) #number of expressed genes overall

tisscounts_table <- data.frame(gene_id=rownames(tisscounts), gene_type=gencode[match(rownames(tisscounts), gencode$gene_id),'gene_type'], stringsAsFactors=FALSE)
table(tisscounts_table$gene_type)

#tisscounts_expressed <- tisscounts[tisscounts[,'number_of_tissues']>=1,]

#tisscounts_uniquegenes <- tisscounts[tisscounts[,'number_of_tissues']==1,-44]
#sort(colSums(tisscounts_uniquegenes))

#sum(tisscounts_expressed[,'number_of_tissues']==43)

#library(ggplot2)

#pdf('/well/got2d/apayne/thesis/chapter3/tissue_expr_dist_full.pdf', width=5.9, height=5.9)
#ggplot(as.data.frame(tisscounts_expressed)) +
#  geom_histogram(aes(x=number_of_tissues),
#    bins=43,
#    col='black',
#    fill='darkorange'
#  ) +
#  labs(x='Number of expressed tissues',
#    y='Count'
#  )
#dev.off()

#############################
# SEPARATE BY LNCRNA CODING #
#############################
library(gridExtra)

lncrna <- c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense",  "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "macro_lncRNA", "bidirectional_promoter_lncrna")
gencode_genes_lncrna <- subset(gencode_genes, gene_type%in%lncrna)
gencode_genes_coding <- subset(gencode_genes, gene_type=='protein_coding')

tisscounts_expressed_lncrna <- tisscounts_expressed[rownames(tisscounts_expressed)%in%gencode_genes_lncrna$gene_id,]
tisscounts_expressed_coding <- tisscounts_expressed[rownames(tisscounts_expressed)%in%gencode_genes_coding$gene_id,]
nrow(tisscounts_expressed_lncrna)
nrow(tisscounts_expressed_coding)

sum(tisscounts_expressed_lncrna[,44]==43)
sum(tisscounts_expressed_coding[,44]==43)

prop.test(matrix(c(691, 6088-691, 11195, 18059-11195), nrow=2, byrow=TRUE))

sum(tisscounts_expressed_lncrna[,44]<=10)
sum(tisscounts_expressed_coding[,44]<=10)

prop.test(matrix(c(3156, 6088-3156, 1717, 18059-1717), nrow=2, byrow=TRUE))

sum(tisscounts_expressed_lncrna[,44]==1)
sum(tisscounts_expressed_coding[,44]==1)

prop.test(matrix(c(966, 6088-966, 303, 18059-303), nrow=2, byrow=TRUE))


pdf('/well/got2d/apayne/thesis/chapter3/tissue_expr_dist_split.pdf', width=5.9, height=8.29)
p1 <- ggplot(as.data.frame(tisscounts_expressed_lncrna)) +
  geom_histogram(aes(x=number_of_tissues),
    bins=43,
    col='black',
    fill='darkblue'
  ) +
  labs(y='lncRNA count')+
  theme(axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank()
  )

p2 <- ggplot(as.data.frame(tisscounts_expressed_coding)) +
  geom_histogram(aes(x=number_of_tissues),
    bins=43,
    col='black',
    fill='darkorange'
  ) +
  labs(x='Number of expressed tissues',
    y='Coding gene count'
  )
grid.arrange(p1, p2, nrow=2)
dev.off()


gencode_genes_lincrna <- subset(gencode_genes, gene_type=='lincRNA')
gencode_genes_antisense <- subset(gencode_genes, gene_type=='antisense')
tisscounts_expressed_lincrna <- tisscounts_expressed[rownames(tisscounts_expressed)%in%gencode_genes_lincrna$gene_id,]
tisscounts_expressed_antisense <- tisscounts_expressed[rownames(tisscounts_expressed)%in%gencode_genes_antisense$gene_id,]

nrow(tisscounts_expressed_lincrna)
nrow(tisscounts_expressed_antisense)

nrow(tisscounts_expressed_lincrna[tisscounts_expressed_lincrna[,44]==43,])
nrow(tisscounts_expressed_antisense[tisscounts_expressed_antisense[,44]==43,])

prop.test(matrix(c(258, 2916-258, 294, 2482-294), nrow=2, byrow=TRUE))

nrow(tisscounts_expressed_lincrna[tisscounts_expressed_lincrna[,44]==1,])
nrow(tisscounts_expressed_antisense[tisscounts_expressed_antisense[,44]==1,])

prop.test(matrix(c(538, 2916-538, 358, 2482-358), nrow=2, byrow=TRUE))



pdf('/well/got2d/apayne/thesis/chapter3/tissue_expr_dist_linc_as.pdf', width=5.9, height=8.29)
p1 <- ggplot(as.data.frame(tisscounts_expressed_lincrna)) +
  geom_histogram(aes(x=number_of_tissues),
    bins=43,
    col='black',
    fill='darkred'
  ) +
  labs(y='lincRNA count')+
  theme(axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank()
  )

p2 <- ggplot(as.data.frame(tisscounts_expressed_antisense)) +
  geom_histogram(aes(x=number_of_tissues),
    bins=43,
    col='black',
    fill='darkgreen'
  ) +
  labs(x='Number of expressed tissues',
    y='Antisense lncRNA count'
  )
grid.arrange(p1, p2, nrow=2)
dev.off()

testdata <- matrix(c(table(tisscounts_expressed_lincrna[,44]), table(tisscounts_expressed_antisense[,44])), byrow=TRUE, nrow=2)
chisq.test(testdata[,-c(1, 2, 42, chisq.test(testdata[,-c(1, 2, 42, 43)])
43)])

##############################
# SPECIFICITY BY EACH TISSUE #
##############################

tisscounts_expressed_uniques <- tisscounts_expressed[tisscounts_expressed[,'number_of_tissues']==1,]
tissue_uniques <- colSums(tisscounts_expressed_uniques)[-44]

samp_sizes <- c()
nexpr <- c()
nlncs <- c()

for (i in 1:length(vooms)) {
  samp_sizes <- c(samp_sizes, ncol(vooms[[i]])-1)
  nexpr <- c(nexpr, nrow(vooms[[i]]))
  nlncs <- c(nlncs, sum(vooms[[i]][,1]%in%subset(gencode, gene_type%in%lncrna)$gene_id))
}

names(nexpr) <- names(nlncs) <- names(samp_sizes) <- names(vooms)

unique_plots <- data.frame(tissue=names(tissue_uniques), stringsAsFactors=FALSE)
unique_plots$n <- samp_sizes[match(unique_plots[,1], names(samp_sizes))]
unique_plots$nexpr <- nexpr[match(unique_plots[,1], names(nexpr))]
unique_plots$counts <- tissue_uniques[match(unique_plots[,1], names(tissue_uniques))]
cor(unique_plots[,2], unique_plots[,3])
corline <- lm(unique_plots[,3] ~ unique_plots[,2])

#pdf('/well/got2d/apayne/thesis/chapter3/unique_genes_per_tissue.pdf', width=5.9, height=5.9)
#ggplot(unique_plots, aes(x=n, y=counts)) +
#  geom_point() +
#  geom_abline(slope=summary(corline)$coefficients[2,1], intercept=summary(corline)$coefficients[1,1], col='red') +
#  labs(y='Count', x='Sample size')
#dev.off()

tisscounts_expressed_uniques_lncrna <- tisscounts_expressed_uniques[rownames(tisscounts_expressed_uniques)%in%gencode_genes_lncrna$gene_id,]
tisscounts_expressed_uniques_coding <- tisscounts_expressed_uniques[rownames(tisscounts_expressed_uniques)%in%gencode_genes_coding$gene_id,]
tissue_uniques_lncrna <- colSums(tisscounts_expressed_uniques_lncrna)[-44]
tissue_uniques_coding <- colSums(tisscounts_expressed_uniques_coding)[-44]

unique_plots_lncrna <- data.frame(tissue=names(tissue_uniques_lncrna), stringsAsFactors=FALSE)
unique_plots_coding <- data.frame(tissue=names(tissue_uniques_coding), stringsAsFactors=FALSE)
unique_plots_lncrna$n <- samp_sizes[match(unique_plots_lncrna[,1], names(samp_sizes))]
unique_plots_coding$n <- samp_sizes[match(unique_plots_coding[,1], names(samp_sizes))]
unique_plots_lncrna$counts <- tissue_uniques_lncrna[match(unique_plots_lncrna[,1], names(tissue_uniques_lncrna))]
unique_plots_coding$counts <- tissue_uniques_coding[match(unique_plots_coding[,1], names(tissue_uniques_coding))]

unique_plots_lncrna$total <- nlncs[match(unique_plots_lncrna$tissue, names(nlncs))]
prop.test(matrix(c(unique_plots_lncrna$counts, unique_plots_lncrna$total-unique_plots_lncrna$counts), byrow=FALSE, nrow=nrow(unique_plots_lncrna)))

library(ggplot2)
library(ggrepel)

unique_plots_lncrna$label <- ''
unique_plots_lncrna[unique_plots_lncrna$counts>50,'label'] <- unique_plots_lncrna[unique_plots_lncrna$counts>50,'tissue']
corline_lncrna <- lm(unique_plots_lncrna[,3] ~ unique_plots_lncrna[,2])

pdf('/well/got2d/apayne/thesis/chapter3/unique_lncrnas_per_tissue.pdf', width=5.9, height=5.9)
ggplot(unique_plots_lncrna, aes(x=n, y=counts)) +
  geom_point() +
  geom_text_repel(aes(label=label), size=2.5) +
  geom_abline(slope=summary(corline_lncrna)$coefficients[2,1], intercept=summary(corline_lncrna)$coefficients[1,1], col='red') +
  labs(y='Count', x='Sample size')
dev.off()

cor.test(unique_plots_lncrna[,2], unique_plots_lncrna[,3])
cor.test(unique_plots_lncrna[,2], unique_plots_lncrna[,3]/unique_plots_lncrna[,4])


cor.test(unique_plots_coding[,2], unique_plots_coding[,3])
corline_lncrna <- lm(unique_plots_lncrna[,3] ~ unique_plots_lncrna[,2])
corline_coding <- lm(unique_plots_coding[,3] ~ unique_plots_coding[,2])


plotout <- data.frame(tissue=unique_plots[,1], n=unique_plots[,2], expressed=unique_plots[,3], unique=unique_plots[,4], lncs=unique_plots_lncrna[,3], coding=unique_plots_coding[,3], stringsAsFactors=FALSE)
write.table(plotout, col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t&\t', eol='\\\\\n\\hline\n', file='/well/got2d/apayne/thesis/chapter3/uniquelncrnatable.txt')


######################################
# Expression quantification patterns #
######################################
library(data.table)
library(rtracklayer)
library(parallel)
library(ggplot2)
library(plyr)

filelist <- list.files('/well/got2d/apayne/GTEx_v7/filtered_final_expr_data_2sex/v_transformed/')
filelist <- substring(filelist, 1, nchar(filelist)-6)
gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode <- gencode[gencode$seqid!='MT',]
gencode_genes <- subset(gencode, type=='gene')
lncrna <- c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense",  "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "macro_lncRNA", "bidirectional_promoter_lncrna")
gencode_genes_lncrna <- subset(gencode_genes, gene_type%in%lncrna)
gencode_genes_coding <- subset(gencode_genes, gene_type=='protein_coding')

vooms <- list()
cpms <- list()

for (i in filelist) {
  vooms[[i]] <- fread(paste('/well/got2d/apayne/GTEx_v7/filtered_final_expr_data_2sex/v_transformed/', i, '.vnorm', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  cpms[[i]] <- fread(paste('/well/got2d/apayne/GTEx_v7/filtered_final_expr_data_2sex/cpms/', i, '.cpms', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
}

fullvoom <- fread('/well/got2d/apayne/thesis/chapter3/all_genes_voom_values_all_samples.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

voomrank <- list()

for (i in filelist) voomrank[[i]] <- vooms[[i]][order(rowMedians(as.matrix(vooms[[i]][,-1])), decreasing=TRUE),]

voomtops <- data.frame(tissue=filelist, stringsAsFactors=FALSE)
voomtops$g1 <- unlist(lapply(filelist, FUN=function(x) gencode_genes[match(voomrank[[x]][1,1], gencode_genes$gene_id),'gene_name']))
voomtops$g2 <- unlist(lapply(filelist, FUN=function(x) gencode_genes[match(voomrank[[x]][2,1], gencode_genes$gene_id),'gene_name']))
voomtops$g3 <- unlist(lapply(filelist, FUN=function(x) gencode_genes[match(voomrank[[x]][3,1], gencode_genes$gene_id),'gene_name']))

write.table(voomtops, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t&\t', eol='\\\\\n\\hline\n', file='/well/got2d/apayne/thesis/chapter3/topexpressedgenes.txt')


genemedians <- data.frame(medians=colMedians(as.matrix(fullvoom), na.rm=TRUE), stringsAsFactors=FALSE)
rownames(genemedians) <- colnames(fullvoom)
genemedians$gene_type <- gencode_genes[match(rownames(genemedians), gencode_genes$gene_id),'gene_type']
genemedians$gene_type_simple <- 'other'
genemedians[genemedians$gene_type=='protein_coding','gene_type_simple'] <- 'Coding'
genemedians[genemedians$gene_type%in%lncrna,'gene_type_simple'] <- 'lncRNA'

pdf('/well/got2d/apayne/thesis/chapter3/expression_voom_dist.pdf', width=5.9, height=5.9)
p1 <- ggplot(genemedians, aes(medians)) +
  geom_density() +
  labs(x='Median expression',
    y='Density'
  ) +
  coord_cartesian(xlim=c(-5, 11.5), ylim=c(0, 0.5))

p2 <- ggplot(subset(genemedians, gene_type_simple=='lncRNA'), aes(medians)) +
  geom_density() +
  labs(x='lncRNA expression',
    y='Density'
  ) +
  coord_cartesian(xlim=c(-5, 11.5), ylim=c(0, 0.5))


p3 <- ggplot(subset(genemedians, gene_type_simple=='Coding'), aes(medians)) +
  geom_density() +
  labs(x='Coding gene expression',
    y='Density'
  ) +
  coord_cartesian(xlim=c(-5, 11.5), ylim=c(0, 0.5))

lay <- matrix(c(NA, 1, 1, NA, NA, 1, 1, NA, 2, 2, 3, 3, 2, 2, 3, 3), byrow=TRUE, nrow=4)

grid.arrange(p1, p2, p3, layout_matrix=lay)
dev.off()

ks.test(genemedians[genemedians$gene_type_simple=='lncRNA','medians'], mean=mean(genemedians[genemedians$gene_type_simple=='lncRNA','medians']), sd=sd(genemedians[genemedians$gene_type_simple=='lncRNA','medians']), y='pnorm', alternative='two.sided')











