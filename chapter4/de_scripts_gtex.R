library(sva)
library(edgeR)
library(limma)
library(data.table)
library(WGCNA)
library(rtracklayer)

filelist <- list.files('/well/got2d/apayne/GTEx_v7/gene_counts/')
gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode <- gencode[!gencode$seqid%in%c('MT'),]

phenotype <- read.table('/well/got2d/GTEx_v7/sample_annotations/GTEx_Analysis_2016-01-15_v7_SubjectPhenotypesDS.txt', header=TRUE, stringsAsFactors=FALSE, quote="\"", sep="\t", comment.char="")

filelist <- filelist[filelist!="Cells_Leukemia.reads.txt"]

worker <- function(currtiss) {
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
  t2d <- phenotype[match(sampids, phenotype$SUBJID),'MHT2D']
  keeps <- t2d!=99

  sex <- sex[keeps]
  t2d <- t2d[keeps]
  expr <- expr[,keeps]

  if (ncol(expr)<50) return(NULL)

  genereads <- expr[rowSums(expr>6)>=10,]
  gc()
  dge <- DGEList(counts=genereads)
  dge <- dge[rowSums(cpm(dge)>1, na.rm=TRUE)>=10,]
  dge$counts <- dge$counts + 1

  if (length(unique(sex))!=1) {
    design <- model.matrix(~ as.factor(sex) + as.factor(t2d))
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot=F)$E

    svars <- sva(dat=v, mod=design)$sv

    colnames(svars) <- paste('SVA', 1:ncol(svars), sep='')

    design <- model.matrix(~ svars + as.factor(sex) + as.factor(t2d))

    colnames(design) <- c('Intercept', paste('SVA', 1:ncol(svars), sep=''), 'sex', 't2d')
    fit <- lmFit(v, design)
    fit2 <- eBayes(fit)
    output <- topTable(fit2, coef=ncol(design), number=nrow(v))
    output$gene_id <- rownames(output)
    output$gene_name <- gencode[match(output$gene_id, gencode$gene_id),'gene_name']
  } else {
    design <- model.matrix(~ as.factor(t2d))
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot=F)$E

    svars <- sva(dat=v, mod=design)$sv

    colnames(svars) <- paste('SVA', 1:ncol(svars), sep='')

    design <- model.matrix(~ svars + as.factor(t2d))

    colnames(design) <- c('Intercept', paste('SVA', 1:ncol(svars), sep=''), 't2d')
    fit <- lmFit(v, design)
    fit2 <- eBayes(fit)
    output <- topTable(fit2, coef=ncol(design), number=nrow(v))
    output$gene_id <- rownames(output)
    output$gene_name <- gencode[match(output$gene_id, gencode$gene_id),'gene_name']
  }

  write.table(output, sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE, file=paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/differential_expression/', substr(currtiss, 1, nchar(currtiss)-10), '.de_by_sex.txt', sep=''))
}

library(parallel)
mclapply(filelist, worker, mc.cores=53)

#Some summary information
library(data.table)
library(rtracklayer)

gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode <- gencode[!gencode$seqid%in%c('MT'),]

filelist <- list.files('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/differential_expression')

resultframe <- NULL

for (currtiss in filelist) {
  outputs <- fread(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/differential_expression/', currtiss, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  outputs$tissue <- rep(substr(currtiss, 1, nchar(currtiss)-14), nrow(outputs))
  resultframe <- rbind(resultframe, outputs)
}

resultframe$CHR <- gencode[match(resultframe$gene_id, gencode$gene_id),'seqid']

resultframe$FDR_overall <- p.adjust(resultframe$P.Value, method='BH')

lncs <- subset(gencode, gene_type%in%c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA"))

resultframe_lncs <- subset(resultframe, gene_id%in%lncs$gene_id)
resultframe_lncs$FDR_lncs <- p.adjust(resultframe_lncs$P.Value, method='BH')

resultframe_pcs <- subset(resultframe, gene_id%in%subset(gencode, gene_type=='protein_coding')$gene_id)
resultframe_pcs$FDR_pcs <- p.adjust(resultframe_pcs$P.Value, method='BH')

library(ggplot2)

resultframe$gene_type <- 'Other'
resultframe[resultframe$gene_id%in%lncs$gene_id,'gene_type'] <- 'lncRNA'
resultframe[resultframe$gene_id%in%subset(gencode, gene_type=='protein_coding')$gene_id,'gene_type'] <- 'Protein coding'

prop.test(matrix(c(90, 1425-90, sum(resultframe$gene_type=='lncRNA'), nrow(resultframe)-sum(resultframe$gene_type=='lncRNA')), nrow=2, byrow=TRUE)) #DE lncRNAs vs total overall lncrnas
prop.test(matrix(c(88, 1, 1173, 1238-1173), nrow=2, byrow=TRUE)) #tissue-specific lncrna prop vs coding propr

pdf('/well/got2d/apayne/thesis/chapter4/gtex_de_results.pdf', width=5.9, height=5.9)
  ggplot(subset(resultframe, adj.P.Val<=0.05), aes(tissue)) +
    geom_bar(aes(fill=gene_type)) +
    scale_fill_manual(name='Gene type', values=c('darkorange', 'darkblue', 'darkgreen')) +
    labs(y='Count') +
    theme(axis.title.x=element_blank(),
    axis.text.x=element_text(angle=70, hjust=1),
    legend.position='top'
    )
dev.off()

des <- subset(resultframe, adj.P.Val<=0.05)
length(unique(des$gene_id))


head(resultframe_lncs[order(resultframe_lncs$P.Value, decreasing=FALSE),])

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

library(data.table)
library(rtracklayer)

gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode <- gencode[!gencode$seqid%in%c('MT'),]

filelist <- list.files('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/differential_expression')

resultframe <- NULL

for (currtiss in filelist) {
  outputs <- fread(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/differential_expression/', currtiss, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  outputs$tissue <- rep(substr(currtiss, 1, nchar(currtiss)-14), nrow(outputs))
  resultframe <- rbind(resultframe, outputs)
}

resultframe$CHR <- gencode[match(resultframe$gene_id, gencode$gene_id),'seqid']
resultframe$FDR_overall <- p.adjust(resultframe$P.Value, method='BH')

lncs <- subset(gencode, gene_type%in%c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA"))

resultframe_lncs <- subset(resultframe, gene_id%in%lncs$gene_id)
resultframe_lncs$FDR_lncs <- p.adjust(resultframe_lncs$P.Value, method='BH')

resultframe_pcs <- subset(resultframe, gene_id%in%subset(gencode, gene_type=='protein_coding')$gene_id)
resultframe_pcs$FDR_pcs <- p.adjust(resultframe_pcs$P.Value, method='BH')

lncs_sig <- subset(resultframe_lncs, adj.P.Val<=0.05)
pcs_sig <- subset(resultframe_pcs, adj.P.Val<=0.05)

lncs_un <- table(lncs_sig$gene_id)
pcs_un <- table(pcs_sig$gene_id)

lncs_un_tiss <- data.frame(gene_id=names(lncs_un), ntiss=c(lncs_un), stringsAsFactors=FALSE)
pcs_un_tiss <- data.frame(gene_id=names(pcs_un), ntiss=c(pcs_un), stringsAsFactors=FALSE)

plotframe <- rbind(lncs_un_tiss, pcs_un_tiss)
plotframe$gtype <- c(rep('lncRNA', nrow(lncs_un_tiss)), rep('Coding', nrow(pcs_un_tiss)))
plotframe$gtype <- factor(plotframe$gtype, levels=c('lncRNA', 'Coding'), labels=c('lncRNA', 'Coding'))

library(ggplot2)

p1 <- ggplot(plotframe) +
  geom_histogram(data=subset(plotframe, gtype=='lncRNA'),
    aes(y=..count../sum(..count..)*100, x=ntiss),
    bins=length(unique(subset(plotframe, gtype=='lncRNA')$ntiss)),
    fill='darkblue'
  ) +
  geom_text(data=subset(plotframe, gtype=='lncRNA'),
    stat='count',
    aes(x=ntiss,
      label=paste('n=', ..count.., sep=''),
      y=..count../sum(..count..)*100
    ),
    vjust=-0.4
  ) +
  scale_fill_discrete(drop=FALSE) +
  scale_x_continuous(breaks=sort(unique(subset(plotframe, gtype=='lncRNA')$ntiss)),
    labels=as.character(unique(subset(plotframe, gtype=='lncRNA')$ntiss))
  ) +
  scale_y_continuous(limits=c(0, 100)) +
  labs(x='Tissues with DE', y='Percentage of lncRNAs')

p2 <- ggplot(plotframe, aes(ntiss, fill=gtype)) +
  geom_histogram(data=subset(plotframe, gtype=='Coding'),
    aes(y=..count../sum(..count..)*100, x=ntiss),
    bins=length(unique(subset(plotframe, gtype=='Coding')$ntiss)),
  ) +
  geom_blank() +
  geom_text(data=subset(plotframe, gtype=='Coding'),
    stat='count',
    aes(x=ntiss,
      label=paste('n=', ..count.., sep=''),
      y=..count../sum(..count..)*100
    ),
    vjust=-0.4
  ) +
  labs(x='Tissues with DE', y='Percentage of coding genes') +
#  theme(axis.title.y=element_blank()) +
  scale_fill_manual(breaks=c('lncRNA', 'Coding'), values=c('darkorange', 'darkblue')) +
  scale_y_continuous(limits=c(0, 100)) +
  guides(fill=guide_legend(title='Gene type'))

pdf('/well/got2d/apayne/thesis/chapter4/gtex_de_results_tissue_sharing.pdf', width=5.9, height=3.5)
  multiplot(p1, p2, cols=2, colprops=c(0.45, 1))
dev.off()


table(plotframe$ntiss, plotframe$gtype)
1173/1238


#########################
# DE correlation with n #
#########################

filelist <- list.files('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/expr_covariates_voomed/')
filelist <- filelist[substring(filelist, nchar(filelist), nchar(filelist))=='m']

de_n <- NULL

for (i in filelist) {
  curr <- read.table(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/expr_covariates_voomed/', i, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, nrows=5)
  curr <- curr[,-1]
  de_n <- rbind(de_n, data.frame(tissue=i, n=ncol(curr), stringsAsFactors=FALSE))
}

de_n$tissue <- substring(de_n$tissue, 1, nchar(de_n$tissue)-4)

filelist <- list.files('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/differential_expression')

resultframe <- NULL

for (currtiss in filelist) {
  outputs <- fread(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/differential_expression/', currtiss, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  outputs$tissue <- rep(substr(currtiss, 1, nchar(currtiss)-14), nrow(outputs))
  resultframe <- rbind(resultframe, outputs)
}

resultframe_sig <- subset(resultframe, adj.P.Val<=0.05)
resultframe_tiss <- table(resultframe_sig$tissue)

de_n$sigde <- 0
for (i in names(resultframe_tiss)) de_n[match(i, de_n$tissue),'sigde'] <- resultframe_tiss[i]

cor(de_n[,2], de_n[,3])

library(ggplot2)

pdf('/well/got2d/apayne/thesis/chapter4/de_vs_n.pdf', width=5.9, height=3.5)
ggplot(de_n, aes(x=n, y=sigde)) +
  scale_y_continuous(limits=c(-30, max(de_n$sigde))) +
  geom_smooth(method='lm', formula=y~x, se=FALSE, col='red') +
  geom_point() +
  geom_text(data=subset(de_n, tissue%in%c('Muscle_Skeletal', 'Nerve_Tibial', 'Testis')),
    aes(x=n,
      y=sigde,
      label=tissue
    ),
    vjust=1.5,
    hjust=1.1
  ) +
  labs(x='Sample size',
    y='Number of differentially expressed genes'
  )
dev.off()






































#Check proportion of DE lncRNAs per tissue
filelist <- list.files('/well/got2d/apayne/GTEx_v7/gene_counts/')
filelist <- filelist[filelist!="Cells_Leukemia.reads.txt"]

summary_table <- NULL
sample_sizes <- NULL

summarizer <- function(currtiss) {
  expr <- fread(paste('/home/apayne/GTEx_v7/gene_counts/', currtiss, sep=''), sep='\t', stringsAsFactors=FALSE, data.table=FALSE)
  rownames(expr) <- expr[,1]
  expr <- expr[,-c(1, 2)]
  expr <- expr[rownames(expr)%in%gencode$gene_id,]
  expr <- as.matrix(expr)

  sampids <- c()

  for (j in 1:ncol(expr)) {
    sampids <- c(sampids, paste(strsplit(colnames(expr), '-', fixed=TRUE)[[j]][1], strsplit(colnames(expr), '-', fixed=TRUE)[[j]][2], sep='-'))
  }

  sex <- phenotype[match(sampids, phenotype$SUBJID),'SEX']
  t2d <- phenotype[match(sampids, phenotype$SUBJID),'MHT2D']
  keeps <- t2d!=99

  sex <- sex[keeps]
  t2d <- t2d[keeps]
  expr <- expr[,keeps]

  if (ncol(expr)<50) return(NULL)

  genereads <- expr[rowSums(expr>6)>=10,]
  gc()
  dge <- DGEList(counts=genereads)
  dge <- dge[rowSums(cpm(dge)>1, na.rm=TRUE)>=10,]

  if (length(unique(sex))!=2) {
    dge <- calcNormFactors(dge)
    v <- voom(dge, plot=F)$E
  } else {
    design <- model.matrix(~as.factor(sex))
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot=F)$E
  }

  currfull <- data.frame(Tissue=rep(currtiss, nrow(dge)), gene_id=rownames(dge), stringsAsFactors=FALSE)
  currfull$median_expr_full <- unlist(apply(v, 1, median))
  currfull$median_expr_t2d <- unlist(apply(v[,t2d==1], 1, median))
  currfull$median_expr_not2d <- unlist(apply(v[,t2d==2], 1, median))

  currsamp <- data.frame(n=length(t2d), n.t2d=sum(t2d==1), n.not2d=sum(t2d==2))
  return(list(currfull, currsamp))
}

library(parallel)
templist <- mclapply(filelist, FUN=summarizer, mc.cores=27)

for (i in 1:length(templist)) {
  summary_table <- rbind(summary_table, templist[[i]][[1]])
  sample_sizes <- rbind(sample_sizes, templist[[i]][[2]])
}

rownames(sample_sizes) <- substring(filelist, 1, nchar(filelist)-10)

cor.test(table(subset(resultframe_lncs, FDR_overall<=0.01)$tissue), sample_sizes$n)




summary_table_lnc <- summary_table[summary_table$gene_id%in%lncs$gene_id,]

total_lncs <- c(table(substring(summary_table_lnc$Tissue, 1, nchar(summary_table_lnc$Tissue)-10)))
total_de_lncs <- c(table(subset(resultframe_lncs, FDR_lncs<=0.01)$tissue))
de_props_lncs <- total_de_lncs/total_lncs
de_props_df <- data.frame(Tissue=names(de_props_lncs), stringsAsFactors=FALSE)
de_props_df$lncs_expressed <- total_lncs
de_props_df$lncs_de <- total_de_lncs
de_props_df$de_prop <- de_props_lncs

cor.test(de_props_df$de_prop, sample_sizes$n)


summary_table_pc <- summary_table[summary_table$gene_id%in%subset(gencode, gene_type=='protein_coding')$gene_id,]

total_pcs <- c(table(substring(summary_table_pc$Tissue, 1, nchar(summary_table_pc$Tissue)-10)))
total_de_pcs <- c(table(subset(resultframe_pcs, FDR_pcs<=0.01)$tissue))
de_props_pcs <- total_de_pcs/total_pcs
de_props_pc <- data.frame(Tissue=names(de_props_pcs), stringsAsFactors=FALSE)
de_props_pc$pcs_expressed <- total_pcs
de_props_pc$pcs_de <- total_de_pcs
de_props_pc$de_prop <- de_props_pcs

cor.test(de_props_pc$de_prop, sample_sizes$n)


library(ggplot2)
library(devEMF)
library(gridExtra)

top_de_props_lncs <- de_props_df[order(de_props_df$de_prop, decreasing=TRUE)[2:11],]
top_de_props_pcs <- de_props_pc[order(de_props_pc$de_prop, decreasing=TRUE)[2:11],]

svg('/well/got2d/apayne/GTEx_v7/sex_specific_analyses/plots/de_props_top10_pc_lnc.svg', width=16, height=9)
p1 <- qplot(factor(as.factor(top_de_props_lncs$Tissue), levels=top_de_props_lncs$Tissue), weight=top_de_props_lncs$de_prop, geom='bar', ylab='Proportion of lncRNAs differentially expressed', xlab='Tissue') + labs(title='Tissues with highest % genes differentially expressed') + theme(axis.text.x=element_text(angle=90, hjust=1, size=12), axis.text.y=element_text(size=12), axis.title=element_text(size=14, face='bold'), title=element_text(size=16, face='bold')) + geom_bar(fill='red')
p2 <- qplot(factor(as.factor(top_de_props_pcs$Tissue), levels=top_de_props_pcs$Tissue), weight=top_de_props_pcs$de_prop, geom='bar', ylab='Proportion of coding genes differentially expressed', xlab='Tissue') + labs(title='') + theme(axis.text.x=element_text(angle=90, hjust=1, size=12), axis.text.y=element_text(size=12), axis.title=element_text(size=14, face='bold'), title=element_text(size=16, face='bold')) + geom_bar(fill='blue')
grid.arrange(p1, p2, ncol=2)
dev.off()





lncprops_top <- de_props_df[order(de_props_df$lncs_expressed, decreasing=TRUE),][1:15,]
cor (de_props_df[,2], sample_sizes[,1])



summary_table_lnc$tiss_num <- unlist(mclapply(summary_table_lnc$gene_id, FUN=function(x) sum(summary_table_lnc$gene_id==x), mc.cores=16))
summary_table_lnc$tissrange <- NA
summary_table_lnc[summary_table_lnc$tiss_num==1,'tissrange'] <- 'Unique'
summary_table_lnc[summary_table_lnc$tiss_num>1&summary_table_lnc$tiss_num<=4,'tissrange'] <- "2-4 tissues"
summary_table_lnc[summary_table_lnc$tiss_num>4&summary_table_lnc$tiss_num<=10,'tissrange'] <- "5-10 tissues"
summary_table_lnc[summary_table_lnc$tiss_num>10,'tissrange'] <- "More than 10 tissues"
summary_table_lnc$tissrange <- factor(as.factor(summary_table_lnc$tissrange), levels=c('Unique', '2-4 tissues', '5-10 tissues', 'More than 10 tissues'))

toptissues <- names(sort(table(summary_table_lnc$Tissue), decreasing=TRUE))[1:15]
summary_table_lnc_tops <- subset(summary_table_lnc, Tissue%in%toptissues)
summary_table_lnc_tops$Tissue <- substring(summary_table_lnc_tops$Tissue, 1, nchar(summary_table_lnc_tops$Tissue)-10)
toptissues <- substring(toptissues, 1, nchar(toptissues)-10)

svg('/well/got2d/apayne/GTEx_v7/sex_specific_analyses/plots/top_tissues_expressed_lncs.svg', width=16, height=9)
qplot(factor(as.factor(Tissue), levels=toptissues), data=summary_table_lnc_tops, fill=tissrange, geom='bar', ylab='Number of expressed lncRNAs', xlab='Tissue') + labs(title='Number of expressed lncRNAs in top tissues') + theme(axis.text.x=element_text(angle=45, hjust=1, size=12), axis.text.y=element_text(size=12), axis.title=element_text(size=14, face='bold'), title=element_text(size=16, face='bold'), legend.text=element_text(size=12))
dev.off()




topdelncs <- sort(table(subset(resultframe_lncs, FDR_overall<=0.01)$gene_name), decreasing=TRUE)[1:11]














