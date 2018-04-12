library(edgeR)
library(limma)
library(sva)
library(rtracklayer)
library(data.table)

#
gencode <- readGFF('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz')
gencode_genes <- subset(gencode, type=='gene')
gencode_transcripts <- subset(gencode, type=='transcript')

ctbind <- match('CTB-12O2.1', gencode_genes$gene_name)
gencode_genes[seq(ctbind-4, ctbind+4, 1),1:15]
subset(gencode_transcripts, gene_name=='CTB-12O2.1')
#
currtiss <- 'Pancreas.reads.txt'
gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode$seqid <- as.character(gencode$seqid)
gencode <- gencode[gencode$seqid!='MT',]
phenotype <- read.table('/well/got2d/GTEx_v7/sample_annotations/GTEx_Analysis_2016-01-15_v7_SubjectPhenotypesDS.txt', header=TRUE, stringsAsFactors=FALSE, quote="\"", sep="\t", comment.char="")
expr <- fread('/well/got2d/apayne/GTEx_v7/gene_counts/Pancreas.reads.txt', sep='\t', stringsAsFactors=FALSE, data.table=FALSE)
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
sampids <- sampids[keeps]

genereads <- expr[rowSums(expr>6)>=10,]
gc()
dge <- DGEList(counts=genereads)
dge <- dge[rowSums(cpm(dge)>1, na.rm=TRUE)>=10,]
dge$counts <- dge$counts + 1

design <- model.matrix(~ as.factor(sex) + as.factor(t2d))
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=F)$E

svars <- sva(dat=v, mod=design)$sv
colnames(svars) <- paste('SVA', 1:ncol(svars), sep='')

pancreas_newvals <- matrix(0, nrow=nrow(v), ncol=ncol(v))
rownames(pancreas_newvals) <- rownames(v)
colnames(pancreas_newvals) <- colnames(v)

if (length(unique(sex))!=1) {
  for (k in 1:nrow(v)) pancreas_newvals[k,] <- resid(lm(v[k,] ~ svars + sex))
} else {
  for (k in 1:nrow(v)) pancreas_newvals[k,] <- resid(lm(v[k,] ~ svars))
}



igencode <- readGFF('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz')
igencode <- igencode[!igencode$seqid%in%c('MT'),]

iphenotype <- read.table('/well/got2d/mvdbunt/InsPIRE/current_analysis/Cov_AllFilesIslets.csv', header=TRUE, stringsAsFactors=FALSE, sep=",")

iexpr <- fread('/well/got2d/apayne/InsPIRE_data/islet.overlapsremoved.gene.reads.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

rownames(iexpr) <- iexpr[,1]
iexpr <- iexpr[,-1]
iexpr <- as.matrix(iexpr)

isex <- iphenotype[match(colnames(iexpr), iphenotype$NewSampleID),'Gender']
it2d <- iphenotype[match(colnames(iexpr), iphenotype$NewSampleID),'T2D']
ikeeps <- !is.na(it2d)

isex <- isex[ikeeps]
it2d <- it2d[ikeeps]
iexpr <- iexpr[,ikeeps]

igenereads <- iexpr[rowSums(iexpr>6)>=10,]
gc()
idge <- DGEList(counts=igenereads)
idge <- idge[rowSums(cpm(idge)>1, na.rm=TRUE)>=10,]
idge$counts <- idge$counts + 1

idesign <- model.matrix(~ as.factor(isex) + as.factor(it2d))
idge <- calcNormFactors(idge)
iv <- voom(idge, idesign, plot=F)$E

isvars <- sva(dat=iv, mod=idesign)$sv

colnames(isvars) <- paste('SVA', 1:ncol(isvars), sep='')

inewvals <- matrix(0, nrow=nrow(iv), ncol=ncol(iv))
rownames(inewvals) <- rownames(iv)
colnames(inewvals) <- colnames(iv)

for (k in 1:nrow(iv)) inewvals[k,] <- resid(lm(iv[k,] ~ isvars + isex))


t2d_loci <- read.table('/well/got2d/apayne/thesis/t2d_loci/t2d_loci_curated_from_anubha.txt', sep='\t', header=FALSE, stringsAsFactors=FALSE)

colnames(t2d_loci) <- 'gene_name'
t2d_loci$gene_id <- gencode[match(t2d_loci$gene_name, gencode$gene_name),'gene_id']
t2d_loci <- na.omit(t2d_loci)

t2d_loci_in_both <- t2d_loci[t2d_loci$gene_id%in%rownames(pancreas_newvals)&t2d_loci$gene_id%in%rownames(inewvals),]


ctbid <- gencode[match('CTB-12O2.1', gencode$gene_name),'gene_id']
library(ggplot2)

pcormat <- data.frame(gene_name=t2d_loci_in_both$gene_name, stringsAsFactors=FALSE)
pcormat$Correlation <- unlist(lapply(t2d_loci_in_both$gene_id, FUN=function(x) cor(pancreas_newvals[ctbid,], pancreas_newvals[x,])))
pcormat <- pcormat[abs(pcormat$Correlation)>=0.05,]
pcormat <- rbind(pcormat, data.frame(gene_name=c('·   ', ' ·   ', '  ·   '), Correlation=rep(0, 3), stringsAsFactors=FALSE))
pcormat$gene_name <- factor(pcormat$gene_name, levels=pcormat[order(pcormat$Correlation, decreasing=FALSE),'gene_name'])

pdf('/well/got2d/apayne/thesis/chapter4/pancreas_ctb_cor_t2d.pdf', width=5.9, height=9.3)
  ggplot(aes(x=gene_name, y=Correlation), data=pcormat) +
    geom_bar(stat='identity', fill='darkblue') +
    theme(axis.text.y=element_text(size=8)) +
    coord_flip()
dev.off()

icormat <- data.frame(gene_name=t2d_loci_in_both$gene_name, stringsAsFactors=FALSE)
icormat$Correlation <- unlist(lapply(t2d_loci_in_both$gene_id, FUN=function(x) cor(inewvals[ctbid,], inewvals[x,])))
icormat <- icormat[abs(icormat$Correlation)>=0.05,]
icormat <- rbind(icormat, data.frame(gene_name=c('·   ', ' ·   ', '  ·   '), Correlation=rep(0, 3), stringsAsFactors=FALSE))
icormat$gene_name <- factor(icormat$gene_name, levels=icormat[order(icormat$Correlation, decreasing=FALSE),'gene_name'])

pdf('/well/got2d/apayne/thesis/chapter4/islet_ctb_cor_t2d.pdf', width=5.9, height=9.3)
  ggplot(aes(x=gene_name, y=Correlation), data=icormat) +
    geom_bar(stat='identity', fill='darkblue') +
    theme(axis.text.y=element_text(size=8)) +
    coord_flip()
dev.off()


##########################
# EXPRESSION ACROSS GTEX #
##########################
library(sva)
library(edgeR)
library(limma)
library(data.table)
library(WGCNA)
library(rtracklayer)

filelist <- list.files('/well/got2d/apayne/GTEx_v7/gene_counts/')
gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode <- gencode[gencode$seqid!='MT',]

phenotype <- read.table('/well/got2d/GTEx_v7/sample_annotations/GTEx_Analysis_2016-01-15_v7_SubjectPhenotypesDS.txt', header=TRUE, stringsAsFactors=FALSE, quote="\"", sep="\t", comment.char="")

filelist <- filelist[filelist!="Cells_Leukemia.reads.txt"]

exfilter <- function(currtiss) {
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
  t2d <- phenotype[match(sampids, phenotype$SUBJID),'MHT2D']
  keeps <- t2d!=99

  sex <- sex[keeps]
  t2d <- t2d[keeps]
  expr <- expr[,keeps]

  if (ncol(expr)<50) return(NULL)

  dge1 <- DGEList(counts=expr)
  cpmvals1 <- cbind(data.frame(gene_id=rownames(cpm(dge1)), stringsAsFactors=FALSE), as.data.frame(cpm(dge1)))
  write.table(cpmvals1, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/expr_covariates_voomed/', substring(currtiss, 1, nchar(currtiss)-10), '.cpm.unfiltered', sep=''))

  genereads <- expr[rowSums(expr>6)>=10,]
  gc()
  dge <- DGEList(counts=genereads)
  dge <- dge[rowSums(cpm(dge)>1, na.rm=TRUE)>=10,]
  cpmvals <- cbind(data.frame(gene_id=rownames(cpm(dge)), stringsAsFactors=FALSE), as.data.frame(cpm(dge)))
  write.table(cpmvals, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/expr_covariates_voomed/', substring(currtiss, 1, nchar(currtiss)-10), '.cpm', sep=''))

  dge$counts <- dge$counts + 1

  if (length(unique(sex))!=1) {
    design <- model.matrix(~ as.factor(sex) + as.factor(t2d))
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot=F)$E

    svars <- sva(dat=v, mod=design)$sv
    colnames(svars) <- paste('SVA', 1:ncol(svars), sep='')

  } else {
    design <- model.matrix(~ as.factor(t2d))
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot=F)$E

    svars <- sva(dat=v, mod=design)$sv
    colnames(svars) <- paste('SVA', 1:ncol(svars), sep='')
  }

  newvals <- matrix(0, nrow=nrow(v), ncol=ncol(v))
  rownames(newvals) <- rownames(v)
  colnames(newvals) <- colnames(v)

  if (length(unique(sex))!=1) {
    for (k in 1:nrow(v)) newvals[k,] <- resid(lm(v[k,] ~ svars + sex))
  } else {
    for (k in 1:nrow(v)) newvals[k,] <- resid(lm(v[k,] ~ svars))
  }

  writevals <- cbind(data.frame(gene_id=rownames(newvals), stringsAsFactors=FALSE), as.data.frame(newvals))
  write.table(writevals, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/expr_covariates_voomed/', substring(currtiss, 1, nchar(currtiss)-10), '.newvals', sep=''))
}

library(parallel)
mclapply(filelist, FUN=exfilter, mc.cores=27)


phenotype <- read.table('/well/got2d/mvdbunt/InsPIRE/current_analysis/Cov_AllFilesIslets.csv', header=TRUE, stringsAsFactors=FALSE, sep=",")

expr <- fread('/well/got2d/apayne/InsPIRE_data/islet.overlapsremoved.gene.reads.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.matrix(expr)

sex <- phenotype[match(colnames(expr), phenotype$NewSampleID),'Gender']
t2d <- phenotype[match(colnames(expr), phenotype$NewSampleID),'T2D']
keeps <- !is.na(t2d)

sex <- sex[keeps]
t2d <- t2d[keeps]
expr <- expr[,keeps]

dge1 <- DGEList(counts=expr)
cpmvals1 <- cbind(data.frame(gene_id=rownames(cpm(dge1)), stringsAsFactors=FALSE), as.data.frame(cpm(dge1)))
write.table(cpmvals1, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/islet_t2d_networks/islet.cpm.unfiltered')


genereads <- expr[rowSums(expr>6)>=10,]
gc()
dge <- DGEList(counts=genereads)
dge <- dge[rowSums(cpm(dge)>1, na.rm=TRUE)>=10,]
cpmvals <- cbind(data.frame(gene_id=rownames(cpm(dge)), stringsAsFactors=FALSE), as.data.frame(cpm(dge)))
write.table(cpmvals, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/islet_t2d_networks/islet.cpm')

dge$counts <- dge$counts + 1

design <- model.matrix(~ as.factor(sex) + as.factor(t2d))
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=F)$E

svars <- sva(dat=v, mod=design)$sv

colnames(svars) <- paste('SVA', 1:ncol(svars), sep='')

newvals <- matrix(0, nrow=nrow(v), ncol=ncol(v))
rownames(newvals) <- rownames(v)
colnames(newvals) <- colnames(v)

for (k in 1:nrow(v)) newvals[k,] <- resid(lm(v[k,] ~ svars + sex))

writevals <- cbind(data.frame(gene_id=rownames(newvals), stringsAsFactors=FALSE), as.data.frame(newvals))
write.table(writevals, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/islet_t2d_networks/islet.newvals')


##############################
# OVERALL EXPRESSION PROFILE #
##############################
library(data.table)
library(ggplot2)

filelist <- list.files('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/expr_covariates_voomed')
filelist <- filelist[substring(filelist, nchar(filelist), nchar(filelist))=='d']
ctbid <- "ENSG00000254226.1"

ctbexpr <- list()

for (i in filelist) {
  curr <- fread(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/expr_covariates_voomed/', i, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  rownames(curr) <- curr[,1]
  curr <- curr[,-1]
  curr <- as.matrix(curr)
  ctbexpr[[i]] <- curr[rownames(curr)==ctbid,]
  print(i)
}

iexpr <- read.table('/well/got2d/apayne/islet_t2d_networks/islet.cpm.unfiltered', sep='\t', header=TRUE, stringsAsFactors=FALSE)
rownames(iexpr) <- iexpr[,1]
iexpr <- iexpr[,-1]
iexpr <- as.matrix(iexpr)

ctbexpr[['Islet.cpm.unfiltered']] <- iexpr[rownames(iexpr)==ctbid,]
rm(iexpr, i, curr)
gc()

plotvals <- NULL

for (i in names(ctbexpr)) {
  plotvals <- rbind(plotvals, data.frame(Tissue=rep(substring(i, 1, nchar(i)-15), length(ctbexpr[[i]])), cpm=ctbexpr[[i]], stringsAsFactors=FALSE))
}

plotvals$log2cpmp1 <- log2(plotvals$cpm+1)

plotvals$color <- rep('darkblue', nrow(plotvals))
plotvals[plotvals$Tissue%in%c('Pancreas', 'Islet'),'color'] <- 'darkorange'

pdf('/well/got2d/apayne/thesis/chapter4/ctb_expr.pdf', width=5.9, height=7)
ggplot(data=plotvals, aes(x=as.factor(Tissue), y=log2cpmp1, fill=color, colour=color)) +
  geom_violin() +
  scale_fill_manual(values=c('darkblue', 'darkorange')) +
  scale_colour_manual(values=c('darkblue', 'darkorange')) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=1, colour='red') +
  scale_y_continuous(breaks=c(0, 1, 2, 3, 4, 5), 
    labels=c('0', '1', '3', '7', '15', '31')
  ) +
  theme(axis.title.y=element_blank()) +
  labs(y='Counts per million reads') +
  guides(fill=FALSE, colour=FALSE) +
  coord_flip()
dev.off()


###############################
# GLRA1, GCG, INS correlation #
###############################

library(data.table)
library(ggplot2)

filelist <- list.files('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/expr_covariates_voomed')
filelist <- filelist[substring(filelist, nchar(filelist), nchar(filelist))=='s']
ctbid <- "ENSG00000254226.1"
glrid <- "ENSG00000145888.6"
gcgid <- 'ENSG00000115263.10'
insid <- 'ENSG00000254647.2'
expr <- list()

for (i in filelist) {
  curr <- fread(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/expr_covariates_voomed/', i, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  rownames(curr) <- curr[,1]
  curr <- curr[,-1]
  curr <- as.matrix(curr)
  expr[[i]] <- curr
  print(i)
}

iexpr <- read.table('/well/got2d/apayne/islet_t2d_networks/islet.newvals', sep='\t', header=TRUE, stringsAsFactors=FALSE)
rownames(iexpr) <- iexpr[,1]
iexpr <- iexpr[,-1]
iexpr <- as.matrix(iexpr)

expr[['Islet.newvals']] <- iexpr

glrc <- c()
gcgc <- c()
insc <- c()

for (i in 1:length(expr)) {
  if(ctbid%in%rownames(expr[[i]])&glrid%in%rownames(expr[[i]])) {
    glrc <- c(glrc, cor(expr[[i]][ctbid,], expr[[i]][glrid,]))
    names(glrc) <- c(names(glrc)[1:length(glrc)-1], names(expr)[i])
  }

  if(ctbid%in%rownames(expr[[i]])&gcgid%in%rownames(expr[[i]])) {
    gcgc <- c(gcgc, cor(expr[[i]][ctbid,], expr[[i]][gcgid,]))
    names(gcgc) <- c(names(gcgc)[1:length(gcgc)-1], names(expr)[i])
  }

  if(ctbid%in%rownames(expr[[i]])&insid%in%rownames(expr[[i]])) {
    insc <- c(insc, cor(expr[[i]][ctbid,], expr[[i]][insid,]))
    names(insc) <- c(names(insc)[1:length(insc)-1], names(expr)[i])
  }
}



#############################
# Number of tissues with DE #
#############################
library(data.table)

filelist <- list.files('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/differential_expression')

resultframe <- NULL

for (currtiss in filelist) {
  outputs <- fread(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/differential_expression/', currtiss, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  outputs$tissue <- rep(substr(currtiss, 1, nchar(currtiss)-14), nrow(outputs))
  resultframe <- rbind(resultframe, outputs)
}

ide <- fread('/well/got2d/apayne/islet_t2d_networks/de_by_t2d.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
ide$tissue <- 'Islet'
resultframe <- rbind(resultframe, ide)

ctbde <- subset(resultframe, gene_name=='CTB-12O2.1')
ctbde$tissue <- factor(ctbde$tissue, ctbde[order(ctbde$P.Value, decreasing=FALSE),'tissue'])

library(ggplot2)

pdf('/well/got2d/apayne/thesis/chapter4/ctb_de_tissue.pdf', width=5.9, height=4)
ggplot(ctbde, aes(x=tissue, y=-log10(P.Value))) +
  geom_point() +
  theme(axis.title.x=element_blank(),
    axis.text.x=element_text(angle=70, hjust=1)
  ) +
#  labs(y=paste('-', expression(log[10]), '(T2D DE p-value)', sep=''))
  labs(y=expression(-log[10]*"(T2D DE p-value)"))
dev.off()


##########################
# Overlap islet and GTEx #
#########################
library(data.table)

filelist <- list.files('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/differential_expression')

resultframe <- NULL

for (currtiss in filelist) {
  outputs <- fread(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/differential_expression/', currtiss, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  outputs$tissue <- rep(substr(currtiss, 1, nchar(currtiss)-14), nrow(outputs))
  resultframe <- rbind(resultframe, outputs)
}

ide <- fread('/well/got2d/apayne/islet_t2d_networks/de_by_t2d.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
ide$tissue <- 'Islet'
resultframe <- rbind(resultframe, ide)
sigs <- subset(resultframe, adj.P.Val<=0.05)

sigsi <- subset(sigs, tissue=='Islet')
sigsg <- subset(sigs, tissue!='Islet')
sigsg_over <- subset(sigsg, gene_id%in%sigsi$gene_id)
length(unique(sigsg_over$gene_id))


###############################
# DE for only T2D or non-T2D? #
###############################








