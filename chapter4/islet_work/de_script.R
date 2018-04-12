library(sva)
library(edgeR)
library(limma)
library(data.table)
library(WGCNA)
library(rtracklayer)

gencode <- readGFF('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz')
gencode <- gencode[!gencode$seqid%in%c('MT'),]

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

design <- model.matrix(~ svars + as.factor(sex) + as.factor(t2d))

colnames(design) <- c('Intercept', paste('SVA', 1:ncol(svars), sep=''), 'sex', 't2d')
fit <- lmFit(v, design)
fit2 <- eBayes(fit)
output <- topTable(fit2, coef=ncol(design), number=nrow(v))
output$gene_id <- rownames(output)
output$gene_name <- gencode[match(output$gene_id, gencode$gene_id),'gene_name']

write.table(output, sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE, file='/well/got2d/apayne/islet_t2d_networks/de_by_t2d.txt')


#Some summary information
library(data.table)
library(rtracklayer)

gencode <- readGFF('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz')
gencode <- gencode[!gencode$seqid%in%c('MT'),]

resultframe <- fread('/well/got2d/apayne/islet_t2d_networks/de_by_t2d.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

resultframe$CHR <- gencode[match(resultframe$gene_id, gencode$gene_id),'seqid']

lncs <- subset(gencode, gene_type%in%c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA"))

resultframe_lncs <- subset(resultframe, gene_id%in%lncs$gene_id)
resultframe_lncs$FDR_lncs <- p.adjust(resultframe_lncs$P.Value, method='BH')

resultframe_pcs <- subset(resultframe, gene_id%in%subset(gencode, gene_type=='protein_coding')$gene_id)
resultframe_pcs$FDR_pcs <- p.adjust(resultframe_pcs$P.Value, method='BH')

resultframe$gene_type_simple <- 'Other'
resultframe[resultframe$gene_id%in%lncs$gene_id,'gene_type_simple'] <- 'lncRNA'
resultframe[resultframe$gene_id%in%subset(gencode, gene_type=='protein_coding')$gene_id,'gene_type_simple'] <- 'Protein coding'
resultframe$gene_type <- gencode[match(resultframe$gene_id, gencode$gene_id),'gene_type']

des <- subset(resultframe, adj.P.Val<=0.05)
dim(des)
table(des$gene_type_simple)
table(des$gene_type)

prop.test(matrix(c(75, nrow(des)-75, sum(resultframe$gene_type_simple=='lncRNA'), nrow(resultframe)-sum(resultframe$gene_type_simple=='lncRNA')), nrow=2, byrow=2))

prop.test(matrix(c(1652, nrow(des)-1652, sum(resultframe$gene_type_simple=='Protein coding'), nrow(resultframe)-sum(resultframe$gene_type_simple=='Protein coding')), nrow=2, byrow=2))

