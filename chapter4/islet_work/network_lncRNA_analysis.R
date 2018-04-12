library(data.table)
library(rtracklayer)
library(ggplot2)

gencode <- readGFF('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz')
gencode <- gencode[!gencode$seqid%in%c('MT'),]

lncs <- subset(gencode, gene_type%in%c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA"))

phenotype <- read.table('/well/got2d/mvdbunt/InsPIRE/current_analysis/Cov_AllFilesIslets.csv', header=TRUE, stringsAsFactors=FALSE, sep=",")

modules <- fread('/well/got2d/apayne/islet_t2d_networks/modules.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
uniquemods <- modules[match(unique(modules$color), modules$color),]
nrow(uniquemods)-1

modsums <- data.frame(genes=nrow(modules), modules=length(unique(modules$color))-1, stringsAsFactors=FALSE)

de <- fread('/well/got2d/apayne/islet_t2d_networks/de_by_t2d.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

de <- de[de$gene_id%in%gencode$gene_id,]
de$FDR_overall <- p.adjust(de$P.Value, method='BH')
de$sig <- de$FDR_overall<=0.01
de$gene_type <- gencode[match(de$gene_id, gencode$gene_id),'gene_type']
de$simple_gene_type <- rep(NA, nrow(de))
de[de$gene_type%in%lncs$gene_type,'simple_gene_type'] <- 'lncRNA'
de[de$gene_type%in%'protein_coding','simple_gene_type'] <- 'protein_coding'
de[is.na(de$simple_gene_type),'simple_gene_type'] <- 'other'

mod_assigned <- subset(modules, module_number!=0)
de_assigned <- subset(de, gene_id%in%mod_assigned$gene_id)

if (sum(de_assigned$P.Value==0)>0) de_assigned[de_assigned$P.Value==0,'P.Value'] <- rep(min(subset(de_assigned, P.Value!=0)$P.Value), sum(de_assigned$P.Value==0))

de_assigned_lnc <- subset(de_assigned, simple_gene_type=='lncRNA')

module_list <- list()

for (i in 1:length(unique(mod_assigned$module_number))) {
  module_list[[subset(mod_assigned, module_number==i)$color[1]]] <- list()
  module_list[[i]][[1]] <- de_assigned[match(subset(mod_assigned, module_number==i)$gene_id, de_assigned$gene_id),c('gene_id', 'gene_name', 'simple_gene_type', 'P.Value', 'FDR_overall')]

  if (nrow(subset(module_list[[i]][[1]], simple_gene_type=='lncRNA'))==0) {
    module_list[[i]][[2]] <- subset(module_list[[i]][[1]], simple_gene_type=='lncRNA')
    module_list[[i]][[3]] <- 0
    module_list[[i]][[4]] <- rep(0, 10000)
    module_list[[i]][[5]] <- 1
    names(module_list[[i]]) <- c('full_module', 'full_lncs', 'fisher_score', 'random_scores', 'module_p')
    next
  }

  module_list[[i]][[2]] <- subset(module_list[[i]][[1]], simple_gene_type=='lncRNA')
  module_list[[i]][[3]] <- -2*sum(log(module_list[[i]][[2]]$P.Value))
  siglncs <- subset(module_list[[i]][[2]], FDR_overall<=0.01)
  module_list[[i]][[4]] <- unlist(lapply(1:10000, FUN=function(x) -2*sum(log(de_assigned_lnc[sample(1:nrow(de_assigned_lnc), nrow(module_list[[i]][[2]]), replace=FALSE),'P.Value']))))
  module_list[[i]][[5]] <- sum(module_list[[i]][[4]]>module_list[[i]][[3]])/10000
  names(module_list[[i]]) <- c('full_module', 'full_lncs', 'fisher_score', 'random_scores', 'module_p')
}

modframe <- NULL

for (i in 1:length(module_list)) modframe <- rbind(modframe, data.frame(module=i, color=modules[match(i, modules$module_number),'color'], pval=module_list[[i]][[5]], stringsAsFactors=FALSE))

modframe$p.adj <- p.adjust(modframe$pval, method='BH')
subset(modframe, p.adj<=0.05)

significant_modules <- NULL

for (i in subset(modframe, p.adj<0.05)$module) {
  significant_modules <- rbind(significant_modules, cbind(module_list[[i]][[1]], data.frame(module=rep(i, nrow(module_list[[i]][[1]])), color=rep(modules[match(i, modules$module_number),'color'], nrow(module_list[[i]][[1]])))))
}

write.table(significant_modules, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE, file='/well/got2d/apayne/thesis/chapter4/islet_work/significant_modules_05_bh/islet.lncRNA_de_modules.txt')


#############################################################
# Overlap of significant module with all other GTEx modules #
#############################################################
islet_mod <- read.table('/well/got2d/apayne/thesis/chapter4/islet_work/significant_modules_05_bh/islet.lncRNA_de_modules.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

tisslist <- list.files('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/coexpression_networks')
tisslist <- tisslist[substring(tisslist, nchar(tisslist)-4, nchar(tisslist)-4)=='s']

mod_overlaps <- NULL

for (i in tisslist) {
  currmods <- read.table(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/coexpression_networks/', i, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  currmods <- currmods[currmods$color!='grey',]

  for (j in unique(currmods$color)) {
    currcol <- subset(currmods, color==j)
    mod_overlaps <- rbind(mod_overlaps, data.frame(tissue=i, module=j, n_genes_mod=nrow(currcol), n_genes_isl=nrow(islet_mod), n_overlap=sum(currcol$gene_id%in%islet_mod$gene_id), stringsAsFactors=FALSE))
  }

  print(i)
}

mod_overlaps$p_of_mod <- mod_overlaps$n_overlap/mod_overlaps$n_genes_mod
mod_overlaps$p_of_isl <- mod_overlaps$n_overlap/mod_overlaps$n_genes_isl

head(mod_overlaps[order(mod_overlaps$p_of_mod, decreasing=TRUE),])
head(mod_overlaps[order(mod_overlaps$p_of_isl, decreasing=TRUE),])
prop.test(matrix(c(25, 54-25, 3, 23-3), nrow=2, byrow=TRUE))

###################################
# COMPARISON WITH PANCREAS MODULE #
###################################
library(edgeR)
library(limma)
library(sva)
library(rtracklayer)
library(data.table)

gencode <- readGFF('/well/got2d/rna-seq/resources/gencode.v19.annotation.gtf')

islet_mod <- read.table('/well/got2d/apayne/thesis/chapter4/islet_work/significant_modules_05_bh/islet.lncRNA_de_modules.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
pancreas_mods <- read.table('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/coexpression_networks/Pancreas.modules.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
pancreas_mods <- subset(pancreas_mods, module_number!=0)

pancreas_overlaps <- matrix(0, nrow=length(unique(pancreas_mods$color)), ncol=3)
colnames(pancreas_overlaps) <- c('n_genes', 'n_overlap', 'prop_overlap')
rownames(pancreas_overlaps) <- unique(pancreas_mods$color)

for (i in rownames(pancreas_overlaps)) {
  pancreas_overlaps[i,1] <- nrow(subset(pancreas_mods, color==i))
  pancreas_overlaps[i,2] <- nrow(subset(pancreas_mods, color==i&gene_id%in%islet_mod$gene_id))
  pancreas_overlaps[i,3] <- pancreas_overlaps[i,2]/pancreas_overlaps[i,1]
}

pancreas_overlaps <- pancreas_overlaps[order(pancreas_overlaps[,3], decreasing=TRUE),]

pancreas_sig <- subset(pancreas_mods, color==rownames(pancreas_overlaps)[1])
pancreas_sig_overlap <- subset(pancreas_sig, gene_id%in%islet_mod$gene_id)

islet_de <- read.table('/well/got2d/apayne/islet_t2d_networks/de_by_t2d.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
islet_de <- subset(islet_de, gene_id%in%islet_mod$gene_id)
dim(subset(islet_de, adj.P.Val<=0.05))
dim(subset(islet_de, adj.P.Val<=0.05&gene_id%in%pancreas_sig_overlap$gene_id))

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


overlap_sums <- data.frame(Gene=unlist(lapply(pancreas_sig_overlap$gene_name, FUN=function(x) rep(x, 4))), T2D=rep(c('T2D', 'No T2D'), 2*nrow(pancreas_sig_overlap)), Tissue=rep(c(rep('Pancreas', 2), rep('Islet', 2)), nrow(pancreas_sig_overlap)), stringsAsFactors=FALSE)

overlap_sums$mean_expr <- 0

for (i in unique(overlap_sums$Gene)) {
  overlap_sums[overlap_sums$Gene==i&overlap_sums$T2D=='T2D'&overlap_sums$Tissue=='Pancreas','mean_expr'] <- mean(pancreas_newvals[pancreas_sig_overlap[match(i, pancreas_sig_overlap$gene_name),'gene_id'],t2d==1])
  overlap_sums[overlap_sums$Gene==i&overlap_sums$T2D=='No T2D'&overlap_sums$Tissue=='Pancreas','mean_expr'] <- mean(pancreas_newvals[pancreas_sig_overlap[match(i, pancreas_sig_overlap$gene_name),'gene_id'],t2d==0])
  overlap_sums[overlap_sums$Gene==i&overlap_sums$T2D=='T2D'&overlap_sums$Tissue=='Islet','mean_expr'] <- mean(inewvals[pancreas_sig_overlap[match(i, pancreas_sig_overlap$gene_name),'gene_id'],it2d=='Yes'])
  overlap_sums[overlap_sums$Gene==i&overlap_sums$T2D=='No T2D'&overlap_sums$Tissue=='Islet','mean_expr'] <- mean(inewvals[pancreas_sig_overlap[match(i, pancreas_sig_overlap$gene_name),'gene_id'],it2d=='No'])
}

overlap_sums$X <- factor(unlist(lapply(1:(nrow(overlap_sums)/2), FUN=function(x) rep(x, 2))), levels=as.character(1:(2*nrow(pancreas_sig_overlap))))

pdf('/well/got2d/apayne/thesis/chapter4/gtex_islet_mod_overlap_genes_expression.pdf', width=5.9, height=5.9)
ggplot(overlap_sums, aes(x=X, y=mean_expr, shape=factor(Tissue), size=factor(Tissue))) +
  geom_point(aes(colour=factor(T2D))) +
  geom_vline(xintercept=seq(2.5, nrow(pancreas_sig_overlap)*2-1.5, 2), size=0.2) +
  scale_shape_manual(values=c(16, 17)) +
  scale_colour_manual(values=c('darkblue', 'darkred')) +
  scale_size_manual(values=c(2, 2)) +
  scale_x_discrete(labels=unique(overlap_sums$Gene), breaks=seq(1, nrow(pancreas_sig_overlap)*2-1, 2)) +
  xlab(NULL) +
  theme(panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text.x=element_text(angle=70, hjust=1),
    axis.ticks.x=element_blank(),
    legend.position='top'
  ) +
  guides(colour=guide_legend(title='Status', override.aes=list(size=c(4, 4))), shape=guide_legend(title='Tissue', override.aes=list(size=c(4, 2))), size='none') +
  labs(y='Mean normalised expression')
dev.off()

panc_diff <- c()
islet_diff <- c()

for (i in seq(1, nrow(overlap_sums), 4)) {
  panc_diff <- c(panc_diff, overlap_sums[i+1,'mean_expr']-overlap_sums[i,'mean_expr'])
  islet_diff <- c(islet_diff, overlap_sums[i+3,'mean_expr']-overlap_sums[i+2,'mean_expr'])
}

sum(sign(islet_diff*panc_diff)==1)  


#######################
# MR WITH GWAS TRAITS #
#######################
library(data.table)
library(rtracklayer)
library(edgeR)
library(MendelianRandomization)
library(glmnet)
library(parallel)

load('/well/got2d/apayne/InsPIRE_data/genotypes_maf05.Rda')
genotypes <- t(na.omit(t(genotypes)))

gc()

filelist <- list.files('/well/got2d/apayne/thesis/chapter4/islet_work/significant_modules_05_bh/')
filelist <- substring(filelist, 1, nchar(filelist)-22)
gencode <- readGFF('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz')
gencode <- gencode[!gencode$seqid=='MT',]
gencode$seqid <- as.character(gencode$seqid)

lncs <- subset(gencode, gene_type%in%c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA"))

phenotype <- read.table('/well/got2d/mvdbunt/InsPIRE/current_analysis/Cov_AllFilesIslets.csv', header=TRUE, stringsAsFactors=FALSE, sep=",")

snplocs <- fread('/well/got2d/apayne/InsPIRE_data/snplocs_full.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

snplocs <- snplocs[match(colnames(genotypes), snplocs$ID),c('CHR', 'POS', 'ID', 'REF', 'ALT')]
colnames(snplocs) <- c('CHROM', 'POS', 'ID147', 'REF', 'ALT')

snplocs <- snplocs[nchar(snplocs$ALT)==1&nchar(snplocs$REF)==1,]
genotypes <- genotypes[,colnames(genotypes)%in%snplocs$ID147]
gc()

bmi <- read.table('/well/got2d/apayne/GTEx_v7/200417_talk/gwas_data/BMI.ACTIVE.ALL.AllAncestry.txt.gz', sep='\t', header=TRUE, stringsAsFactors=FALSE)
diagram <- read.table('/well/got2d/apayne/GTEx_v7/200417_talk/gwas_data/DIAGRAMv3.2012DEC17.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
hdl <- read.table('/well/got2d/apayne/GTEx_v7/200417_talk/gwas_data/jointGwasMc_HDL.txt.gz', sep='\t', header=TRUE, stringsAsFactors=FALSE)
ldl <- read.table('/well/got2d/apayne/GTEx_v7/200417_talk/gwas_data/jointGwasMc_LDL.txt.gz', sep='\t', header=TRUE, stringsAsFactors=FALSE)
tc <- read.table('/well/got2d/apayne/GTEx_v7/200417_talk/gwas_data/jointGwasMc_TC.txt.gz', sep='\t', header=TRUE, stringsAsFactors=FALSE)
tg <- read.table('/well/got2d/apayne/GTEx_v7/200417_talk/gwas_data/jointGwasMc_TG.txt.gz', sep='\t', header=TRUE, stringsAsFactors=FALSE)
glucose <- read.table('/well/got2d/apayne/GTEx_v7/200417_talk/gwas_data/MAGIC_Manning_et_al_FastingGlucose_MainEffect.txt.gz', sep='\t', header=TRUE, stringsAsFactors=FALSE)
insulin <- read.table('/well/got2d/apayne/GTEx_v7/200417_talk/gwas_data/MAGIC_Manning_et_al_lnFastingInsulin_MainEffect.txt.gz', sep='\t', header=TRUE, stringsAsFactors=FALSE)

bmi <- subset(bmi, rsid%in%snplocs$ID147)
bmi$ALT <- toupper(bmi$Effect_allele)
bmi[bmi$ALT!=snplocs[match(bmi$rsid, snplocs$ID147),'ALT'],'Effect'] <- -bmi[bmi$ALT!=snplocs[match(bmi$rsid, snplocs$ID147),'ALT'],'Effect']
bmi <- bmi[,c('rsid', 'Effect', 'Stderr')]
colnames(bmi) <- c('ID', 'beta', 'se')

diagram <- subset(diagram, SNP%in%snplocs$ID147)
diagram$effect <- log(diagram[,'OR'])
diagram$se <- (log(diagram[,'OR_95U'])-log(diagram[,'OR_95L']))/(2*1.96)
diagram[diagram$RISK_ALLELE!=snplocs[match(diagram$SNP, snplocs$ID147),'ALT'],'effect'] <- -diagram[diagram$RISK_ALLELE!=snplocs[match(diagram$SNP, snplocs$ID147),'ALT'],'effect']
diagram <- diagram[,c('SNP', 'effect', 'se')]
colnames(diagram) <- c('ID', 'beta', 'se')

hdl <- subset(hdl, rsid%in%snplocs$ID147)
hdl$ALT <- toupper(hdl$A1)
hdl[hdl$ALT!=snplocs[match(hdl$rsid, snplocs$ID147),'ALT'],'beta'] <- -hdl[hdl$ALT!=snplocs[match(hdl$rsid, snplocs$ID147),'ALT'],'beta']
hdl <- hdl[,c('rsid', 'beta', 'se')]
colnames(hdl) <- c('ID', 'beta', 'se')

ldl <- subset(ldl, rsid%in%snplocs$ID147)
ldl$ALT <- toupper(ldl$A1)
ldl[ldl$ALT!=snplocs[match(ldl$rsid, snplocs$ID147),'ALT'],'beta'] <- -ldl[ldl$ALT!=snplocs[match(ldl$rsid, snplocs$ID147),'ALT'],'beta']
ldl <- ldl[,c('rsid', 'beta', 'se')]
colnames(ldl) <- c('ID', 'beta', 'se')

tc <- subset(tc, rsid%in%snplocs$ID147)
tc$ALT <- toupper(tc$A1)
tc[tc$ALT!=snplocs[match(tc$rsid, snplocs$ID147),'ALT'],'beta'] <- -tc[tc$ALT!=snplocs[match(tc$rsid, snplocs$ID147),'ALT'],'beta']
tc <- tc[,c('rsid', 'beta', 'se')]
colnames(tc) <- c('ID', 'beta', 'se')

tg <- subset(tg, rsid%in%snplocs$ID147)
tg$ALT <- toupper(tg$A1)
tg[tg$ALT!=snplocs[match(tg$rsid, snplocs$ID147),'ALT'],'beta'] <- -tg[tg$ALT!=snplocs[match(tg$rsid, snplocs$ID147),'ALT'],'beta']
tg <- tg[,c('rsid', 'beta', 'se')]
colnames(tg) <- c('ID', 'beta', 'se')

glucose <- subset(glucose, Snp%in%snplocs$ID147)
glucose$ALT <- toupper(glucose$effect_allele)
glucose[glucose$ALT!=snplocs[match(glucose$Snp, snplocs$ID147),'ALT'],'MainEffects'] <- -glucose[glucose$ALT!=snplocs[match(glucose$Snp, snplocs$ID147),'ALT'],'MainEffects']
glucose <- glucose[,c('Snp', 'MainEffects', 'MainSE')]
colnames(glucose) <- c('ID', 'beta', 'se')

insulin <- subset(insulin, Snp%in%snplocs$ID147)
insulin$ALT <- toupper(insulin$effect_allele)
insulin[insulin$ALT!=snplocs[match(insulin$Snp, snplocs$ID147),'ALT'],'MainEffects'] <- -insulin[insulin$ALT!=snplocs[match(insulin$Snp, snplocs$ID147),'ALT'],'MainEffects']
insulin <- insulin[,c('Snp', 'MainEffects', 'MainSE')]
colnames(insulin) <- c('ID', 'beta', 'se')

fullsnps <- unique(c(bmi$ID, diagram$ID, ldl$ID, hdl$ID, tc$ID, tg$ID, insulin$ID, glucose$ID))
fullsnps <- fullsnps[fullsnps%in%bmi$ID&fullsnps%in%diagram$ID&fullsnps%in%ldl$ID&fullsnps%in%hdl$ID&fullsnps%in%tc$ID&fullsnps%in%tg$ID&fullsnps%in%insulin$ID&fullsnps%in%glucose$ID]

genotypes <- genotypes[,colnames(genotypes)%in%fullsnps]

rm(gencode, lncs, snplocs, fullsnps)
gc()

eigensums <- NULL
curreigen <- read.table('/well/got2d/apayne/islet_t2d_networks/EG.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
currenmods <- read.table('/well/got2d/apayne/thesis/chapter4/islet_work/significant_modules_05_bh/islet.lncRNA_de_modules.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

sampids <- curreigen[,1]
rownames(curreigen) <- sampids
curreigen <- curreigen[,-1]

sampids <- sampids[sampids%in%rownames(genotypes)]
curreigen <- curreigen[match(sampids, rownames(curreigen)),]
currgenotypes <- genotypes[match(sampids, rownames(genotypes)),]
currphenotypes <- phenotype[match(sampids, phenotype$NewSampleID),]
gc()

t2d <- currphenotypes$T2D

for (i in unique(currenmods$color)) {
  print(paste('Starting work for eigengene ', i, sep=''))
  eigenmodel <- lm(curreigen[,i] ~ as.factor(t2d))

  lassofit <- cv.glmnet(x=currgenotypes, y=curreigen[,i], alpha=1, nfolds=10, keep=TRUE, dfmax=192)
  gc()

  MSEs <- lassofit$cvm
  nsnps <- lassofit$nzero
  ls <- lassofit$lambda

  names(MSEs) <- names(nsnps) <- ls

  lmin <- as.numeric(names(which.min(MSEs[nsnps>0])))

  coefs <- as.matrix(coef(lassofit, s=lmin))
  coefs <- coefs[coefs[,1]!=0,]

  if (length(coefs)==1) {
    current_summaries <- data.frame(tissue=currtiss, module_colour=i, cor.t2d=cor.test(curreigen[,i], t2d)$estimate, cor.t2d.p=cor.test(curreigen[,i], t2d)$p.value, t2d.beta=summary(eigenmodel)$coefficients[2,1], t2d.beta.se=summary(eigenmodel)$coefficients[2,2], t2d.beta.p=summary(eigenmodel)$coefficients[2,4], grs.eigen.r2=0)
    current_summaries[,1] <- as.character(current_summaries[,1])
    current_summaries[,2] <- as.character(current_summaries[,2])

    for (j in c('bmi', 'diagram', 'hdl', 'ldl', 'tc', 'tg', 'glucose', 'insulin')) {
      current_summaries$a <- 0
      current_summaries$b <- 1
      current_summaries$c <- 1
      colnames(current_summaries)[length(current_summaries)-(2:0)] <- paste(j, '.', c('mr.beta', 'mr.se', 'mr.p'), sep='')
    }

    return(current_summaries)
  }

  if (length(coefs)==2) {
    coefgenotypes <- currgenotypes[,names(coefs)[-1]]
    final_snps <- names(coefs)[-1]
  } else {
    coefgenotypes <- currgenotypes[,names(coefs)[-1]]
    lm_inter <- lm(curreigen[,i] ~ coefgenotypes)
    intersnps <- names(lm_inter$coefficients)[!is.na(lm_inter$coefficients)][-1]
    final_snps <- substring(intersnps, 14, nchar(intersnps))
    coefgenotypes <- coefgenotypes[,final_snps]
  }

  current_summaries <- data.frame(tissue='islet', module_colour=i, cor.t2d=cor.test(curreigen[,i], as.numeric(as.factor(t2d)))$estimate, cor.t2d.p=cor.test(curreigen[,i], as.numeric(as.factor(t2d)))$p.value, t2d.beta=summary(eigenmodel)$coefficients[2,1], t2d.beta.se=summary(eigenmodel)$coefficients[2,2], t2d.beta.p=summary(eigenmodel)$coefficients[2,4])
  current_summaries[,1] <- as.character(current_summaries[,1])
  current_summaries[,2] <- as.character(current_summaries[,2])

  lm_risk <- lm(curreigen[,i] ~ coefgenotypes)

  bx <- summary(lm_risk)$coefficients[-1,1]
  bxse <- summary(lm_risk)$coefficients[-1,2]

  current_summaries$grs.eigen.r2 <- summary(lm_risk)$adj.r.squared

  for (j in c('bmi', 'diagram', 'hdl', 'ldl', 'tc', 'tg', 'glucose', 'insulin')) {
    gc()

    mr_object <- mr_input(bx=bx, bxse=bxse, by=get(j)[match(final_snps, get(j)$ID),'beta'], byse=get(j)[match(final_snps, get(j)$ID),'se'], exposure=paste('eigen_', i, sep=''), outcome=j, snps=final_snps)
    currmr <- mr_ivw(mr_object, weights='delta')

    current_summaries$a <- currmr$Estimate
    current_summaries$b <- currmr$StdError
    current_summaries$c <- currmr$Pvalue

    colnames(current_summaries)[length(current_summaries)-(2:0)] <- paste(j, '.', c('mr.beta', 'mr.se', 'mr.p'), sep='')
  }

}

write.table(current_summaries, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE, file='/well/got2d/apayne/thesis/chapter4/islet_work/mr_and_correlation_2smaple/islet.mr_results.txt')

current_summaries_ps <- c(current_summaries[1,seq(11, 32, 3)])
current_summaries_ps_adj <- p.adjust(current_summaries_ps[1,], method='BH')


##########################
# GO ENRICHMENT ANALYSIS #
##########################
library("org.Hs.eg.db")
library(GOstats)
library(edgeR)
library(rtracklayer)
library(data.table)

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
environensg <- unlist(strsplit(rownames(dge), '.', fixed=TRUE))[seq(1, 2*nrow(dge), 2)]
environmap <- select(org.Hs.eg.db, keys=environensg, columns=c('ENSEMBL', 'ENTREZID'), keytype='ENSEMBL')
universeentrez <- environmap[match(environensg, environmap$ENSEMBL),'ENTREZID']
universeentrez <- unique(universeentrez[!is.na(universeentrez)])
modules <- read.table('/well/got2d/apayne/thesis/chapter4/islet_work/significant_modules_05_bh/islet.lncRNA_de_modules.txt', stringsAsFactors=FALSE, header=TRUE)

go_final <- NULL

for (j in unique(modules$module)) {
  selectedensg <- subset(modules, module==j)$gene_id
  selectedensg <- unlist(strsplit(selectedensg, '.', fixed=TRUE))[seq(1, 2*length(selectedensg), 2)]
  selectedentrez <- environmap[match(selectedensg, environmap$ENSEMBL),'ENTREZID']
  selectedentrez <- unique(selectedentrez[!is.na(selectedentrez)])
  params <- new('GOHyperGParams', geneIds=selectedentrez, universeGeneIds=universeentrez, annotation='org.Hs.eg.db', ontology='BP', pvalueCutoff=1, conditional=FALSE, testDirection='over')
  hgover <- summary(hyperGTest(params))
  hgover$module <- rep(j, nrow(hgover))    
  go_final <- rbind(go_final, hgover)
}

write.table(go_final, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter4/islet_work/go_enrichments/islet.go.txt')

full_GO_enrichment <- read.table('/well/got2d/apayne/thesis/chapter4/islet_work/go_enrichments/islet.go.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

full_GO_enrichment$p.adj <- p.adjust(full_GO_enrichment$Pvalue, method='BH')

full_GO_significant <- subset(full_GO_enrichment, p.adj<=0.05)
length(unique(full_GO_significant$Term))

insulin_inds <- c()
glucose_inds <- c()
glucagon_inds <- c()

for (i in 1:nrow(full_GO_significant)) {
  if (grepl('insulin', full_GO_significant[i,'Term'])) insulin_inds <- c(insulin_inds, i)
#  if (grepl('Insulin', full_GO_significant[i,'Term'])) insulin_inds <- c(insulin_inds, i)
  if (grepl('glucose', full_GO_significant[i,'Term'])) glucose_inds <- c(glucose_inds, i)
#  if (grepl('Glucose', full_GO_significant[i,'Term'])) glucose_inds <- c(glucose_inds, i)
  if (grepl('glucagon', full_GO_significant[i,'Term'])) glucagon_inds <- c(glucagon_inds, i)
}

t2d_GO <- full_GO_significant[unique(c(insulin_inds, glucose_inds)),]


#######################
# lncRNA MR on module #
#######################
library(data.table)
library(rtracklayer)
library(limma)
library(edgeR)
library(sva)

resultframe <- fread('/well/got2d/apayne/islet_t2d_networks/de_by_t2d.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

fullmodlist <- read.table('/well/got2d/apayne/thesis/chapter4/islet_work/significant_modules_05_bh/islet.lncRNA_de_modules.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
sigmodlist <- fullmodlist[match(unique(fullmodlist$module), fullmodlist$module),]

sigmod <- fullmodlist
sigmod_des <- subset(resultframe, gene_id%in%sigmod$gene_id)

phenotype <- read.table('/well/got2d/mvdbunt/InsPIRE/current_analysis/Cov_AllFilesIslets.csv', header=TRUE, stringsAsFactors=FALSE, sep=",")

#######################
degene <- 'ENSG00000254226.1'
#######################
# degene <- 'ENSG00000152254.6' #G6PC2
# degene <- 'ENSG00000010282.10' #HHATL
# degene <- 'ENSG00000138796.11' #HADH
# degene <- 'ENSG00000121351.3' #IAPP
# degene <- 'ENSG00000135447.12' #PPP1R1A

gencode <- readGFF('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz')
gencode$seqid <- as.character(gencode$seqid)
gencode <- gencode[gencode$seqid!='MT',]

expr <- fread('/well/got2d/apayne/InsPIRE_data/islet.overlapsremoved.gene.reads.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- expr[rownames(expr)%in%gencode$gene_id,]
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

newvals <- matrix(0, nrow=nrow(v), ncol=ncol(v))
rownames(newvals) <- rownames(v)
colnames(newvals) <- colnames(v)

for (k in 1:nrow(v)) newvals[k,] <- resid(lm(v[k,] ~ svars + sex))

load('/well/got2d/apayne/InsPIRE_data/genotypes_maf05.Rda')
genotypes <- genotypes[match(colnames(newvals), rownames(genotypes)),]
genotypes <- genotypes[!is.na(rownames(genotypes)),]
gc()
genotypes <- t(na.omit(t(genotypes)))
gc()

snplocs <- fread('/well/got2d/apayne/InsPIRE_data/snplocs_full.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

snplocs <- snplocs[match(colnames(genotypes), snplocs$ID),c('CHR', 'POS', 'ID', 'REF', 'ALT')]
colnames(snplocs) <- c('CHROM', 'POS', 'ID147', 'REF', 'ALT')

snplocs <- snplocs[nchar(snplocs$ALT)==1&nchar(snplocs$REF)==1,]
genotypes <- genotypes[,colnames(genotypes)%in%snplocs$ID147]
gc()

lncstart <- gencode[match(degene, gencode$gene_id),'start']
lncend <- gencode[match(degene, gencode$gene_id),'end']
lncchr <- gencode[match(degene, gencode$gene_id),'seqid']

snplocs$CHROM <- paste('chr', snplocs$CHROM, sep='')
psnps <- subset(snplocs, CHROM==lncchr)

lncgenos <- genotypes[,match(psnps$ID147, colnames(genotypes))]
testvals <- newvals[,match(rownames(lncgenos), colnames(newvals))]

#rm(design, dge, expr, gencode, genereads, genotypes, j, k, keeps, phenotype, psnps, sex, snplocs, svars, t2d, v)
gc()

library(glmnet)

MSEs <- NULL
fitr2 <- NULL
nsnps <- NULL
cvfit <- cv.glmnet(x=lncgenos, y=testvals[degene,], alpha=1, nfolds=10, keep=TRUE, dfmax=192)
MSEs <- cbind(MSEs, cvfit$cvm)
ls <- cvfit$lambda
fitr2 <- cbind(fitr2, unlist(lapply(1:length(ls), FUN=function(x) summary(lm(testvals[degene,] ~ cvfit$fit.preval[,x]))$r.squared)))
nsnps <- cbind(nsnps, cvfit$nzero)

rownames(MSEs) <- ls
rownames(fitr2) <- ls
rownames(nsnps) <- ls

lmin <- as.numeric(names(which.min(rowMeans(data.frame(MSEs[nsnps[,1]>0,]), na.rm=TRUE))))
gc()

coefs <- as.matrix(coef(cvfit, s=lmin))
coefs <- coefs[coefs[,1]!=0,]

testsnps <- names(coefs)[-1]
lncmodel <- lm(testvals[degene,] ~ lncgenos[,testsnps])
lncb <- summary(lncmodel)$coefficients[-1,1]
lncse <- summary(lncmodel)$coefficients[-1,2]

library(MendelianRandomization)

full_results <- matrix(0, nrow=nrow(subset(sigmod_des, adj.P.Val<=0.05)), ncol=1)
rownames(full_results) <- c(subset(sigmod_des, gene_id!=degene&adj.P.Val<=0.05)$gene_id, 'eigengene')
colnames(full_results) <- c('2SLS')

library(AER)

for (i in subset(sigmod_des, gene_id!=degene&adj.P.Val<=0.05)$gene_id) {
  if (i%in%rownames(testvals)) {
    full_results[i,] <- summary(ivreg(testvals[i,] ~ testvals[degene,] | lncgenos[,testsnps]))$coefficients[2,4]
  }
}

eigen <- read.table('/well/got2d/apayne/islet_t2d_networks/EG.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
sigmod_eigen <- eigen[,sigmod[1,'color']]
names(sigmod_eigen) <- eigen$samp_id
sigmod_eigen_gt <- sigmod_eigen[match(colnames(testvals), names(sigmod_eigen))]

full_results['eigengene','2SLS'] <- summary(ivreg(sigmod_eigen_gt ~ testvals[degene,] | lncgenos[,testsnps]))$coefficients[2,4]
full_results <- as.data.frame(full_results)
full_results$p.adj <- p.adjust(full_results[,1], method='BH')


corvals <- matrix(0, nrow=nrow(subset(sigmod_des, adj.P.Val<=0.05)), ncol=3)
rownames(corvals) <- c(subset(sigmod_des, gene_id!=degene&adj.P.Val<=0.05)$gene_id, 'eigengene')
colnames(corvals) <- c('r', 'pval', 'p.adj')

for (i in subset(sigmod_des, gene_id!=degene&adj.P.Val<=0.05)$gene_id) {
  corvals[i,'r'] <- cor.test(newvals[i,], newvals[degene,])$estimate
  corvals[i,'pval'] <- cor.test(newvals[i,], newvals[degene,])$p.value
}

corvals[,'p.adj'] <- p.adjust(corvals[,2], method='BH')

#colnames(newvals) <- sampids
#cor(newvals[degene,], phenotype[match(sampids, phenotype$SUBJID),'AGE'])
#cor(newvals[degene,], phenotype[match(sampids, phenotype$SUBJID),'HGHT'])
#or(newvals[degene,], phenotype[match(sampids, phenotype$SUBJID),'WGHT'])
#or(newvals[degene,], phenotype[match(sampids, phenotype$SUBJID),'BMI'])
#ummary(glm(as.factor(phenotype[match(sampids, phenotype$SUBJID),'MHT2D']) ~ newvals[degene,], family=binomial(link='logit')))
#ummary(lm(newvals[degene,] ~ as.factor(phenotype[match(sampids, phenotype$SUBJID),'MHT2D'])))
#ummary(lm(sigmod_eigen ~ as.factor(phenotype[match(sampids, phenotype$SUBJID),'MHT2D'])))





















library(data.table)
library(rtracklayer)

snplocs <- fread('/well/got2d/apayne/InsPIRE_data/snplocs_full.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
gencode <- readGFF('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz')
gencode <- gencode[!gencode$seqid%in%c('MT'),]
gencode[,1] <- as.character(gencode[,1])
gencode[,3] <- as.character(gencode[,3])

lncstart <- gencode[match('ENSG00000254226.1', gencode$gene_id),'start']
lncend <- gencode[match('ENSG00000254226.1', gencode$gene_id),'end']
lncchr <- gencode[match('ENSG00000254226.1', gencode$gene_id),'seqid']

load('/well/got2d/apayne/InsPIRE_data/genotypes_maf05.Rda')

psnps <- subset(snplocs, CHR==lncchr&POS>(lncstart-1000000)&POS<(lncend+1000000))
psnps <- subset(psnps, ID%in%colnames(genotypes))

finalinds <- colnames(newvals)[colnames(newvals)%in%rownames(genotypes)]

lncgenos <- genotypes[match(finalinds, rownames(genotypes)),match(psnps$ID, colnames(genotypes))]
lncgenos <- t(na.omit(t(lncgenos)))
testvals <- newvals[,finalinds]

library(glmnet)

MSEs <- NULL
fitr2 <- NULL
cvfit <- cv.glmnet(x=lncgenos, y=testvals['ENSG00000254226.1',], alpha=0.8, nfolds=5, keep=TRUE)
MSEs <- cbind(MSEs, cvfit$cvm)
ls <- cvfit$lambda
fitr2 <- cbind(fitr2, unlist(lapply(1:length(ls), FUN=function(x) summary(lm(testvals['ENSG00000254226.1',] ~ cvfit$fit.preval[,x]))$r.squared)))
counter <- 0

while (counter<200) {
  print(paste(i, ' iteration ', counter, sep=''))
  counter <- counter + 1
  cvfit <- cv.glmnet(x=lncgenos, y=testvals['ENSG00000254226.1',], alpha=0.8, nfolds=5, keep=TRUE)
  MSEs <- cbind(MSEs, cvfit$cvm[match(ls, cvfit$lambda)])
  fitr2 <- cbind(fitr2, unlist(lapply(match(ls, cvfit$lambda), FUN=function(x) if (is.na(x)) return(NA) else summary(lm(testvals['ENSG00000254226.1',] ~ cvfit$fit.preval[,x]))$r.squared)))
}

rownames(MSEs) <- ls
rownames(fitr2) <- ls
lmin <- as.numeric(names(which.min(rowMeans(MSEs, na.rm=TRUE))))
gc()

coefs <- as.matrix(coef(cvfit, s=lmin))
coefs <- coefs[coefs[,1]!=0,]

testsnps <- names(coefs)[-1]
lncmodel <- lm(testvals['ENSG00000254226.1',] ~ lncgenos[,testsnps])
lncb <- summary(lncmodel)$coefficients[-1,1]
lncse <- summary(lncmodel)$coefficients[-1,2]

library(MendelianRandomization)

#full_results <- matrix(0, nrow=length(sigmod_des$gene_id[sigmod_des$gene_id!='ENSG00000254226.1']), ncol=1)
full_results <- matrix(0, nrow=nrow(subset(sigmod_des, FDR_overall<=0.01)), ncol=1)
#rownames(full_results) <- sigmod_des$gene_id[sigmod_des$gene_id!='ENSG00000254226.1']
rownames(full_results) <- subset(sigmod_des, FDR_overall<=0.01)$gene_id
colnames(full_results) <- c('2SLS')

library(AER)

#for (i in sigmod_des$gene_id[sigmod_des$gene_id!='ENSG00000254226.1']) {
for (i in subset(sigmod_des, FDR_overall<=0.01)$gene_id) {

  if (i%in%rownames(testvals)) {
    full_results[i,] <- summary(ivreg(testvals[i,] ~ testvals['ENSG00000254226.1',] | lncgenos[,testsnps]))$coefficients[2,4]
#  currmodel <- lm(testvals[i,] ~ lncgenos[,finalsnps])
#  currb <- summary(currmodel)$coefficients[-1,1]
#  currse <- summary(currmodel)$coefficients[-1,2]
#  mrin <- mr_input(bx=lncb, bxse=lncse, by=currb, byse=currse)
#  ress <- mr_allmethods(mrin)
#  full_results[i,] <- ress$Values[1:3,'P-value']
#  lassofit2 <- cv.glmnet(x=lncgenos, y=testvals[i,], alpha=1, nfolds=10, keep=TRUE, dfmax=min(50, length(sampids)-1))
#  gc()
#  coefs2 <- as.matrix(coef(lassofit2, s=lassofit2$lambda.min))
#  coefs2 <- coefs2[coefs2[,1]!=0,]
#  testsnps2 <- names(coefs2)[-1]
#  lncmodel2 <- lm(testvals[i,] ~ lncgenos[,testsnps2])
#  lncb2 <- summary(lncmodel2)$coefficients[-1,1]
#  lncse2 <- summary(lncmodel2)$coefficients[-1,2]
#  finalsnps <- c(unlist(strsplit(rownames(summary(lncmodel2)$coefficients[-1,]), ']', fixed=TRUE))[seq(2, 2*nrow(summary(lncmodel2)$coefficients[-1,]), 2)])
  }

}

summary(ivreg(sigmod_eigen[finalinds] ~ testvals['ENSG00000254226.1',] | lncgenos[,testsnps]))$coefficients[2,4]


gencode[match(c('ENSG00000123836.10', 'ENSG00000058335.11', 'ENSG00000182759.3', 'ENSG00000019505.3', 'ENSG00000254647.2', 'ENSG00000083067.18', 'ENSG00000215644.5'), gencode$gene_id),]


t2dcors <- c()

for (i in subset(sigmod_des, adj.P.Val<=0.05)$gene_id) {
  t2dcors <- c(t2dcors, cor(newvals[i,], t2d))
}

names(t2dcors) <- subset(sigmod_des, adj.P.Val<=0.05)$gene_name














































#######################
# Analysis of Go work #
#######################
library(data.table)
library(sva)
library(edgeR)
library(limma)
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

newvals <- matrix(0, nrow=nrow(v), ncol=ncol(v))
rownames(newvals) <- rownames(v)
colnames(newvals) <- colnames(v)

for (k in 1:nrow(v)) newvals[k,] <- resid(lm(v[k,] ~ svars + sex))


thresh <- 0.01

modules <- fread(paste('/well/got2d/apayne/islet_t2d_networks/lncRNA_de_enriched_modules_', thresh, '.txt', sep=''), stringsAsFactors=FALSE, data.table=FALSE, header=TRUE)
GO_enrichment <- fread(file=paste('/well/got2d/apayne/islet_t2d_networks/lncRNA_de_enriched_modules_go_', thresh, '.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)




gwas <- read.table('/well/got2d/apayne/gwas_datasets/gwas_catalog_v1.0-associations_e89_r2017-08-15.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE, quote='', comment.char='', as.is=TRUE)

t2d_relevant <- c("Type 2 diabetes", "Type 2 diabetes and other traits", "Type 2 diabetes (young onset) and obesity", "Type 2 diabetes (dietary heme iron intake interaction)", "Type 2 diabetes (age of onset)", "Diabetes related insulin traits")

c("Two-hour glucose challenge", "Type 2 diabetes (young onset) and obesity", "Type 2 diabetes (dietary heme iron intake interaction)", "Type 2 diabetes (age of onset)", "Diabetes related insulin traits", "Proinsulin levels", "Fasting insulin-related traits (interaction with BMI)", "Fasting insulin-related traits", "Homeostasis model assessment of insulin resistance (interaction)", "Fasting insulin (interaction)", "Peak insulin response", "Fasting plasma glucose", "Fasting glucose-related traits (interaction with BMI)", "Fasting glucose-related traits", "Two-hour glucose challenge", "Fasting plasma glucose (childhood)", "Fasting blood glucose", "Insulin resistance/response", "Insulin-like growth factors", "Insulin-related traits", "Insulin secretion rate", "Insulin disposition index", "Modified Stumvoll Insulin Sensitivity Index", "Modified Stumvoll Insulin Sensitivity Index (BMI interaction)", "Glucose homeostasis traits", "Glycated hemoglobin levels", "Glycemic traits")

t2d_genes <- gwas[gwas$DISEASE.TRAIT%in%t2d_relevant,'MAPPED_GENE']

GO_enrichment$p.adjusted <- p.adjust(GO_enrichment$Pvalue, method='BH')

GO_significant <- subset(GO_enrichment, p.adjusted<=thresh)
GO_plotting <- subset(GO_significant, Term%in%names(sort(table(GO_significant$Term), decreasing=TRUE))[1:10])

insulin_inds <- c()
glucose_inds <- c()

for (i in 1:nrow(GO_significant)) {
  if (grepl('insulin', GO_significant[i,'Term'])) insulin_inds <- c(insulin_inds, i)
  if (grepl('glucose', GO_significant[i,'Term'])) glucose_inds <- c(glucose_inds, i)
}

GO_sig_ins <- GO_significant[c(insulin_inds, glucose_inds),]

enrich_mods <- list()
enrich_mods_t2d_genes <- c()

for (i in unique(GO_sig_ins$color)) {
  enrich_mods[[i]] <- list()
  enrich_mods[[i]][['full_module']] <- subset(modules, color==i)
  enrich_mods[[i]][['de_lncs']] <- subset(modules, color==i&simple_gene_type=='lncRNA'&FDR_overall<=thresh)
  enrich_mods[[i]][['de_pcs']] <- subset(modules, color==i&simple_gene_type=='protein_coding'&FDR_overall<=thresh)
#  enrich_mods[[i]][['t2d_genes']] <- subset(modules, color==i&gene_name%in%t2d_loci)
}

#######################
# lncRNA MR on module #
#######################
library(data.table)
library(rtracklayer)

snplocs <- fread('/well/got2d/apayne/InsPIRE_data/snplocs_full.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode[,1] <- as.character(gencode[,1])
gencode[,3] <- as.character(gencode[,3])

lncstart <- gencode[match('ENSG00000254226.1', gencode$gene_id),'start']
lncend <- gencode[match('ENSG00000254226.1', gencode$gene_id),'end']
lncchr <- gencode[match('ENSG00000254226.1', gencode$gene_id),'seqid']

load('/well/got2d/apayne/InsPIRE_data/genotypes_maf05.Rda')

psnps <- subset(snplocs, CHR==lncchr&POS>(lncstart-1000000)&POS<(lncend+1000000))
psnps <- subset(psnps, ID%in%colnames(genotypes))

finalinds <- colnames(newvals)[colnames(newvals)%in%rownames(genotypes)]

lncgenos <- genotypes[match(finalinds, rownames(genotypes)),match(psnps$ID, colnames(genotypes))]
lncgenos <- t(na.omit(t(lncgenos)))
testvals <- newvals[,finalinds]

library(glmnet)

MSEs <- NULL
fitr2 <- NULL
cvfit <- cv.glmnet(x=lncgenos, y=testvals['ENSG00000254226.1',], alpha=0.8, nfolds=5, keep=TRUE)
MSEs <- cbind(MSEs, cvfit$cvm)
ls <- cvfit$lambda
fitr2 <- cbind(fitr2, unlist(lapply(1:length(ls), FUN=function(x) summary(lm(testvals['ENSG00000254226.1',] ~ cvfit$fit.preval[,x]))$r.squared)))
counter <- 0

while (counter<200) {
  print(paste(i, ' iteration ', counter, sep=''))
  counter <- counter + 1
  cvfit <- cv.glmnet(x=lncgenos, y=testvals['ENSG00000254226.1',], alpha=0.8, nfolds=5, keep=TRUE)
  MSEs <- cbind(MSEs, cvfit$cvm[match(ls, cvfit$lambda)])
  fitr2 <- cbind(fitr2, unlist(lapply(match(ls, cvfit$lambda), FUN=function(x) if (is.na(x)) return(NA) else summary(lm(testvals['ENSG00000254226.1',] ~ cvfit$fit.preval[,x]))$r.squared)))
}

rownames(MSEs) <- ls
rownames(fitr2) <- ls
lmin <- as.numeric(names(which.min(rowMeans(MSEs, na.rm=TRUE))))
gc()

coefs <- as.matrix(coef(cvfit, s=lmin))
coefs <- coefs[coefs[,1]!=0,]

testsnps <- names(coefs)[-1]
lncmodel <- lm(testvals['ENSG00000254226.1',] ~ lncgenos[,testsnps])
lncb <- summary(lncmodel)$coefficients[-1,1]
lncse <- summary(lncmodel)$coefficients[-1,2]

library(MendelianRandomization)

#full_results <- matrix(0, nrow=length(sigmod_des$gene_id[sigmod_des$gene_id!='ENSG00000254226.1']), ncol=1)
full_results <- matrix(0, nrow=nrow(subset(sigmod_des, FDR_overall<=0.01)), ncol=1)
#rownames(full_results) <- sigmod_des$gene_id[sigmod_des$gene_id!='ENSG00000254226.1']
rownames(full_results) <- subset(sigmod_des, FDR_overall<=0.01)$gene_id
colnames(full_results) <- c('2SLS')

library(AER)

#for (i in sigmod_des$gene_id[sigmod_des$gene_id!='ENSG00000254226.1']) {
for (i in subset(sigmod_des, FDR_overall<=0.01)$gene_id) {

  if (i%in%rownames(testvals)) {
    full_results[i,] <- summary(ivreg(testvals[i,] ~ testvals['ENSG00000254226.1',] | lncgenos[,testsnps]))$coefficients[2,4]
#  currmodel <- lm(testvals[i,] ~ lncgenos[,finalsnps])
#  currb <- summary(currmodel)$coefficients[-1,1]
#  currse <- summary(currmodel)$coefficients[-1,2]
#  mrin <- mr_input(bx=lncb, bxse=lncse, by=currb, byse=currse)
#  ress <- mr_allmethods(mrin)
#  full_results[i,] <- ress$Values[1:3,'P-value']
#  lassofit2 <- cv.glmnet(x=lncgenos, y=testvals[i,], alpha=1, nfolds=10, keep=TRUE, dfmax=min(50, length(sampids)-1))
#  gc()
#  coefs2 <- as.matrix(coef(lassofit2, s=lassofit2$lambda.min))
#  coefs2 <- coefs2[coefs2[,1]!=0,]
#  testsnps2 <- names(coefs2)[-1]
#  lncmodel2 <- lm(testvals[i,] ~ lncgenos[,testsnps2])
#  lncb2 <- summary(lncmodel2)$coefficients[-1,1]
#  lncse2 <- summary(lncmodel2)$coefficients[-1,2]
#  finalsnps <- c(unlist(strsplit(rownames(summary(lncmodel2)$coefficients[-1,]), ']', fixed=TRUE))[seq(2, 2*nrow(summary(lncmodel2)$coefficients[-1,]), 2)])
  }

}

summary(ivreg(sigmod_eigen[finalinds] ~ testvals['ENSG00000254226.1',] | lncgenos[,testsnps]))$coefficients[2,4]


gencode[match(c('ENSG00000123836.10', 'ENSG00000058335.11', 'ENSG00000182759.3', 'ENSG00000019505.3', 'ENSG00000254647.2', 'ENSG00000083067.18', 'ENSG00000215644.5'), gencode$gene_id),]


t2dcors <- c()

for (i in subset(sigmod_des, adj.P.Val<=0.05)$gene_id) {
  t2dcors <- c(t2dcors, cor(newvals[i,], t2d))
}

names(t2dcors) <- subset(sigmod_des, adj.P.Val<=0.05)$gene_name





















