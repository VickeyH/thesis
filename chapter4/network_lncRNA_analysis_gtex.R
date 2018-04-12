library(data.table)
library(rtracklayer)
library(edgeR)
library(ggplot2)

filelist <- list.files('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/coexpression_networks')
tisslist <- unique(unlist(strsplit(filelist, '.', fixed=TRUE))[seq(1, 3*length(filelist), 3)])
gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode <- gencode[!gencode$seqid%in%c('MT'),]

lncs <- subset(gencode, gene_type%in%c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA"))

phenotype <- read.table('/well/got2d/GTEx_v7/sample_annotations/GTEx_Analysis_2016-01-15_v7_SubjectPhenotypesDS.txt', header=TRUE, stringsAsFactors=FALSE, quote="\"", sep="\t", comment.char="")

#OPTING TO EXCLUDE BREAST FOR NOW, AS MORE THAN 50% OF GENES ARE DE
modules <- NULL
uniquemods <- NULL

for (i in tisslist) {
  currmod <- read.table(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/coexpression_networks/', i, '.modules.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  modules <- rbind(modules, cbind(currmod, rep(i, nrow(currmod))))
  uniquemods <- rbind(uniquemods, currmod[match(unique(currmod$color), currmod$color),])
}

colnames(modules) <- c('gene_id', 'module_number', 'color', 'gene_name', 'tissue')
modules$tissue <- as.character(modules$tissue)

modsums <- NULL

for (i in unique(modules$tissue)) {
  cursums <- subset(modules, color!='grey'&tissue==i)
  modsums <- rbind(modsums, data.frame(tissue=i, genes=nrow(cursums), modules=length(unique(cursums$color)), stringsAsFactors=FALSE))
}



mod_plot <- do.call('rbind', lapply(unique(modules$tissue), FUN=function(currtiss) {
  currmods <- subset(modules, tissue==currtiss&module_number!=0)
  currtab <- table(currmods$color)
  return(data.frame(Tissue=rep(currtiss, length(currtab)), Color=names(currtab), N_genes=c(currtab), stringsAsFactors=FALSE))
}))

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
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
    pushViewport(viewport(layout = grid.layout(nrow=1, ncol=2, widths=c(1.44, 0.6))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

p2 <- ggplot(mod_plot, aes(x=Tissue, y=N_genes)) +
  geom_boxplot(fill='darkorange') +
  scale_y_log10() +
  theme(axis.title.y=element_blank(),
  axis.text.y=element_blank()) +
  labs(y='Genes per modules') +
  coord_flip()

p1 <- ggplot(modsums, aes(x=tissue, y=modules)) +
  geom_bar(stat='identity', fill='darkblue') +
  theme(axis.title.y=element_blank()) +
  labs(y='Number of modules') +
  coord_flip()

pdf('/well/got2d/apayne/thesis/chapter4/gtex_mod_sizes.pdf', width=5.9, height=7)
multiplot(p1, p2, cols=2)
dev.off()


de <- NULL

for (i in tisslist) {
  outputs <- fread(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/differential_expression/', i, '.de_by_sex.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  outputs$tissue <- rep(i, nrow(outputs))
  de <- rbind(de, outputs)
}

de <- de[de$gene_id%in%gencode$gene_id,]
de$FDR_overall <- p.adjust(de$P.Value, method='BH')
de$sig <- de$FDR_overall<=0.01
de$gene_type <- gencode[match(de$gene_id, gencode$gene_id),'gene_type']
de$simple_gene_type <- rep(NA, nrow(de))
de[de$gene_type%in%lncs$gene_type,'simple_gene_type'] <- 'lncRNA'
de[de$gene_type%in%'protein_coding','simple_gene_type'] <- 'protein_coding'
de[is.na(de$simple_gene_type),'simple_gene_type'] <- 'other'



mod_enrichments <- function(currtiss) {
  currde <- subset(de, tissue==currtiss)
  currmods <- subset(modules, tissue==currtiss)

  mod_assigned <- subset(currmods, module_number!=0)
  de_assigned <- subset(currde, gene_id%in%mod_assigned$gene_id)

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

  for (i in 1:length(module_list)) modframe <- rbind(modframe, data.frame(tissue=currtiss, module=i, color=currmods[match(i, currmods$module_number),'color'], pval=module_list[[i]][[5]], stringsAsFactors=FALSE))

  modframe$p.adj <- p.adjust(modframe$pval, method='BH')

  significant_modules <- NULL

  if (nrow(subset(modframe, p.adj<0.05))>0) {

    for (i in subset(modframe, p.adj<0.05)$module) {
      significant_modules <- rbind(significant_modules, cbind(module_list[[i]][[1]], data.frame(module=rep(i, nrow(module_list[[i]][[1]])), color=rep(currmods[match(i, currmods$module_number),'color'], nrow(module_list[[i]][[1]])))))
    }

    write.table(significant_modules, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE, file=paste('/well/got2d/apayne/thesis/chapter4/significant_modules_05_bh/', currtiss, '.lncRNA_de_modules.txt', sep=''))
  }

}

library(parallel)
mclapply(unique(modules$tissue), FUN=mod_enrichments, mc.cores=48)


filelist <- list.files('/well/got2d/apayne/thesis/chapter4/significant_modules_05_bh/')

fullmods <- do.call('rbind', lapply(filelist, FUN=function(i) {
  a <- read.table(paste('/well/got2d/apayne/thesis/chapter4/significant_modules_05_bh/', i, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  a$tissue <- rep(i, nrow(a))
  return(a)
}))

uniques <- do.call('rbind', lapply(unique(fullmods$tissue), FUN=function(x) {
  b <- subset(fullmods, tissue==x)
  return(b[match(unique(b$module), b$module),c('module', 'color', 'tissue')])
}))

dim(uniques)


######################
# MR WTH GWAS TRAITS #
######################
library(data.table)
library(rtracklayer)
library(edgeR)
library(MendelianRandomization)
library(glmnet)
library(parallel)

load('/well/got2d/apayne/GTEx_v7/vcf_maf05/genotypes.Rda')
is.na(genotypes) <- genotypes==(-1)
genotypes <- genotypes[,!is.na(colnames(genotypes))]

gc()

filelist <- list.files('/well/got2d/apayne/thesis/chapter4/significant_modules_05_bh/')
filelist <- substring(filelist, 1, nchar(filelist)-22)
gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode <- gencode[!gencode$seqid=='MT',]
gencode$seqid <- as.character(gencode$seqid)

lncs <- subset(gencode, gene_type%in%c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA"))

phenotype <- read.table('/well/got2d/GTEx_v7/sample_annotations/GTEx_Analysis_2016-01-15_v7_SubjectPhenotypesDS.txt', header=TRUE, stringsAsFactors=FALSE, quote="\"", sep="\t", comment.char="")

snplocs <- fread('/well/got2d/apayne/GTEx_v7/vcf_maf05/snplocs_and_ids_maf_0.05.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

lookup_table <- fread('/well/got2d/GTEx_v7/genotypes/WGS/variant_calls/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_support_files/GTEx_Analysis_2016-01-15_v7_WGS_652ind_VarID_Lookup_Table.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

snplocs$REF <- lookup_table[match(snplocs$ID147, lookup_table$RS_ID_dbSNP147_GRCh37p13),'Ref_b37']
snplocs$ALT <- lookup_table[match(snplocs$ID147, lookup_table$RS_ID_dbSNP147_GRCh37p13),'ALT']

snplocs <- snplocs[match(colnames(genotypes), snplocs$ID147),c('CHROM', 'POS', 'ID147', 'REF', 'ALT')]

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

rm(gencode, lncs, snplocs, lookup_table, fullsnps)
gc()

for (currtiss in filelist[1:length(filelist)]) {
  gc()
  eigensums <- NULL
  curreigen <- read.table(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/coexpression_networks/', currtiss, '.EG.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  currenmods <- read.table(paste('/well/got2d/apayne/thesis/chapter4/significant_modules_05_bh/', currtiss, '.lncRNA_de_modules.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)

  sampids <- curreigen[,1]
  splitlen <- length(unlist(strsplit(sampids[1], '-', fixed=TRUE)))
  sampids <- paste(unlist(strsplit(sampids, '-', fixed=TRUE))[seq(1, splitlen*length(sampids), splitlen)], '-', unlist(strsplit(sampids, '-', fixed=TRUE))[seq(2, splitlen*length(sampids), splitlen)], sep='')
  rownames(curreigen) <- sampids
  curreigen <- curreigen[,-1]

  sampids <- sampids[sampids%in%rownames(genotypes)]
  curreigen <- curreigen[match(sampids, rownames(curreigen)),]
  currgenotypes <- genotypes[match(sampids, rownames(genotypes)),]
  currphenotypes <- phenotype[match(sampids, phenotype$SUBJID),]
  gc()
  nacols <- colSums(is.na(currgenotypes))>0
  gc()
  currgenotypes <- currgenotypes[,!nacols]

  t2d <- currphenotypes$MHT2D

  gc()

  modtest <- function(i) {#  for (i in unique(currenmods$color)) {
    print(paste('Starting work for eigengene ', i, ' in ', currtiss, '.', sep=''))
    gc()
    eigenmodel <- lm(curreigen[,i] ~ t2d)

    lassofit <- cv.glmnet(x=currgenotypes, y=curreigen[,i], alpha=1, nfolds=10, keep=TRUE, dfmax=min(50, length(sampids)-1))
    gc()
    coefs <- as.matrix(coef(lassofit, s=lassofit$lambda.min))
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

    current_summaries <- data.frame(tissue=currtiss, module_colour=i, cor.t2d=cor.test(curreigen[,i], t2d)$estimate, cor.t2d.p=cor.test(curreigen[,i], t2d)$p.value, t2d.beta=summary(eigenmodel)$coefficients[2,1], t2d.beta.se=summary(eigenmodel)$coefficients[2,2], t2d.beta.p=summary(eigenmodel)$coefficients[2,4])
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

    return(current_summaries)

  }

  print(currtiss)
  eigensums <- do.call('rbind', mclapply(unique(currenmods$color), FUN=modtest, mc.cores=4))

  write.table(eigensums, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE, file=paste('/well/got2d/apayne/thesis/chapter4/mr_and_correlation_2smaple/', currtiss, '.mr_results.txt', sep=''))

}


filelist <- list.files('/well/got2d/apayne/thesis/chapter4/mr_and_correlation_2smaple/')

fullmr <- NULL

for (i in filelist) {
  if (length(scan(paste('/well/got2d/apayne/thesis/chapter4/mr_and_correlation_2smaple/', i, sep=''), sep='\t', 'character'))==0) next
  curr <- read.table(paste('/well/got2d/apayne/thesis/chapter4/mr_and_correlation_2smaple/', i, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  fullmr <- rbind(fullmr, curr)
}

rownames(fullmr) <- paste(fullmr[,1], fullmr[,2],  sep=':')
mr_ps <- fullmr[,c(seq(11, ncol(fullmr), 3))]

adj.ps <- matrix(p.adjust(unlist(mr_ps), method='BH'), nrow=nrow(mr_ps), byrow=FALSE)
rownames(adj.ps) <- rownames(mr_ps)
colnames(adj.ps) <- colnames(mr_ps)

adj.ps <- cbind(as.matrix(fullmr$cor.t2d), adj.ps)
colnames(adj.ps)[1] <- 't2d.cor'

sig_mrs <- NULL

for (i in 1:nrow(adj.ps)) {

  for (j in 2:ncol(adj.ps)) {

    if (adj.ps[i,j]<=0.05&!is.na(adj.ps[i,j])) {
      curr <- data.frame(tiss=unlist(strsplit(rownames(adj.ps)[i], ':', fixed=TRUE))[1], mod=unlist(strsplit(rownames(adj.ps)[i], ':', fixed=TRUE))[2], t2dcor=adj.ps[i,1], trait=unlist(strsplit(colnames(adj.ps)[j], '.', fixed=TRUE))[1], adj.p=adj.ps[i,j], stringsAsFactors=FALSE)
      sig_mrs <- rbind(sig_mrs, curr)
    }

  }

}


##########################
# GO ENRICHMENT ANALYSIS #
##########################

library("org.Hs.eg.db")
library(GOstats)
library(edgeR)
library(rtracklayer)
library(data.table)

gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode <- gencode[!gencode$seqid%in%c('MT'),]

filelist <- list.files('/well/got2d/apayne/thesis/chapter4/significant_modules_05_bh/')

system2('mkdir', '/well/got2d/apayne/thesis/chapter4/go_enrichments/')

for (i in filelist) {
  curri <- substring(i, 1, nchar(i)-22)
  currexpr <- fread(paste('/well/got2d/apayne/GTEx_v7/gene_counts/', curri, '.reads.txt', sep=''), sep='\t', stringsAsFactors=FALSE, data.table=FALSE)
  rownames(currexpr) <- currexpr[,1]
  currexpr <- currexpr[,-c(1, 2)]
  currexpr <- as.matrix(currexpr)
  currexpr <- currexpr[rownames(currexpr)%in%gencode$gene_id,]
  currexpr <- currexpr[rownames(currexpr)%in%gencode$gene_id,]
  currgenereads <- currexpr[rowSums(currexpr>6)>=10,]
  gc()
  currdge <- DGEList(counts=currgenereads)
  currdge <- currdge[rowSums(cpm(currdge)>1, na.rm=TRUE)>=10,]
  environensg <- unlist(strsplit(rownames(currdge), '.', fixed=TRUE))[seq(1, 2*nrow(currdge), 2)]
  environmap <- select(org.Hs.eg.db, keys=environensg, columns=c('ENSEMBL', 'ENTREZID'), keytype='ENSEMBL')
  universeentrez <- environmap[match(environensg, environmap$ENSEMBL),'ENTREZID']
  universeentrez <- unique(universeentrez[!is.na(universeentrez)])
  currmodules <- read.table(paste('/well/got2d/apayne/thesis/chapter4/significant_modules_05_bh/', i, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)

  for (j in unique(currmodules$module)) {
    selectedensg <- subset(currmodules, module==j)$gene_id
    selectedensg <- unlist(strsplit(selectedensg, '.', fixed=TRUE))[seq(1, 2*length(selectedensg), 2)]
    selectedentrez <- environmap[match(selectedensg, environmap$ENSEMBL),'ENTREZID']
    selectedentrez <- unique(selectedentrez[!is.na(selectedentrez)])
    params <- new('GOHyperGParams', geneIds=selectedentrez, universeGeneIds=universeentrez, annotation='org.Hs.eg.db', ontology='BP', pvalueCutoff=1, conditional=FALSE, testDirection='over')
    hgover <- summary(hyperGTest(params))
    hgover$tissue <- rep(curri, nrow(hgover))
    hgover$module <- rep(j, nrow(hgover))    

    if (nrow(hgover)>0) write.table(hgover, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t', file=paste('/well/got2d/apayne/thesis/chapter4/go_enrichments/', curri, '.go.txt', sep=''))
  }
}

golist <- list.files('/well/got2d/apayne/thesis/chapter4/go_enrichments/')

full_GO_enrichment <- NULL

for (i in golist) full_GO_enrichment <- rbind(full_GO_enrichment, read.table(paste('/well/got2d/apayne/thesis/chapter4/go_enrichments/', i, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE))

full_GO_enrichment$p.adj <- 1

for (j in unique(full_GO_enrichment$tissue)) full_GO_enrichment[full_GO_enrichment$tissue==j,'p.adj'] <- p.adjust(full_GO_enrichment[full_GO_enrichment$tissue==j,'Pvalue'], method='BH')


full_GO_significant <- subset(full_GO_enrichment, p.adj<=0.05)
length(unique(full_GO_significant$Term))
length(unique(full_GO_significant$tissue))

insulin_inds <- c()
glucose_inds <- c()

for (i in 1:nrow(full_GO_significant)) {
  if (grepl('insulin', full_GO_significant[i,'Term'])) insulin_inds <- c(insulin_inds, i)
  if (grepl('Insulin', full_GO_significant[i,'Term'])) insulin_inds <- c(insulin_inds, i)
  if (grepl('glucose', full_GO_significant[i,'Term'])) glucose_inds <- c(glucose_inds, i)
  if (grepl('Glucose', full_GO_significant[i,'Term'])) glucose_inds <- c(glucose_inds, i)
}

t2d_GO <- full_GO_significant[unique(c(insulin_inds, glucose_inds)),]

filelist <- list.files('/well/got2d/apayne/thesis/chapter4/significant_modules_05_bh/')

sigmodlist <- NULL
fullmodlist <- NULL

for (currtiss in filelist) {
  modules <- read.table(paste('/well/got2d/apayne/thesis/chapter4/significant_modules_05_bh/', currtiss, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  modules$tissue <- rep(currtiss, nrow(modules))
  fullmodlist <- rbind(fullmodlist, modules)
  modules <- modules[match(unique(modules$module), modules$module),]
  sigmodlist <- rbind(sigmodlist, modules)
}

filelist <- list.files('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/differential_expression')
lncs <- subset(gencode, gene_type%in%c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA"))

resultframe <- NULL

for (currtiss in filelist) {
  outputs <- fread(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/differential_expression/', currtiss, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  outputs$tissue <- rep(substr(currtiss, 1, nchar(currtiss)-14), nrow(outputs))
  resultframe <- rbind(resultframe, outputs)
}

de <- resultframe

de$gene_type <- gencode[match(de$gene_id, gencode$gene_id),'gene_type']
de$simple_gene_type <- rep(NA, nrow(de))
de[de$gene_type%in%lncs$gene_type,'simple_gene_type'] <- 'lncRNA'
de[de$gene_type%in%'protein_coding','simple_gene_type'] <- 'protein_coding'
de[is.na(de$simple_gene_type),'simple_gene_type'] <- 'other'
de$CHR <- gencode[match(de$gene_id, gencode$gene_id),'seqid']

t2d_GO_mods <- do.call('rbind', lapply(match(unique(paste(t2d_GO[,8], t2d_GO[,9])), paste(t2d_GO[,8], t2d_GO[,9])), FUN=function(x) {
  currtiss <- t2d_GO[x,8]
  currmod <- t2d_GO[x,9]
  currde <- subset(de, tissue==currtiss)
  currgenes <- subset(fullmodlist, tissue==paste(currtiss, '.lncRNA_de_modules.txt', sep='')&module==currmod)
  currgenes <- cbind(currgenes, currde[match(currgenes$gene_id, currde$gene_id),c(4:5, 11:12)])
  return(currgenes)
}))

subset(t2d_GO_mods, adj.P.Val<=0.05)


###############################################
#### PANCREAS SIGNIFICANT SO PANCREAS HERE ####
###############################################

thresh <- 0.01

filelist <- list.files(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/lncRNA_de_enriched_modules_', thresh, '/', sep=''))

sigmodlist <- NULL
fullmodlist <- NULL

for (currtiss in filelist) {
  modules <- read.table(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/lncRNA_de_enriched_modules_', thresh, '/', currtiss, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  modules$tissue <- rep(currtiss, nrow(modules))
  fullmodlist <- rbind(fullmodlist, modules)
  modules <- modules[match(unique(modules$module), modules$module),]
  sigmodlist <- rbind(sigmodlist, modules)
}

sigmod <- subset(fullmodlist, tissue=='Pancreas.lncRNA_de_modules.txt'&module==22)
sigmod_des <- subset(resultframe, tissue=='Pancreas'&gene_id%in%sigmod$gene_id)
sigmod_eigen <- read.table('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/coexpression_networks/Pancreas.EG.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
sigmod_eigen <- sigmod_eigen[,sigmod[1,'color']]

insec <- read.table('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/ashg/0030073.txt', sep='\t', header=FALSE, stringsAsFactors=FALSE)[,2]

currtiss <- 'Pancreas.reads.txt'

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

if (length(unique(sex))!=1) {
  for (k in 1:nrow(v)) newvals[k,] <- resid(lm(v[k,] ~ svars + sex))
} else {
  for (k in 1:nrow(v)) newvals[k,] <- resid(lm(v[k,] ~ svars))
}

library(ggplot2)

cormat <- cor(t(newvals[match(sigmod$gene_id, rownames(newvals)),]))
colnames(cormat) <- rownames(cormat) <- gencode[match(colnames(cormat), gencode$gene_id),'gene_name']
cormat.d <- as.data.frame(cormat)
ord <- hclust(dist(cormat.d, method='euclidean'), method='ward.D')$order
cormat.d$Gene1 <- rownames(cormat.d)
cormat.m <- melt(cormat.d, id.vars='Gene1', variable.name='Gene2')

cormat.m$Gene2 <- factor(cormat.m$Gene2, levels=colnames(cormat)[ord], labels=colnames(cormat))
cormat.m$Gene1 <- factor(cormat.m$Gene1, levels=rownames(cormat)[ord], labels=rownames(cormat))

pdf('/well/got2d/apayne/thesis/chapter4/pancreas_module_heat.pdf', width=5.9, height=6.5)

ggplot(cormat.m, aes(Gene1, Gene2)) +
  geom_tile(colour='white',
    aes(fill = value)
  ) +
  scale_fill_gradient2(limits=c(-1, 1),
    low='darkblue',
    high='darkred'
  ) +
  theme(axis.text.x=element_text(size=8,
      angle=90, 
      hjust=1
    ), axis.text.y=element_text(size=8),
    legend.position='top',
    axis.title=element_blank()
  )

dev.off()


################################
# MR for lncRNA on other genes #
################################
library(data.table)
library(rtracklayer)
library(limma)
library(edgeR)
library(sva)

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

if (length(unique(sex))!=1) {
  for (k in 1:nrow(v)) newvals[k,] <- resid(lm(v[k,] ~ svars + sex))
} else {
  for (k in 1:nrow(v)) newvals[k,] <- resid(lm(v[k,] ~ svars))
}

snplocs <- fread('/well/got2d/apayne/GTEx_v7/vcf_maf05/snplocs_and_ids_maf_0.05.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode$seqid <- as.character(gencode$seqid)

lncstart <- gencode[match('ENSG00000254226.1', gencode$gene_id),'start']
lncend <- gencode[match('ENSG00000254226.1', gencode$gene_id),'end']
lncchr <- gencode[match('ENSG00000254226.1', gencode$gene_id),'seqid']

load('/well/got2d/apayne/GTEx_v7/vcf_maf05/genotypes.Rda')

psnps <- subset(snplocs, CHROM==lncchr)
psnps <- subset(snplocs, ID142%in%colnames(genotypes))

finalinds <- sampids[keeps]%in%rownames(genotypes)

lncgenos <- genotypes[match(sampids[keeps][finalinds], rownames(genotypes)),match(psnps$ID142, colnames(genotypes))]
lncgenos <- t(na.omit(t(lncgenos)))
testvals <- newvals[,finalinds]

rm(design, dge, expr, gencode, genereads, genotypes, j, k, keeps, MSEs, phenotype, psnps, sampids, sex, snplocs, svars, t2d, v)
gc()

library(glmnet)

MSEs <- NULL
fitr2 <- NULL
cvfit <- cv.glmnet(x=lncgenos, y=testvals['ENSG00000254226.1',], alpha=1, nfolds=10, keep=TRUE, dfmax=217)
MSEs <- cbind(MSEs, cvfit$cvm)
ls <- cvfit$lambda
fitr2 <- cbind(fitr2, unlist(lapply(1:length(ls), FUN=function(x) summary(lm(testvals['ENSG00000254226.1',] ~ cvfit$fit.preval[,x]))$r.squared)))
counter <- 0

while (counter<20) {
  gc()
  print(paste(counter, ' iteration ', counter, sep=''))
  counter <- counter + 1
  cvfit <- cv.glmnet(x=lncgenos, y=testvals['ENSG00000254226.1',], alpha=1, nfolds=10, keep=TRUE, dfmax=217)
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

#finalsnps <- c(unlist(strsplit(rownames(summary(lncmodel)$coefficients[-1,]), ']', fixed=TRUE))[seq(2, 2*nrow(summary(lncmodel)$coefficients[-1,]), 2)])

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













library(ggplot2)

emf('/well/got2d/apayne/GTEx_v7/sex_specific_analyses/plots/top_GO_terms_pie.emf', width=6, height=6)
ggplot(full_GO_plotting, aes(x=1, fill=Term)) + geom_bar(colour='black') + coord_polar(theta='y') + labs(title='Ten most common over-represented GO terms', x='', y='') + theme(axis.text=element_blank(), axis.ticks=element_blank())
dev.off()


#######################
# T2D module overlaps #
#######################

library(data.table)
library(rtracklayer)
library(ggplot2)

filelist <- list.files('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/coexpression_networks')
tisslist <- unique(unlist(strsplit(filelist, '.', fixed=TRUE))[seq(1, 3*length(filelist), 3)])
gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode <- gencode[!gencode$seqid%in%c('MT'),]

lncs <- subset(gencode, gene_type%in%c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA"))

de <- NULL
thresh <- 0.01

for (i in tisslist) {
  outputs <- fread(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/differential_expression/', i, '.de_by_sex.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  outputs$tissue <- rep(i, nrow(outputs))
  de <- rbind(de, outputs)
}

de <- de[de$gene_id%in%gencode$gene_id,]
de$FDR_overall <- p.adjust(de$P.Value, method='BH')
de$sig <- de$FDR_overall<=0.01
de$gene_type <- gencode[match(de$gene_id, gencode$gene_id),'gene_type']
de$simple_gene_type <- rep(NA, nrow(de))
de[de$gene_type%in%lncs$gene_type,'simple_gene_type'] <- 'lncRNA'
de[de$gene_type%in%'protein_coding','simple_gene_type'] <- 'protein_coding'
de[is.na(de$simple_gene_type),'simple_gene_type'] <- 'other'

de_sig <- subset(de, FDR_overall<=thresh)

thresh <- 0.01

filelist <- list.files(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/lncRNA_de_enriched_modules_', thresh, '/', sep=''))

sigmodlist <- NULL
fullmodlist <- NULL

for (currtiss in filelist) {
  modules <- read.table(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/lncRNA_de_enriched_modules_', thresh, '/', currtiss, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  modules$tissue <- rep(currtiss, nrow(modules))
  fullmodlist <- rbind(fullmodlist, modules)
  modules <- modules[match(unique(modules$module), modules$module),]
  sigmodlist <- rbind(sigmodlist, modules)
}

t2d_genes <- read.table('/well/got2d/apayne/230117_diagram_for_martijn/t2d_locus_names.txt', sep='\t', header=FALSE, stringsAsFactors=FALSE)


#########################################
# MR WITH GWAS #
#########################################

library(data.table)
library(rtracklayer)
library(edgeR)
#library(ggplot2)
library(MendelianRandomization)
library(glmnet)
library(parallel)

thresh <- 0.01

load('/well/got2d/apayne/GTEx_v7/vcf_maf05/genotypes.Rda')
is.na(genotypes) <- genotypes==(-1)
genotypes <- genotypes[,!is.na(colnames(genotypes))]

gc()

filelist <- list.files(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/lncRNA_de_enriched_modules_', thresh, sep=''))
filelist <- substring(filelist, 1, nchar(filelist)-22)
#gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
#gencode <- gencode[!gencode$seqid%in%c('MT', 'Y', 'X'),]
gencode <- gencode[!gencode$seqid=='MT',]
gencode$seqid <- as.character(gencode$seqid)

lncs <- subset(gencode, gene_type%in%c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA"))

phenotype <- read.table('/well/got2d/GTEx_v7/sample_annotations/GTEx_Analysis_2016-01-15_v7_SubjectPhenotypesDS.txt', header=TRUE, stringsAsFactors=FALSE, quote="\"", sep="\t", comment.char="")

snplocs <- fread('/well/got2d/apayne/GTEx_v7/vcf_maf05/snplocs_and_ids_maf_0.05.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

lookup_table <- fread('/well/got2d/GTEx_v7/genotypes/WGS/variant_calls/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_support_files/GTEx_Analysis_2016-01-15_v7_WGS_652ind_VarID_Lookup_Table.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

snplocs$REF <- lookup_table[match(snplocs$ID147, lookup_table$RS_ID_dbSNP147_GRCh37p13),'Ref_b37']
snplocs$ALT <- lookup_table[match(snplocs$ID147, lookup_table$RS_ID_dbSNP147_GRCh37p13),'ALT']

snplocs <- snplocs[match(colnames(genotypes), snplocs$ID147),c('CHROM', 'POS', 'ID147', 'REF', 'ALT')]

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

rm(gencode, lncs, snplocs, lookup_table, fullsnps)
gc()

for (currtiss in filelist[1:length(filelist)]) {
  gc()
  eigensums <- NULL
  curreigen <- read.table(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/coexpression_networks/', currtiss, '.EG.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  currenmods <- read.table(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/lncRNA_de_enriched_modules_', thresh, '/', currtiss, '.lncRNA_de_modules.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
#  colormatch <- read.table(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/coexpression_networks/', currtiss, '.modules.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, colClasses=c('NULL', 'integer', 'character', 'NULL'))
#  currenmods$color <- colormatch[match(currenmods$module, colormatch$module_number),'color']

  sampids <- curreigen[,1]
  splitlen <- length(unlist(strsplit(sampids[1], '-', fixed=TRUE)))
  sampids <- paste(unlist(strsplit(sampids, '-', fixed=TRUE))[seq(1, splitlen*length(sampids), splitlen)], '-', unlist(strsplit(sampids, '-', fixed=TRUE))[seq(2, splitlen*length(sampids), splitlen)], sep='')
  rownames(curreigen) <- sampids
  curreigen <- curreigen[,-1]

  sampids <- sampids[sampids%in%rownames(genotypes)]
  curreigen <- curreigen[match(sampids, rownames(curreigen)),]
  currgenotypes <- genotypes[match(sampids, rownames(genotypes)),]
  currphenotypes <- phenotype[match(sampids, phenotype$SUBJID),]
  gc()
  nacols <- colSums(is.na(currgenotypes))>0
  gc()
  currgenotypes <- currgenotypes[,!nacols]

  t2d <- currphenotypes$MHT2D

  gc()

  modtest <- function(i) {#  for (i in unique(currenmods$color)) {
    print(paste('Starting work for eigengene ', i, ' in ', currtiss, '.', sep=''))
    gc()
    eigenmodel <- lm(curreigen[,i] ~ t2d)

    lassofit <- cv.glmnet(x=currgenotypes, y=curreigen[,i], alpha=1, nfolds=10, keep=TRUE, dfmax=min(50, length(sampids)-1))
    gc()
    coefs <- as.matrix(coef(lassofit, s=lassofit$lambda.min))
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

    current_summaries <- data.frame(tissue=currtiss, module_colour=i, cor.t2d=cor.test(curreigen[,i], t2d)$estimate, cor.t2d.p=cor.test(curreigen[,i], t2d)$p.value, t2d.beta=summary(eigenmodel)$coefficients[2,1], t2d.beta.se=summary(eigenmodel)$coefficients[2,2], t2d.beta.p=summary(eigenmodel)$coefficients[2,4])
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

    return(current_summaries)

  }

  print(currtiss)
  eigensums <- do.call('rbind', mclapply(unique(currenmods$color), FUN=modtest, mc.cores=4))

  write.table(eigensums, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE, file=paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/mr_and_correlation_per_tissue_and_module_2sample_', thresh, '/', currtiss, '.mr_results.txt', sep=''))

}

thresh <- 0.01

filelist <- list.files(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/mr_and_correlation_per_tissue_and_module_2sample_', thresh, sep=''))

fullmr <- NULL

for (i in filelist) {
  if (length(scan(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/mr_and_correlation_per_tissue_and_module_2sample_', thresh, '/', i, sep=''), sep='\t', 'character'))==0) next
  curr <- read.table(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/mr_and_correlation_per_tissue_and_module_2sample_', thresh, '/', i, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  fullmr <- rbind(fullmr, curr)
}

rownames(fullmr) <- paste(fullmr[,1], fullmr[,2],  sep=':')
mr_ps <- fullmr[,c(seq(11, ncol(fullmr), 3))]

adj.ps <- matrix(p.adjust(unlist(mr_ps), method='BH'), nrow=nrow(mr_ps), byrow=FALSE)
rownames(adj.ps) <- rownames(mr_ps)
colnames(adj.ps) <- colnames(mr_ps)

adj.ps <- cbind(as.matrix(fullmr$cor.t2d), adj.ps)
colnames(adj.ps)[1] <- 't2d.cor'

sig_mrs <- NULL

for (i in 1:nrow(adj.ps)) {

  for (j in 2:ncol(adj.ps)) {

    if (adj.ps[i,j]<=thresh&!is.na(adj.ps[i,j])) {
      curr <- data.frame(tiss=unlist(strsplit(rownames(adj.ps)[i], ':', fixed=TRUE))[1], mod=unlist(strsplit(rownames(adj.ps)[i], ':', fixed=TRUE))[2], t2dcor=adj.ps[i,1], trait=unlist(strsplit(colnames(adj.ps)[j], '.', fixed=TRUE))[1], adj.p=adj.ps[i,j], stringsAsFactors=FALSE)
      sig_mrs <- rbind(sig_mrs, curr)
    }

  }

}

sort(table(sig_mrs$trait))

sig_nosex <- subset(sig_mrs, abs(sexcor)<0.1)

########################
# Shows nothing useful #
########################
library(ggplot2)

mrplot <- expand.grid(Module=rownames(mr_ps), GWAS=c('BMI', 'T2D', 'HDL', 'LDL', 'TC', 'TG', 'GLUCOSE', 'INSULIN'))
mrplot$pvalue <- unlist(log(mr_ps))

ggplot(data=mrplot, aes(x=GWAS, y=Module)) +
  geom_tile(aes(fill=pvalue))
########################
########################
########################











































##########################
# Various eigengene work #
##########################

library(data.table)
#library(rtracklayer)
library(edgeR)
#library(ggplot2)
library(MendelianRandomization)
library(glmnet)
thresh <- 0.01

load('/well/got2d/apayne/GTEx_v7/vcf_maf05/genotypes.Rda')
is.na(genotypes) <- genotypes==(-1)
genotypes <- genotypes[,!is.na(colnames(genotypes))]

gc()

filelist <- list.files(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/lncRNA_de_enriched_modules_', thresh, sep=''))
filelist <- substring(filelist, 1, nchar(filelist)-22)
#gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode <- fread('/well/got2d/apayne/GTEx_v7/gencode.v19.genes.v7.patched_contigs.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
gencode <- gencode[!gencode$seqid%in%c('MT'),]

lncs <- subset(gencode, gene_type%in%c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA"))

phenotype <- read.table('/well/got2d/GTEx_v7/sample_annotations/GTEx_Analysis_2016-01-15_v7_SubjectPhenotypesDS.txt', header=TRUE, stringsAsFactors=FALSE, quote="\"", sep="\t", comment.char="")

binary_covariates <- colnames(phenotype)[substring(colnames(phenotype), 1, 2)%in%c('MH', 'LB')]
binary_covariates <- binary_covariates[!binary_covariates%in%c('MHBLDDNDR', 'MHGENCMT', 'MHSRC', 'MHTTCMT')]
continuous_covariates <- c('BMI', 'HGHT', 'WGHT')

gc()


for (currtiss in filelist) {
  gc()
  eigensums <- NULL
  curreigen <- read.table(paste('/well/got2d/apayne/GTEx_v7/coexpression_networks/coexpression_networks_keeping_sex_effect/', currtiss, '.EG.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  currenmods <- read.table(paste('/well/got2d/apayne/GTEx_v7/sex_specific_analyses/lncRNA_de_enriched_modules_', thresh, '/', currtiss, '.lncRNA_de_modules.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  colormatch <- read.table(paste('/well/got2d/apayne/GTEx_v7/coexpression_networks/coexpression_networks_keeping_sex_effect/', currtiss, '.modules.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, colClasses=c('NULL', 'integer', 'character', 'NULL'))
  currenmods$color <- colormatch[match(currenmods$module, colormatch$module_number),'color']

  sampids <- curreigen[,1]
  splitlen <- length(unlist(strsplit(sampids[1], '-', fixed=TRUE)))
  sampids <- paste(unlist(strsplit(sampids, '-', fixed=TRUE))[seq(1, splitlen*length(sampids), splitlen)], '-', unlist(strsplit(sampids, '-', fixed=TRUE))[seq(2, splitlen*length(sampids), splitlen)], sep='')
  rownames(curreigen) <- sampids
  curreigen <- curreigen[,-1]

  sampids <- sampids[sampids%in%rownames(genotypes)]
  curreigen <- curreigen[match(sampids, rownames(curreigen)),]
  currgenotypes <- genotypes[match(sampids, rownames(genotypes)),]
  currphenotypes <- phenotype[match(sampids, phenotype$SUBJID),]
  gc()
  nacols <- colSums(is.na(currgenotypes))>0
  gc()
  currgenotypes <- currgenotypes[,!nacols]

  sex <- currphenotypes$SEX

  gc()  

  for (i in unique(currenmods$color)) {
    print(paste('Starting work for eigengene ', i, ' in ', currtiss, '.', sep=''))
    gc()
    eigenmodel <- lm(curreigen[,i] ~ sex)
    current_summaries <- data.frame(tissue=currtiss, module_colour=i, cor.sex=cor.test(curreigen[,i], sex)$estimate, cor.sex.p=cor.test(curreigen[,i], sex)$p.value, sex.beta=summary(eigenmodel)$coefficients[2,1], sex.beta.se=summary(eigenmodel)$coefficients[2,2], sex.beta.p=summary(eigenmodel)$coefficients[2,4])
    current_summaries[,1] <- as.character(current_summaries[,1])
    current_summaries[,2] <- as.character(current_summaries[,2])

    lassofit <- cv.glmnet(x=currgenotypes, y=curreigen[,i], alpha=1, nfolds=10, keep=TRUE, dfmax=min(50, length(sampids)-1))
    gc()
    coefs <- as.matrix(coef(lassofit, s=lassofit$lambda.min))
    coefs <- coefs[coefs[,1]!=0,]

    if (length(coefs)==1) next

    coefgenotypes <- currgenotypes[,names(coefs)[-1]]

    lm_inter <- lm(curreigen[,i] ~ coefgenotypes)
    intersnps <- names(lm_inter$coefficients)[!is.na(lm_inter$coefficients)][-1]
    final_snps <- substring(intersnps, 14, nchar(intersnps))

    coefgenotypes <- coefgenotypes[,final_snps]
    lm_risk <- lm(curreigen[,i] ~ coefgenotypes)

    bx <- summary(lm_risk)$coefficients[-1,1]
    bxse <- summary(lm_risk)$coefficients[-1,2]
    final_snps <- substring(names(bx), 14, nchar(names(bx)))

    current_summaries$grs.eigen.r2 <- summary(lm_risk)$adj.r.squared

    for (j in binary_covariates) {
      gc()
      currcov <- currphenotypes[,j]
      currcov_keep <- currcov%in%c(0, 1)
      currcov_outcome <- currcov[currcov_keep]
      currcov_genotypes <- currgenotypes[currcov_keep,final_snps]
      gc()

      if (length(unique(currcov_outcome))==1|length(currcov_outcome)==0) {
        current_summaries$a <- 0
        current_summaries$b <- 1
        current_summaries$c <- 1
        colnames(current_summaries)[length(current_summaries)-(2:0)] <- paste(j, '.', c('mr.beta', 'mr.se', 'mr.p'), sep='')
      } else if (ncol(currcov_genotypes)>=nrow(currcov_genotypes)) {
        current_summaries$a <- 0
        current_summaries$b <- 1
        current_summaries$c <- 1
        colnames(current_summaries)[length(current_summaries)-(2:0)] <- paste(j, '.', c('mr.beta', 'mr.se', 'mr.p'), sep='')
      } else {
        lm_outcome <- glm(currcov_outcome ~ currcov_genotypes, family='binomial', control=list(maxit=100))
        by <- summary(lm_outcome)$coefficients[-1,1]
        byse <- summary(lm_outcome)$coefficients[-1,2]

        outcome_snps <- substring(names(by), 18, nchar(names(by)))
        by <- by[match(final_snps, outcome_snps)]
        byse <- byse[match(final_snps, outcome_snps)]

        mr_object <- mr_input(bx=bx, bxse=bxse, by=by, byse=byse, exposure=paste('eigen_', i, sep=''), outcome=j, snps=final_snps)
        currmr <- mr_ivw(mr_object, weights='delta', psi=cor(curreigen[currcov_keep,i], currcov_outcome))

        current_summaries$a <- currmr$Estimate
        current_summaries$b <- currmr$StdError
        current_summaries$c <- currmr$Pvalue

        colnames(current_summaries)[length(current_summaries)-(2:0)] <- paste(j, '.', c('mr.beta', 'mr.se', 'mr.p'), sep='')
      }

    }

    for (j in continuous_covariates) {
      gc()
      currcov <- currphenotypes[,j]
      currcov_keep <- !is.na(currcov)
      currcov_outcome <- currcov[currcov_keep]
      currcov_genotypes <- currgenotypes[currcov_keep,final_snps]
      gc()

      lm_outcome <- lm(currcov_outcome ~ currcov_genotypes)
      by <- summary(lm_outcome)$coefficients[-1,1]
      byse <- summary(lm_outcome)$coefficients[-1,2]

      outcome_snps <- substring(names(by), 18, nchar(names(by)))
      by <- by[match(final_snps, outcome_snps)]
      byse <- byse[match(final_snps, outcome_snps)]

      mr_object <- mr_input(bx=bx, bxse=bxse, by=by, byse=byse, exposure=paste('eigen_', i, sep=''), outcome=j, snps=final_snps)
      currmr <- mr_ivw(mr_object, weights='delta', psi=cor(curreigen[currcov_keep,i], currcov_outcome))

      current_summaries$a <- currmr$Estimate
      current_summaries$b <- currmr$StdError
      current_summaries$c <- currmr$Pvalue

      colnames(current_summaries)[length(current_summaries)-(2:0)] <- paste(j, '.', c('mr.beta', 'mr.se', 'mr.p'), sep='')
    }

    eigensums <- rbind(eigensums, current_summaries)
  }

  write.table(eigensums, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE, file=paste('/well/got2d/apayne/GTEx_v7/sex_specific_analyses/mr_and_correlation_per_tissue_and_module/', currtiss, '.mr_results.txt', sep=''))
  print(paste(currtiss, ' completed.', sep=''))
}





thresh <- 0.01

filelist <- list.files('/well/got2d/apayne/GTEx_v7/sex_specific_analyses/mr_and_correlation_per_tissue_and_module/')

fullmr <- NULL

for (i in filelist) {
  curr <- read.table(paste('/well/got2d/apayne/GTEx_v7/sex_specific_analyses/mr_and_correlation_per_tissue_and_module/', i, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  fullmr <- rbind(fullmr, curr)
}

rownames(fullmr) <- paste(fullmr[,1], fullmr[,2], round(fullmr[,3], 2),  sep=':')
mr_ps <- fullmr[,c(seq(11, ncol(fullmr), 3))]

adj.ps <- matrix(p.adjust(unlist(mr_ps), method='BH'), nrow=nrow(mr_ps), byrow=FALSE)
rownames(adj.ps) <- rownames(mr_ps)
colnames(adj.ps) <- colnames(mr_ps)


sig_mrs <- NULL

for (i in 1:nrow(adj.ps)) {

  for (j in 1:ncol(adj.ps)) {

    if (adj.ps[i,j]<=thresh) {
      curr <- data.frame(tiss_mod_sex=rownames(adj.ps)[i], trait=colnames(adj.ps)[j], adj.p=adj.ps[i,j], stringsAsFactors=FALSE)
      sig_mrs <- rbind(sig_mrs, curr)
    }

  }

}

sort(table(sig_mrs$trait))

sig_nosex <- sig_mrs[sig_mrs$adj.p<=0.05&abs(as.numeric(unlist(strsplit(sig_mrs$tiss_mod_sex, ':', fixed=TRUE))[seq(3, 3*nrow(sig_mrs), 3)]))<0.8,]



library(data.table)
library(rtracklayer)
library(edgeR)
library(ggplot2)
library(MendelianRandomization)
library(glmnet)
thresh <- 0.05

gc()

filelist <- list.files(paste('/well/got2d/apayne/GTEx_v7/sex_specific_analyses/lncRNA_de_enriched_modules_', thresh, sep=''))
filelist <- substring(filelist, 1, nchar(filelist)-22)
gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode <- gencode[!gencode$seqid%in%c('MT'),]

load('/well/got2d/apayne/GTEx_v7/vcf_maf05/genotypes.Rda')
is.na(genotypes) <- genotypes==(-1)
genotypes <- genotypes[,!is.na(colnames(genotypes))]

lncs <- subset(gencode, gene_type%in%c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA"))

phenotype <- read.table('/well/got2d/GTEx_v7/sample_annotations/GTEx_Analysis_2016-01-15_v7_SubjectPhenotypesDS.txt', header=TRUE, stringsAsFactors=FALSE, quote="\"", sep="\t", comment.char="")

binary_covariates <- colnames(phenotype)[substring(colnames(phenotype), 1, 2)%in%c('MH', 'LB')]
binary_covariates <- binary_covariates[!binary_covariates%in%c('MHBLDDNDR', 'MHGENCMT', 'MHSRC', 'MHTTCMT')]
continuous_covariates <- c('BMI', 'HGHT', 'WGHT')

gc()

currtiss <- 'Thyroid'

curreigen <- read.table(paste('/well/got2d/apayne/GTEx_v7/coexpression_networks/coexpression_networks_keeping_sex_effect/', currtiss, '.EG.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
currenmods <- read.table(paste('/well/got2d/apayne/GTEx_v7/sex_specific_analyses/lncRNA_de_enriched_modules_', thresh, '/', currtiss, '.lncRNA_de_modules.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
colormatch <- read.table(paste('/well/got2d/apayne/GTEx_v7/coexpression_networks/coexpression_networks_keeping_sex_effect/', currtiss, '.modules.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, colClasses=c('NULL', 'integer', 'character', 'NULL'))
currenmods$color <- colormatch[match(currenmods$module, colormatch$module_number),'color']

sampids <- curreigen[,1]
splitlen <- length(unlist(strsplit(sampids[1], '-', fixed=TRUE)))
sampids <- paste(unlist(strsplit(sampids, '-', fixed=TRUE))[seq(1, splitlen*length(sampids), splitlen)], '-', unlist(strsplit(sampids, '-', fixed=TRUE))[seq(2, splitlen*length(sampids), splitlen)], sep='')
rownames(curreigen) <- sampids
curreigen <- curreigen[,-1]

curreigen <- curreigen[match(sampids, rownames(curreigen)),]
currgenotypes <- genotypes[match(sampids, rownames(genotypes)),]
currphenotypes <- phenotype[match(sampids, phenotype$SUBJID),]
nacols <- colSums(is.na(currgenotypes))>0
currgenotypes <- currgenotypes[,!nacols]

sex <- currphenotypes$SEX

gc()


















##############################################
##################################################
##############################################################
################################
# Gene set enrichment analysis #
################################

library(KEGGprofile)
library(biomaRt)

thresh <- 0.01

filelist <- list.files('/well/got2d/apayne/thesis/chapter4/significant_modules_05_bh/')

sigmodlist <- NULL
fullmodlist <- NULL

for (currtiss in filelist) {
  modules <- read.table(paste('/well/got2d/apayne/thesis/chapter4/significant_modules_05_bh/', currtiss, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  modules$tissue <- rep(currtiss, nrow(modules))
  fullmodlist <- rbind(fullmodlist, modules)
  modules <- modules[match(unique(modules$module), modules$module),]
  sigmodlist <- rbind(sigmodlist, modules)
}


filelist <- list.files('/well/got2d/apayne/thesis/chapter4/significant_modules_05_bh/')
mart <- useDataset('hsapiens_gene_ensembl', useEnsembl(biomart='ensembl', GRCh=37))

#for (currtiss in filelist) {
kegger <- function(currtiss) {
  modules <- read.table(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/lncRNA_de_enriched_modules_', thresh, '/', currtiss, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  ensembl_ids <- unlist(strsplit(modules$gene_id, '.', fixed=TRUE))[seq(1, 2*nrow(modules), 2)]
  convert_id <- getBM(filters='ensembl_gene_id', attributes=c('ensembl_gene_id', 'entrezgene'), values=ensembl_ids, mart=mart)
  modules$entrez <- convert_id[match(ensembl_ids, convert_id$ensembl_gene_id),'entrezgene']
  pathways_to_output <- NULL

  for (j in unique(modules$module)) {
    currmod <- subset(modules, module==j)
    curr_pathways <- find_enriched_pathway(currmod$entrez, species='hsa', , returned_pvalue=1, download_latest=TRUE)[[1]]
    curr_pathways$module <- rep(j, nrow(curr_pathways))
    curr_pathways$color <- rep(currmod[1,'color'], nrow(curr_pathways))
    pathways_to_output <- rbind(pathways_to_output, curr_pathways)
  }

  if (nrow(pathways_to_output)>0)  write.table(pathways_to_output, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, file=paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/lncRNA_de_enriched_modules_pathways_', thresh, '/', substring(currtiss, 1, nchar(currtiss)-22), '.pathways.txt', sep=''))

}

mclapply(filelist, FUN=kegger, mc.cores=48)

gsea_files <- list.files(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/lncRNA_de_enriched_modules_pathways_', thresh, '/', sep=''))

full_gsea <- NULL

for (currtiss in gsea_files) {
  currgsea <- read.table(paste('/well/got2d/apayne/GTEx_v7/t2d_specific_analyses/lncRNA_de_enriched_modules_pathways_', thresh, '/', currtiss, sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  currgsea$tissue <- rep(substring(currtiss, 1, nchar(currtiss)-13), nrow(currgsea))
  full_gsea <- rbind(full_gsea, currgsea)
}

full_gsea$fdr_overall <- p.adjust(full_gsea$pvalue, method='BH')



