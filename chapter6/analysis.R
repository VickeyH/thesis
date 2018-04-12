library(data.table)

fi <- fread('zcat /well/got2d/apayne/thesis/chapter6/gwas_datasets/MAGIC_Manning_et_al_lnFastingInsulin_MainEffect.txt.gz', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
fg <- fread('zcat /well/got2d/apayne/thesis/chapter6/gwas_datasets/MAGIC_Manning_et_al_FastingGlucose_MainEffect.txt.gz', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
t2d <- fread('/well/got2d/apayne/thesis/chapter6/gwas_datasets/METAANALYSIS_DIAGRAM_SE1.txt', sep='\t', header=FALSE, stringsAsFactors=FALSE, data.table=FALSE, skip=1)
colnames(t2d) <- c('chr_pos', 'Allele1', 'Allele2', 'Effect', 'StdErr', 'Pvalue', 'TotalSampleSize')
t2d$chr <- unlist(strsplit(t2d$chr_pos, ':', fixed=TRUE))[seq(1, 2*nrow(t2d), 2)]
t2d$pos <- unlist(strsplit(t2d$chr_pos, ':', fixed=TRUE))[seq(2, 2*nrow(t2d), 2)]
t2d$pos <- as.integer(t2d$pos)

load('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/filtered.Rda')

t2d_with_rsids <- merge(t2d, snplocs, by.x=c('chr', 'pos'), by.y=c('CHR', 'POS'))

fisnp <- fi$Snp
fgsnp <- fg$Snp
t2dsnp <- t2d_with_rsids$ID

final_snp_list <- fisnp[fisnp%in%fgsnp&fisnp%in%t2dsnp]

gt <- gt[,colnames(gt)%in%final_snp_list]
snplocs <- snplocs[snplocs$ID%in%final_snp_list,]

fi <- fi[match(snplocs$ID, fi$Snp),]
fg <- fg[match(snplocs$ID, fg$Snp),]
t2d <- t2d_with_rsids[match(snplocs$ID, t2d_with_rsids$ID),]

fi$effect_allele <- toupper(fi$effect_allele)
fi$other_allele <- toupper(fi$other_allele)
fg$effect_allele <- toupper(fg$effect_allele)
fg$other_allele <- toupper(fg$other_allele)

fi_all_matches <- fi[(fi$effect_allele==snplocs$REF|fi$effect_allele==snplocs$ALT)&(fi$other_allele==snplocs$REF|fi$other_allele==snplocs$ALT),]
fg_all_matches <- fg[(fg$effect_allele==snplocs$REF|fg$effect_allele==snplocs$ALT)&(fg$other_allele==snplocs$REF|fg$other_allele==snplocs$ALT),]
t2d_all_matches <- t2d[(t2d$Allele1==snplocs$REF|t2d$Allele1==snplocs$ALT)&(t2d$Allele2==snplocs$REF|t2d$Allele2==snplocs$ALT),]

fisnp <- fi_all_matches$Snp
fgsnp <- fg_all_matches$Snp
t2dsnp <- t2d_all_matches$ID

final_snp_list <- fisnp[fisnp%in%fgsnp&fisnp%in%t2dsnp]

gt <- gt[,colnames(gt)%in%final_snp_list]
snplocs <- snplocs[snplocs$ID%in%final_snp_list,]

fi <- fi[match(snplocs$ID, fi$Snp),]
fg <- fg[match(snplocs$ID, fg$Snp),]
t2d <- t2d_with_rsids[match(snplocs$ID, t2d_with_rsids$ID),]

fi$NewEffects <- fi$MainEffects
fi[fi$effect_allele!=snplocs$ALT,'NewEffects'] <- -fi[fi$effect_allele!=snplocs$ALT,'NewEffects']

fg$NewEffects <- fg$MainEffects
fg[fg$effect_allele!=snplocs$ALT,'NewEffects'] <- -fg[fg$effect_allele!=snplocs$ALT,'NewEffects']

t2d$NewEffects <- t2d$Effect
t2d[t2d$Allele1!=snplocs$ALT,'NewEffects'] <- -t2d[t2d$Allele1!=snplocs$ALT,'NewEffects']

fiwrite <- fi[,c('Snp', 'NewEffects', 'MainSE', 'MainP')]
fgwrite <- fg[,c('Snp', 'NewEffects', 'MainSE', 'MainP')]
t2dwrite <- t2d[,c('ID', 'NewEffects', 'StdErr', 'Pvalue')]

colnames(fiwrite) <- colnames(fgwrite) <- colnames(t2dwrite) <- c('ID', 'EFFECT', 'SE', 'P')

save(covs, final_genes, gt, snplocs, v, file='/well/got2d/apayne/thesis/chapter6/filtered.Rda')

fwrite(fiwrite, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/fi.txt')
fwrite(fgwrite, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/fg.txt')
fwrite(t2dwrite, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/t2d.txt')


##################
# EQTL FILTERING #
##################
eqtl <- fread('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/fullmatrixresults_residualised_expr_input.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
load('/well/got2d/apayne/thesis/chapter6/filtered.Rda')
eqtl <- subset(eqtl, SNP%in%colnames(gt))
eqtl$FDR_overall <- p.adjust(eqtl$"p-value", method='BH')
fwrite(eqtl, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/fulleqtls.txt')

topmodeller <- function(i) {
  fitmod <- lm(unlist(GE[GE[,1]==i,-1]) ~ unlist(SNP[SNP[,1]==tops[tops$gene==i,'SNP'],-1]))
  return(list(data.frame(gene=i, snp=tops[tops$gene==i,'SNP'], coefs=summary(fitmod)$coefficients[2,1], stringsAsFactors=FALSE),
    data.frame(gene=i, nsnps=1, fullinds.r2=summary(fitmod)$r.squared, fullinds.p=pf(summary(fitmod)$fstatistic[1], summary(fitmod)$fstatistic[2], summary(fitmod)$fstatistic[3], lower.tail=F))
  ))
}

fulltops <- eqtl[match(unique(eqtl$gene), eqtl$gene),]

for (j in as.character(1:22)) {
  tops <- subset(fulltops, CHR==j)
  GE <- fread(paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', j, '/GE_resid.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  SNP <- fread(paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', j, '/SNP.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  geneloc <- fread(paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', j, '/geneloc.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  snpsloc <- fread(paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', j, '/snpsloc.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

  outframes <- mclapply(tops$gene, FUN=topmodeller, mc.cores=56)
  coefframe <- rbindlist(mclapply(1:length(outframes), FUN=function(x) outframes[[x]][[1]], mc.cores=48))
  genesumframe <- rbindlist(mclapply(1:length(outframes), FUN=function(x) outframes[[x]][[2]], mc.cores=48))
  write.table(coefframe, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter6/lasso/eqtl_coef_sums_temp/chr_', j, '_top_eqtl_snps_and_coefficients.txt', sep=''))
  write.table(genesumframe, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter6/lasso/eqtl_coef_sums_temp/chr_', j, '_top_eqtl_gene_model_summaries.txt', sep=''))
  print(paste('Chromosome ', j, ' completed.', sep=''))
}

coefs <- NULL
sums <- NULL

for (j in as.character(1:22)) {
  currcoef <- read.table(paste('/well/got2d/apayne/thesis/chapter6/lasso/eqtl_coef_sums_temp/chr_', j, '_top_eqtl_snps_and_coefficients.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  currsums <- read.table(paste('/well/got2d/apayne/thesis/chapter6/lasso/eqtl_coef_sums_temp/chr_', j, '_top_eqtl_gene_model_summaries.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  print(paste('Chromosome ', j, ' completed.', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  coefs <- rbind(coefs, currcoef)
  sums <- rbind(sums, currsums)
}

write.table(coefs, col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t', file='/well/got2d/apayne/thesis/chapter6/eqtl_model_coefficients.txt')
write.table(sums, col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t', file='/well/got2d/apayne/thesis/chapter6/eqtl_model_summaries.txt')


#####################
# RUN LASSO 500 ITS #
#####################
library(parallel)
load('/well/got2d/apayne/thesis/chapter6/filtered.Rda')
snplocs <- snplocs[,-c(4:5)]
final_genes$seqid <- as.character(final_genes$seqid)
final_genes$seqid <- substring(final_genes$seqid, 4, nchar(final_genes$seqid))
modeller <- function(i) {
  gc()
  testexpr <- v[i,]
  gencodetemp <- final_genes[final_genes[,'gene_id']==i,]
  snpmin <- gencodetemp[,'start']
  snpmax <- gencodetemp[,'end']
  snps <- snplocs[snplocs[,'CHR']==as.integer(gencodetemp['seqid'])&((snplocs[,'POS']>(snpmin-1000000))&(snplocs[,'POS']<(snpmax+1000000))),'ID']
  testexpr <- resid(lm(testexpr ~ as.matrix(covs)))
  testgenotypes <- gt[,snps]
  gc()

  if ((ncol(t(t(testgenotypes)))>=1)&(length(testexpr)>0)) {

    if (ncol(t(t(testgenotypes)))==1) {
      testgenotypes <- matrix(c(rep(1, nrow(gt)), testgenotypes), ncol=2, byrow=FALSE)
      colnames(testgenotypes) <- c('dummy_snp', snps)
      rownames(testgenotypes) <- rownames(gt)
    }

    MSEs <- NULL
    fitr2 <- NULL
    cvfit <- cv.glmnet(testgenotypes, testexpr, nfolds=10, type.measure='mse', alpha=1, nlambda=100, keep=TRUE)
    MSEs <- cbind(MSEs, cvfit$cvm)
    ls <- cvfit$lambda
    fitr2 <- cbind(fitr2, unlist(lapply(1:length(ls), FUN=function(x) summary(lm(testexpr ~ cvfit$fit.preval[,x]))$r.squared)))
    counter <- 0

    while (counter<500) {
      print(paste(i, ' iteration ', counter, sep=''))
      counter <- counter + 1
      cvfit <- cv.glmnet(testgenotypes, testexpr, nfolds=10, type.measure='mse', alpha=1, nlambda=100, keep=TRUE)
      MSEs <- cbind(MSEs, cvfit$cvm[match(ls, cvfit$lambda)])
      fitr2 <- cbind(fitr2, unlist(lapply(match(ls, cvfit$lambda), FUN=function(x) if (is.na(x)) return(NA) else summary(lm(testexpr ~ cvfit$fit.preval[,x]))$r.squared)))
    }

    rownames(MSEs) <- ls
    rownames(fitr2) <- ls
    lambda.min <- as.numeric(names(which.min(rowMeans(MSEs, na.rm=TRUE))))
    return(list(cvfit, lambda.min, MSEs, fitr2))
  } else return(NA)

}

save(covs, v, final_genes, gt, snplocs, modeller, file='/well/got2d/apayne/thesis/chapter6/lasso/filtered.Rda')

basedir <- "/well/got2d/apayne/thesis/chapter6/lasso/" #Base dii of analysis
Rtemplate <- "templates/r_template.R" #Filename of analysis script (in base dir)
Subtemplate <- "templates/submit_template.sh" #Filename of submission script (in base dir)
Sourcetemplate <- "templates/r_source_template.R"
errordir <- "errors" #Name of error dir (in base dir)
outputdir <- "outputs" #Name of output dir (in base dir)
scriptdir <- "R_scripts" #Name of dir for final analysis scripts (in base dir)
scriptdirtemp <- 'R_scripts_temp'
submissiondir <- "submission_scripts" #Name of dir for submission scripts (in base dir)
submissiondirtemp <- 'submission_scripts_temp'
rechodir <- "R_echos"
rechodirtemp <- 'R_echos_temp'
Repstring <- "ASDF" #Name of string to be replaced
Chrstring <- 'JKL'

load('/well/got2d/apayne/thesis/chapter6/lasso/filtered.Rda')
snplocs_full <- snplocs
gt_full <- gt
final_genes_full <- final_genes
v_full <- v
rm(snplocs, gt, final_genes, v)
gc()

for (j in unique(final_genes_full$seqid)) {
  gc()
  snplocs <- subset(snplocs_full, CHR==j)
  gt <- gt_full[,match(snplocs$ID, colnames(gt_full))]
  final_genes <- subset(final_genes_full, seqid==j)
  v <- v_full[rownames(v_full)%in%final_genes$gene_id,]
  system2('mkdir', paste(basedir, 'chr', j, sep=''))
  system2('mkdir', paste(basedir, 'chr', j, '/', errordir, sep=''))
  system2('mkdir', paste(basedir, 'chr', j, '/', outputdir, sep=''))
  system2('mkdir', paste(basedir, 'chr', j, '/', scriptdir, sep=''))
  system2('mkdir', paste(basedir, 'chr', j, '/', scriptdirtemp, sep=''))
  system2('mkdir', paste(basedir, 'chr', j, '/', submissiondir, sep=''))
  system2('mkdir', paste(basedir, 'chr', j, '/', submissiondirtemp, sep=''))
  system2('mkdir', paste(basedir, 'chr', j, '/', rechodir, sep=''))
  system2('mkdir', paste(basedir, 'chr', j, '/', rechodirtemp, sep=''))
  system2('mkdir', paste(basedir, 'chr', j, '/final_data_files', sep=''))

  save(snplocs, gt, covs, modeller, v, final_genes, file=paste(basedir, 'chr', j, '/filtered.Rda', sep=''))

  for (i in 1:ceiling(nrow(v)/5)) {
    system2('sed', c(paste('s/', Chrstring, '/', j, '/g', sep=''), paste(basedir, Rtemplate, sep=""), '>', paste(basedir, 'chr', j, '/', scriptdirtemp, '/Rscript_', i, '.R', sep='')))
    system2('sed', c(paste('s/', Repstring, '/', i, '/g', sep=''), paste(basedir, 'chr', j, '/', scriptdirtemp, '/Rscript_', i, '.R', sep=''), '>', paste(basedir, 'chr', j, '/', scriptdir, '/Rscript_', i, '.R', sep='')))

    system2('sed', c(paste('s/', Chrstring, '/', j, '/g', sep=''), paste(basedir, Subtemplate, sep=""), '>', paste(basedir, 'chr', j, '/', submissiondirtemp, '/submissionscript_', i, '.R', sep='')))
    system2('sed', c(paste('s/', Repstring, '/', i, '/g', sep=''), paste(basedir, 'chr', j, '/', submissiondirtemp, '/submissionscript_', i, '.R', sep=''), '>', paste(basedir, 'chr', j, '/', submissiondir, '/submissionscript_', i, '.R', sep='')))

    system2('sed', c(paste('s/', Chrstring, '/', j, '/g', sep=''), paste(basedir, Sourcetemplate, sep=""), '>', paste(basedir, 'chr', j, '/', rechodirtemp, '/Recho_', i, '.R', sep='')))
    system2('sed', c(paste('s/', Repstring, '/', i, '/g', sep=''),  paste(basedir, 'chr', j, '/', rechodirtemp, '/Recho_', i, '.R', sep=''), '>', paste(basedir, 'chr', j, '/', rechodir, '/Recho_', i, '.R', sep='')))
  }

  system2('rm', c('-r', paste(basedir, 'chr', j, '/', scriptdirtemp, sep='')))
  system2('rm', c('-r', paste(basedir, 'chr', j, '/', rechodirtemp, sep='')))
  system2('rm', c('-r', paste(basedir, 'chr', j, '/', submissiondirtemp, sep='')))

  print(j)
}

for (j in as.character(1:22)) {
  submitdir <- paste('/well/got2d/apayne/thesis/chapter6/lasso/chr', j, '/submission_scripts/', sep='') #Location of submission files
  files <- list.files(submitdir)
  for (file in files) system2('qsub', paste(submitdir, file, sep=""))
}

##########################
#                        #
#      ##       ##       #
#      ##       ##       #
#                        #
#    ###         ###     #
#     ####     ####      #
#       #########        #
#                        #
#   RESTART WORK HERE    #
##########################

full_lasso_model <- list()

for (i in as.character(1:22)) {
  lassolist <- list.files(paste('/well/got2d/apayne/thesis/chapter6/lasso/chr', i, '/final_data_files', sep=''))

  for (j in lassolist) {
    load(paste('/well/got2d/apayne/thesis/chapter6/lasso/chr', i, '/final_data_files/', j, sep=''))
    full_lasso_model <- c(full_lasso_model, lasso_model)
    print(paste(i, ':', j, sep=''))
  }

}

lasso_model <- full_lasso_model
save(lasso_model, file='/well/got2d/apayne/thesis/chapter6/full_lasso_all_genes.Rda')


#####################################
# PROCESS DATA FROM 500 ITERATIONS  #
#####################################
library(glmnet)
library(parallel)
library(data.table)

load('/well/got2d/apayne/thesis/chapter6/full_lasso_all_genes.Rda')
load('/well/got2d/apayne/thesis/chapter6/lasso/filtered.Rda')
rm(final_genes, gt, snplocs)
gc()
covs <- as.matrix(covs)

modeller <- function(i) {
  gc()

  if (is.na(match(as.character(lasso_model[[i]][[2]]), as.character(lasso_model[[i]][[1]]$lambda)))) {
    coefs <- as.matrix(coef(lasso_model[[i]][[1]], s=lasso_model[[i]][[1]]$lambda.min))
    lamda_missing=1
  } else {
    lambda_missing=0
    coefs <- as.matrix(coef(lasso_model[[i]][[1]], s=lasso_model[[i]][[1]]$lambda[match(as.character(lasso_model[[i]][[2]]), as.character(lasso_model[[i]][[1]]$lambda))]))
  }

  coefs <- coefs[coefs[,1]!=0,]
  testexpr <- resid(lm(v[i,] ~ covs))

  if (length(coefs)%in%c(0,1)) {
    return(NA)
  } else {
    snps <- names(coefs)[-1]

    if (is.na(match(as.character(lasso_model[[i]][[2]]), as.character(lasso_model[[i]][[1]]$lambda)))) {
      res <- summary(lm(testexpr ~ lasso_model[[i]][[1]]$fit.preval[,match(lasso_model[[i]][[1]]$lambda.min, lasso_model[[i]][[1]]$lambda)]))
    } else {
      res <- summary(lm(testexpr ~ lasso_model[[i]][[1]]$fit.preval[,match(as.character(lasso_model[[i]][[2]]), as.character(lasso_model[[i]][[1]]$lambda))]))
    }
    r2 <- res$r.squared
    qval <- pf(res$fstatistic[1], res$fstatistic[2], res$fstatistic[3], lower.tail=F)

    return(list(data.frame(gene=rep(i, length(snps)), snp=snps, coefs=coefs[-1], stringsAsFactors=FALSE),
      data.frame(gene=i, nsnps=length(snps), preval.r2=r2, preval.p=qval)
    ))
  }

}

gc()
outframes <- mclapply(names(lasso_model), modeller, mc.cores=60)
navector <- unlist(lapply(1:length(outframes), FUN=function(x) !is.list(outframes[[x]])))
goodinds <- (1:length(outframes))[!navector]
cleanframes <- outframes[goodinds]
coefframe <- data.frame(rbindlist(mclapply(1:length(cleanframes), FUN=function(x) cleanframes[[x]][[1]], mc.cores=24)))
genesumframe <- data.frame(rbindlist(mclapply(1:length(cleanframes), FUN=function(x) cleanframes[[x]][[2]], mc.cores=24)))
write.table(coefframe, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/lasso/full_lasso_snps_and_coefficients.txt')
write.table(genesumframe, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/lasso/full_lasso_gene_model_summaries_temp.txt')

library(data.table)
library(parallel)

coefs <- fread('/well/got2d/apayne/thesis/chapter6/lasso/full_lasso_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
sums <- fread('/well/got2d/apayne/thesis/chapter6/lasso/full_lasso_gene_model_summaries_temp.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

load('/well/got2d/apayne/thesis/chapter6/lasso/filtered.Rda')
rm(final_genes, snplocs)
gt <- gt[,colnames(gt)%in%coefs$snp]
gc()
covs <- as.matrix(covs)

r2_calc <- function(i, whichdf, whichsums) {
  testexpr <- resid(lm(v[i,] ~ covs))
  currdf <- subset(whichdf, gene==i)
  snps <- currdf$snp
  xmat <- matrix(c(1, currdf[match(snps, currdf$snp),'coefs']), ncol=1)
  ymat <- as.matrix(cbind(data.frame(intcpt=rep(1, nrow(gt)), stringsAsFactors=FALSE), gt[,match(snps, colnames(gt))]))
  resvec <- ymat%*%xmat
  res <- summary(lm(testexpr ~ resvec))
  r2 <- res$r.squared

  if (is.null(res$fstatistic)) {
    qval <- 1
  } else {
    qval <- pf(res$fstatistic[1], res$fstatistic[2], res$fstatistic[3], lower.tail=F)
  }

  return(cbind(whichsums[match(i, whichsums$gene),], data.frame(fullinds.r2=r2, fullinds.p=qval, stringsAsFactors=FALSE)))
}

finalr2 <- rbindlist(mclapply(sums$gene, FUN=r2_calc, coefs, sums, mc.cores=60))

write.table(finalr2, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/lasso/full_lasso_gene_model_summaries.txt')


####################
# FILTER MV MODELS #
####################
library(parallel)
library(glmnet)
coef_table <- read.table('/well/got2d/apayne/thesis/chapter6/lasso/full_lasso_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
load('/well/got2d/apayne/thesis/chapter6/lasso/filtered.Rda')
gt <- gt[,colnames(gt)%in%coef_table$snp]
covs <- as.matrix(covs)

modeller <- function(i) {
  gc()
  testexpr <- resid(lm(v[i,] ~ covs))
  currcoefs <- subset(coef_table, gene==i)

  snps <- currcoefs$snp
  MSEs <- NULL
  testgenotypes <- gt[,snps]

  if (length(snps)>1) {
    fullmodel <- lm(testexpr ~ testgenotypes) 
    fullr2 <- summary(fullmodel)$r.squared
    r2per <- 0
    counter <- 0

    while (r2per<0.95) {
      counter <- counter + 1
      currsnps <- substr(rownames(summary(fullmodel)$coefficients)[-1], 14, nchar(rownames(summary(fullmodel)$coefficients)[-1]))[((summary(fullmodel)$coefficients)[-1,4])%in%(sort(summary(fullmodel)$coefficients[-1,4])[1:counter])]
      currmodel <- lm(testexpr ~ gt[,currsnps])
      r2per <- summary(currmodel)$r.squared/fullr2
    }

  } else {
    currsnps <- snps
  }

  testgenotypes <- gt[,currsnps]

  if (length(currsnps)==1) {
    testgenotypes <- matrix(c(rep(1, nrow(gt)), testgenotypes), ncol=2, byrow=FALSE)
    colnames(testgenotypes) <- c('dummy_snp', currsnps)
    rownames(testgenotypes) <- rownames(gt)
  }

  MSEs <- NULL
  cvfit <- cv.glmnet(testgenotypes, testexpr, nfolds=10, type.measure='mse', alpha=0, nlambda=100, keep=TRUE)
  MSEs <- cbind(MSEs, cvfit$cvm)
  ls <- cvfit$lambda
  counter2 <- 0

  while (counter2<100) {
    counter2 <- counter2 + 1
    cvfit <- cv.glmnet(testgenotypes, testexpr, nfolds=10, type.measure='mse', alpha=0, nlambda=100, keep=TRUE)
    MSEs <- cbind(MSEs, cvfit$cvm[match(ls, cvfit$lambda)])
  }

  rownames(MSEs) <- ls
  lambda.min <- as.numeric(names(which.min(rowMeans(MSEs, na.rm=TRUE))))

  write.table(c('a', 'a'), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter6/lasso/temp_iteration_counter/', i, sep=''))

  return(list(cvfit, lambda.min, currsnps))
}

its <- unique(coef_table$gene)
names(its) <- its
gc()

finalmodels <- mclapply(its, FUN=modeller, mc.cores=60)
save(finalmodels, file='/well/got2d/apayne/thesis/chapter6/lasso/finalmodels.Rda')

#################################################
# Extract coef and calculate r2 filtered models #
#################################################
library(glmnet)
library(parallel)
library(data.table)

load('/well/got2d/apayne/thesis/chapter6/lasso/finalmodels.Rda')
load('/well/got2d/apayne/thesis/chapter6/lasso/filtered.Rda')
rm(final_genes, gt, snplocs)
gc()
covs <- as.matrix(covs)

modeller <- function(i) {
  gc()

  if (is.na(match(as.character(finalmodels[[i]][[2]]), as.character(finalmodels[[i]][[1]]$lambda)))) {
    coefs <- as.matrix(coef(finalmodels[[i]][[1]], s=finalmodels[[i]][[1]]$lambda.min))
    lamda_missing=1
  } else {
    lambda_missing=0
    coefs <- as.matrix(coef(finalmodels[[i]][[1]], s=finalmodels[[i]][[1]]$lambda[match(as.character(finalmodels[[i]][[2]]), as.character(finalmodels[[i]][[1]]$lambda))]))
  }

  coefs <- coefs[coefs[,1]!=0,]
  testexpr <- resid(lm(v[i,] ~ covs))

  if (length(coefs)%in%c(0,1)) {
    return(NA)
  } else {
    snps <- names(coefs)[-1]

    if (is.na(match(as.character(finalmodels[[i]][[2]]), as.character(finalmodels[[i]][[1]]$lambda)))) {
      res <- summary(lm(testexpr ~ finalmodels[[i]][[1]]$fit.preval[,match(finalmodels[[i]][[1]]$lambda.min, finalmodels[[i]][[1]]$lambda)]))
    } else {
      res <- summary(lm(testexpr ~ finalmodels[[i]][[1]]$fit.preval[,match(as.character(finalmodels[[i]][[2]]), as.character(finalmodels[[i]][[1]]$lambda))]))
    }
    r2 <- res$r.squared
    qval <- pf(res$fstatistic[1], res$fstatistic[2], res$fstatistic[3], lower.tail=F)

    return(list(data.frame(gene=rep(i, length(snps)), snp=snps, coefs=coefs[-1], stringsAsFactors=FALSE),
      data.frame(gene=i, nsnps=length(snps), preval.r2=r2, preval.p=qval)
    ))
  }

}

outframes <- mclapply(names(finalmodels), modeller, mc.cores=60)
navector <- unlist(lapply(1:length(outframes), FUN=function(x) !is.list(outframes[[x]])))
goodinds <- (1:length(outframes))[!navector]
cleanframes <- outframes[goodinds]
coefframe <- data.frame(rbindlist(mclapply(1:length(cleanframes), FUN=function(x) cleanframes[[x]][[1]], mc.cores=24)))
genesumframe <- data.frame(rbindlist(mclapply(1:length(cleanframes), FUN=function(x) cleanframes[[x]][[2]], mc.cores=24)))
write.table(coefframe, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_snps_and_coefficients.txt')
write.table(genesumframe, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_gene_model_summaries_temp.txt')

library(parallel)
library(data.table)
coefs <- read.table('/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
gsums <- read.table('/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_gene_model_summaries_temp.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
load('/well/got2d/apayne/thesis/chapter6/lasso/filtered.Rda')
rm(final_genes, snplocs)
gt <- gt[,colnames(gt)%in%coefs$snp]
gc()
covs <- as.matrix(covs)

r2_calc <- function(i, whichdf, whichsums) {
  testexpr <- resid(lm(v[i,] ~ covs))
  currdf <- subset(whichdf, gene==i)
  snps <- currdf$snp
  xmat <- matrix(c(1, currdf[match(snps, currdf$snp),'coefs']), ncol=1)
  ymat <- as.matrix(cbind(data.frame(intcpt=rep(1, nrow(gt)), stringsAsFactors=FALSE), gt[,match(snps, colnames(gt))]))
  resvec <- ymat%*%xmat
  res <- summary(lm(testexpr ~ resvec))
  r2 <- res$r.squared

  if (is.null(res$fstatistic)) {
    qval <- 1
  } else {
    qval <- pf(res$fstatistic[1], res$fstatistic[2], res$fstatistic[3], lower.tail=F)
  }

  return(cbind(whichsums[match(i, whichsums$gene),], data.frame(fullinds.r2=r2, fullinds.p=qval, stringsAsFactors=FALSE)))
}

finalr2 <- rbindlist(mclapply(gsums$gene, FUN=r2_calc, coefs, gsums, mc.cores=62))

write.table(finalr2, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_gene_model_summaries.txt')


#######################
# LD OLD LASSO VS NEW #
#######################
library(parallel)
library(genetics)
lasso_coefs <- read.table('/well/got2d/apayne/thesis/chapter6/lasso/full_lasso_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
final_coefs <- read.table('/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
load('/well/got2d/apayne/thesis/chapter6/lasso/filtered.Rda')

gt <- gt[,colnames(gt)%in%lasso_coefs$snp]
rm(covs, final_genes, modeller, snplocs, v)
gc()

gt_matrix <- matrix(unlist(mclapply(1:ncol(gt), FUN=function(x) as.genotype.allele.count(round(gt[,x], 0), alleles=c('a','b')), mc.cores=56)), byrow=FALSE, ncol=ncol(gt))
colnames(gt_matrix) <- colnames(gt)

ld_checker_lasso <- function(gene_i) {
  lassosnps <- subset(lasso_coefs, gene==gene_i)$snp

  if (length(lassosnps)==1) return(NA)

  currgts <- data.frame(lapply(lassosnps, FUN=function(x) genotype(gt_matrix[,x])))
  names(currgts) <- lassosnps

  return((LD(currgts)$r^2))
}

its <- unique(lasso_coefs$gene)
names(its) <- its
ld_list_lasso <- mclapply(its, FUN=ld_checker_lasso, mc.cores=62)
names(ld_list_lasso) <- its


ld_checker_final <- function(gene_i) {
  finalsnps <- subset(final_coefs, gene==gene_i)$snp

  if (length(finalsnps)==1) return(NA)

  currgts <- data.frame(lapply(finalsnps, FUN=function(x) genotype(gt_matrix[,x])))
  names(currgts) <- finalsnps

  return((LD(currgts)$r^2))
}

its <- unique(final_coefs$gene)
names(its) <- its
ld_list_final <- mclapply(its, FUN=ld_checker_final, mc.cores=62)
names(ld_list_final) <- its

ldsums_lasso <- matrix(0, ncol=20, nrow=length(ld_list_lasso))
rownames(ldsums_lasso) <- names(ld_list_lasso)
colnames(ldsums_lasso) <- as.character(seq(0.05, 1, 0.05))

for (i in 1:length(ld_list_lasso)) {
  ldsums_lasso[names(ld_list_lasso)[i],] <- unlist(lapply(seq(0.05, 1, 0.05), FUN=function(x) sum(ld_list_lasso[[i]]>=x, na.rm=TRUE)))
}

ldwrite_lasso <- cbind(rownames(ldsums_lasso), ldsums_lasso)
colnames(ldwrite_lasso) <- c('gene_id', colnames(ldsums_lasso))
write.table(ldwrite_lasso, quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t', file='/well/got2d/apayne/thesis/chapter6/lasso/ld_r2_summaries_lasso_models.txt')

ldsums_final <- matrix(0, ncol=20, nrow=length(ld_list_final))
rownames(ldsums_final) <- names(ld_list_final)
colnames(ldsums_final) <- as.character(seq(0.05, 1, 0.05))

for (i in 1:length(ld_list_final)) {
  ldsums_final[names(ld_list_final)[i],] <- unlist(lapply(seq(0.05, 1, 0.05), FUN=function(x) sum(ld_list_final[[i]]>=x, na.rm=TRUE)))
}

ldwrite_final <- cbind(rownames(ldsums_final), ldsums_final)
colnames(ldwrite_final) <- c('gene_id', colnames(ldsums_final))
write.table(ldwrite_final, quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t', file='/well/got2d/apayne/thesis/chapter6/lasso/ld_r2_summaries_final_models.txt')



###################################
###################################
### ACTUAL THESIS ANALYSIS WORK ###
###################################
###################################
library(data.table)
library(ggplot2)
library(gridExtra)
matrix_results <- fread('/well/got2d/apayne/thesis/chapter6/fulleqtls.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
load('/well/got2d/apayne/thesis/chapter6/lasso/filtered.Rda')
gc()

totalgenes <- unique(matrix_results$gene)
egenes <- unique(subset(matrix_results, FDR_overall<=0.05)$gene)
length(unique(subset(matrix_results, FDR_overall<=0.01)$gene)) #compare to inspire
snps_per_gene <- table(matrix_results$gene)
gc()

dim(matrix_results)
median(snps_per_gene)
min(snps_per_gene)
max(snps_per_gene)

matrix_coefs <- read.table('/well/got2d/apayne/thesis/chapter6/eqtl_model_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
matrix_sums <- read.table('/well/got2d/apayne/thesis/chapter6/eqtl_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

matrix_sums_sig <- matrix_sums[matrix_sums$gene%in%egenes,]

median(matrix_sums_sig$fullinds.r2)
max(matrix_sums_sig$fullinds.r2)
min(matrix_sums_sig$fullinds.r2)





genesums <- read.table('/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
eqtlsums <- read.table('/well/got2d/apayne/thesis/chapter6/eqtl_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
genesums <- genesums[genesums$gene%in%eqtlsums$gene,]
eqtlsums <- eqtlsums[match(genesums$gene, eqtlsums$gene),]

plotframe <- as.data.frame(cbind(eqtlsums$fullinds.r2, genesums$fullinds.r2))
colnames(plotframe) <- c('eqtl', 'final')

pdf('/well/got2d/apayne/thesis/chapter6/compare_final_to_top_eqtl_r2.pdf', width=5.9, height=5.9)
ggplot(plotframe, aes(eqtl, final)) +
  geom_point() +
  xlim(0, round(max(plotframe), 1)) +
  ylim(0, round(max(plotframe), 1)) +
  geom_abline(slope=1, intercept=0, col='red') +
  labs(x=bquote('Top eQTL' ~ italic('R')^2),
    y=bquote('Final multi-variant model ' ~ italic('R')^2)
  )
dev.off()


####################
#  EQTL GWAS COLOC # pull from eqtl_gwas_overlap scripts from way back
####################
library(data.table)
library(parallel)
eqtl <- fread('/well/got2d/apayne/thesis/chapter6/fulleqtls.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
load('/well/got2d/apayne/thesis/chapter6/lasso/filtered.Rda')


for (currgwas in c('t2d', 'fi', 'fg')) {
gc()
for (thresh in c(5E-5, 1E-5, 5E-6, 1E-6, 5E-7, 1E-7, 5E-8)) {
  gwas_table <- fread(paste('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/', currgwas, '.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

  gwas_table <- subset(gwas_table, ID%in%eqtl$SNP)
  gwas_table$CHR <- eqtl[match(gwas_table$ID, eqtl$SNP),'CHR']
  gwas_table$POS <- eqtl[match(gwas_table$ID, eqtl$SNP),'snp_pos']
  gwas_sig <- subset(gwas_table, P<=thresh)

  gwas_pruned <- NULL

  while (nrow(gwas_sig)>0) {
    gwas_current <- subset(gwas_sig, CHR==gwas_sig[1,'CHR']&POS<(gwas_sig[1,'POS']+500000))
    gwas_sig <- gwas_sig[-(1:nrow(gwas_current)),]
    gwas_current <- gwas_current[match(min(gwas_current$P), gwas_current$P),]
    gwas_pruned <- rbind(gwas_pruned, gwas_current)
  }

  for (k in unique(gwas_pruned$CHR)) {
    numlds <- 1

    while (numlds>0) {
      gt_sigs <- t(t(gt[,colnames(gt)%in%subset(gwas_pruned, CHR==k)$ID]))

      if (ncol(gt_sigs)==1) {
        numlds <- 0
        next
      }

      gt_matrix <- matrix(unlist(lapply(1:ncol(gt_sigs), FUN=function(x) as.genotype.allele.count(round(gt_sigs[,x], 0), alleles=c('a','b')))), byrow=FALSE, ncol=ncol(gt_sigs))
      colnames(gt_matrix) <- colnames(gt_sigs)
      currgts <- data.frame(lapply(colnames(gt_sigs), FUN=function(x) genotype(gt_matrix[,x])))
      names(currgts) <- colnames(gt_sigs)
      currld <- LD(currgts)$r^2

      if (sum(currld>0.8, na.rm=TRUE)>0) {
        currpair <- rownames(currld)[which(currld==max(currld, na.rm=TRUE), arr.ind=TRUE)]
        highp <- currpair[match(max(gwas_pruned[match(currpair, gwas_pruned$ID),'P']), gwas_pruned[match(currpair, gwas_pruned$ID),'P'])]
        gwas_pruned <- gwas_pruned[-match(highp, gwas_pruned$ID),]
        print(paste('k', ': ', currpair[1], ' and ', currpair[2], sep=''))
      } else {
        numlds <- 0
      }

    }

  }

  eqtlmatch <- function(i) {
    curreqtls <- subset(eqtl, SNP==gwas_pruned[i,'ID'])
    curreqtls <- curreqtls[match(min(curreqtls$"p-value"), curreqtls$"p-value"),]
    return(curreqtls)
  }

  paired_eqtls <- do.call('rbind', mclapply(1:nrow(gwas_pruned), FUN=eqtlmatch, mc.cores=32))
  eqtls_sig <- subset(paired_eqtls, FDR_overall<=0.05)
  eqtls_sig <- cbind(eqtls_sig, gwas_pruned[match(eqtls_sig$SNP, gwas_pruned$ID),c('EFFECT', 'SE', 'P')])
  write.table(eqtls_sig, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter6/eqtl_gwas_coloc/', currgwas, '.', as.character(thresh), '.txt', sep=''))
  write.table(gwas_pruned, file=paste('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/lead_variants/', currgwas, '.', as.character(thresh), '.lead.variants', sep=''), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
  print(paste(currgwas, ' completed at ', as.character(thresh), sep=''))
}
}

gwas_sig_number <- data.frame(Threshold=as.character(c(5E-5, 1E-5, 5E-6, 1E-6, 5E-7, 1E-7, 5E-8)),
  DIAGRAM=unlist(lapply(c(5E-5, 1E-5, 5E-6, 1E-6, 5E-7, 1E-7, 5E-8), function(thresh) as.integer(strsplit(system2('wc', c('-l', paste('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/lead_variants/t2d.', as.character(thresh), '.lead.variants', sep='')), stdout=TRUE), ' ', fixed=TRUE)[[1]][1])-1)),
  FG=unlist(lapply(c(5E-5, 1E-5, 5E-6, 1E-6, 5E-7, 1E-7, 5E-8), function(thresh) as.integer(strsplit(system2('wc', c('-l', paste('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/lead_variants/fg.', as.character(thresh), '.lead.variants', sep='')), stdout=TRUE), ' ', fixed=TRUE)[[1]][1])-1)),
  FI=unlist(lapply(c(5E-5, 1E-5, 5E-6, 1E-6, 5E-7, 1E-7, 5E-8), function(thresh) as.integer(strsplit(system2('wc', c('-l', paste('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/lead_variants/fi.', as.character(thresh), '.lead.variants', sep='')), stdout=TRUE), ' ', fixed=TRUE)[[1]][1])-1)),
  stringsAsFactors=FALSE
)

write.table(gwas_sig_number, sep='\t&\t', quote=FALSE, col.names=FALSE, row.names=FALSE, eol='\\\\\n\\hline\n', file='/well/got2d/apayne/thesis/chapter6/significant_variants_per_gwas.latexformat')


#TEST for unique signals
t2d <- read.table('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/lead_variants/t2d.5e-08.lead.variants', sep='\t', header=TRUE, stringsAsFactors=FALSE)
fi <- read.table('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/lead_variants/fi.5e-08.lead.variants', sep='\t', header=TRUE, stringsAsFactors=FALSE)
fg <- read.table('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/lead_variants/fg.5e-08.lead.variants', sep='\t', header=TRUE, stringsAsFactors=FALSE)
gwas_eqtls <- subset(eqtl, SNP%in%c(t2d$ID, fg$ID, fi$ID))
dim(gwas_eqtls)
length(unique(gwas_eqtls$gene))

t2d_coloc <- subset(gwas_eqtls, SNP%in%t2d$ID)
fg_coloc <- subset(gwas_eqtls, SNP%in%fg$ID)
fi_coloc <- subset(gwas_eqtls, SNP%in%fi$ID)
t2d_coloc <- subset(t2d_coloc, FDR<=0.05)
fg_coloc <- subset(fg_coloc, FDR<=0.05)
fi_coloc <- subset(fi_coloc, FDR<=0.05)
t2d_gwas <- read.table('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/t2d.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
fg_gwas <- read.table('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/fg.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
fi_gwas <- read.table('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/fi.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
t2d_coloc$P <- t2d_gwas[match(t2d_coloc$SNP, t2d_gwas$ID),'P']
fg_coloc$P <- fg_gwas[match(fg_coloc$SNP, fg_gwas$ID),'P']
fi_coloc$P <- fi_gwas[match(fi_coloc$SNP, fi_gwas$ID),'P']
write.table(t2d_coloc, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/eqtl_gwas_coloc/full_eqtls_not_top_only/t2d.txt')
write.table(fg_coloc, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/eqtl_gwas_coloc/full_eqtls_not_top_only/fg.txt')
write.table(fi_coloc, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/eqtl_gwas_coloc/full_eqtls_not_top_only/fi.txt')

v1 <- 'rs7944584'
v2 <- 'rs1483121'
gt_matrix <- matrix(unlist(lapply(c(v1, v2), FUN=function(x) as.genotype.allele.count(round(gt[,x], 0), alleles=c('a','b')))), byrow=FALSE, ncol=2)
colnames(gt_matrix) <- c(v1, v2)
currgts <- data.frame(lapply(colnames(gt_matrix), FUN=function(x) genotype(gt_matrix[,x])))
names(currgts) <- colnames(gt_matrix)
LD(currgts)$r^2

overlaps <- rbind(overlaps, data.frame(fi=NA, fg='rs2191349', t2d='rs2215383', stringsAsFactors=FALSE))

t2d <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
fg <- c(1, 1, 1, 1, 1, 1)
fi <- c(1)
t2dfg <- c(1, 1, 1)
t2dfi <- c()
fgfi <- c(1)
all <- c()

overlaps <- NULL
overlaps <- rbind(overlaps, data.frame(fi='rs780094', fg='rs780093', t2d=NA, stringsAsFactors=FALSE))
overlaps <- rbind(overlaps, data.frame(fi='rs9987289', fg='rs9987289', t2d=NA, stringsAsFactors=FALSE))
overlaps <- rbind(overlaps, data.frame(fi=NA, fg='rs11708067', t2d='rs11708067', stringsAsFactors=FALSE))
overlaps <- rbind(overlaps, data.frame(fi=NA, fg='rs2191349', t2d='rs2215383', stringsAsFactors=FALSE))
overlaps <- rbind(overlaps, data.frame(fi=NA, fg='rs11558471', t2d='rs3802177', stringsAsFactors=FALSE))
overlaps <- rbind(overlaps, data.frame(fi=NA, fg='rs11603334', t2d='rs11603334', stringsAsFactors=FALSE))

#Results for coloc
fullrestable <- NULL

currtab <- read.table('/well/got2d/apayne/thesis/chapter6/eqtl_gwas_coloc/full_eqtls_not_top_only/t2d.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
currtab$gwas <- 'T2D'
fullrestable <- rbind(fullrestable, currtab[,c('SNP', 'CHR', 'snp_pos', 'gwas', 'P', 'gene_name', 'p.value')])
currtab <- read.table('/well/got2d/apayne/thesis/chapter6/eqtl_gwas_coloc/full_eqtls_not_top_only/fi.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
currtab$gwas <- 'FI'
fullrestable <- rbind(fullrestable, currtab[,c('SNP', 'CHR', 'snp_pos', 'gwas', 'P', 'gene_name', 'p.value')])
currtab <- read.table('/well/got2d/apayne/thesis/chapter6/eqtl_gwas_coloc/full_eqtls_not_top_only/fg.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
currtab$gwas <- 'FG'
fullrestable <- rbind(fullrestable, currtab[,c('SNP', 'CHR', 'snp_pos', 'gwas', 'P', 'gene_name', 'p.value')])
fullrestable$chrpos <- paste(as.character(fullrestable$CHR), as.character(fullrestable$snp_pos), sep=':')
write.table(fullrestable, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/eqtl_gwas_coloc/full_significant_results.txt')


#######################
# DRAFT 2 LD checking #
#######################
library(data.table)
library(parallel)
library(genetics)
eqtl <- fread('/well/got2d/apayne/thesis/chapter6/fulleqtls.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
load('/well/got2d/apayne/thesis/chapter6/lasso/filtered.Rda')
fullrestable <- read.table('/well/got2d/apayne/thesis/chapter6/eqtl_gwas_coloc/full_significant_results.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
fullrestable$chrpos <- paste(as.character(fullrestable$CHR), as.character(fullrestable$snp_pos), sep=':')

eqtl <- subset(eqtl, gene_name%in%fullrestable$gene_name)
gt <- gt[,colnames(gt)%in%eqtl$SNP]
gc()

ld_checker <- function(i) {
  currtest <- fullrestable[i,]
  curreqtl <- subset(eqtl, gene_name==currtest$gene_name&(snp_pos>(currtest$snp_pos-11000)|snp_pos<(currtest$snp_pos+11000)))
  topeqtl <- curreqtl[1,]

  if (topeqtl$SNP==currtest$SNP) {
    outmat <- matrix(c(NA, NA, 1, NA), byrow=FALSE, nrow=2)
    colnames(outmat) <- rownames(outmat) <- rep(topeqtl$SNP, 2)
  }

  gt_matrix <- matrix(unlist(lapply(match(c(currtest$SNP, topeqtl$SNP), colnames(gt)), FUN=function(x) as.genotype.allele.count(round(gt[,x], 0), alleles=c('a', 'b')))), byrow=FALSE, ncol=2)
  colnames(gt_matrix) <- c(currtest$SNP, topeqtl$SNP)
  currgts <- data.frame(lapply(colnames(gt_matrix), FUN=function(x) genotype(gt_matrix[,x])))
  names(currgts) <- colnames(gt_matrix)

  return(LD(currgts)$r^2)
}

fullrestable$r2 <-  unlist(lapply(1:nrow(fullrestable), function(x) lapply(x, ld_checker)[[1]][1,2]))

write.table(fullrestable, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/eqtl_gwas_coloc/results_after_ld/full_significant_results_ld_with_top_eqtl.txt')

fullrestable <- subset(fullrestable, r2>=0.8)

fullrestable <- fullrestable[,c('SNP', 'chrpos', 'gwas', 'P', 'gene_name', 'p.value')]

fullrestable$P <- signif(fullrestable$P, 3)
fullrestable$p.value <- signif(fullrestable$p.value, 3)
colnames(fullrestable) <- c('SNP', 'chrpos', 'Trait', 'Trait P', 'eGene', 'eQTL P')

plottable <- fullrestable
colnames(plottable) <- c('SNP', 'chrpos', 'trait', 'traitp', 'egene', 'eqtlp')

publishedloci <- read.table('/well/got2d/apayne/thesis/chapter6/published_eqtl_gwas_coloc.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
publishedloci <- publishedloci[,-6]
publishedloci <- publishedloci[rowSums(publishedloci[,-1])>0,]
publishedloci <- publishedloci[order(publishedloci$ottosson, publishedloci$taneera, publishedloci$fadista, publishedloci$vandebunt, publishedloci$varshney, decreasing=TRUE),]
publishedloci <- publishedloci[,c(1, 3, 4, 2, 5, 6)]
publishedloci[publishedloci==1] <- 'X'
publishedloci[publishedloci==0] <- ''
#publishedloci[,1] <- paste(publishedloci[,1], '}', sep='')
#write.table(publishedloci, sep='\t&\t', eol='\\\\\n\\hline\n\\emph{', col.names=FALSE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/published_eqtl_gwas_coloc.latex')

plottable$egene[plottable$egene%in%publishedloci$gene_name]
length(unique(plottable$egene))
nrow(fullrestable[match(unique(fullrestable$SNP), fullrestable$SNP),])

T2D <- 'rs13094957', 'rs4402960', 'rs11257655', 'rs2292626', 'rs7903146', 'rs1061810', 'rs7178572'
FG <- 'rs10830963', 'rs7944584', 'rs174541'
FI <- 'rs6963105'
T2D.FG <- 'rs11708067', 'rs2215383', 'rs11603334', 

length(unique(plottable[plottable$published,]$egene))
unique(plottable[plottable$published,]$egene)
unique(plottable[!plottable$published,]$egene)
plottable[!plottable$published,]

subset(fullrestable, SNP%in%names(table(fullrestable$SNP)[table(fullrestable$SNP)>1]))
plottable$published <- plottable$egene%in%publishedloci$gene_name

inspire <- c('ADCY5', 'CAMK1D', 'STARD10', 'DGKB', 'HMG20A', 'IGF2BP2', 'KLH42', 'ITGB6', 'TCF7L2', 'UBE2E2', 'PDE8B', 'DGKB', 'FADS1', 'CTNNAL1', 'TCF7L2')
sum(plottable$egene%in%inspire)

#pdf('/well/got2d/apayne/thesis/chapter6/eqtl_coloc_gwas_eqtl_p.pdf', width=5.9, height=5.9)
pdf('/well/got2d/apayne/thesis/chapter6/eqtl_coloc_gwas_eqtl_p_ld.pdf', width=5.9, height=5.9)
ggplot(plottable, aes(x=-log10(traitp), y=-log10(eqtlp), label=egene)) +
  geom_point(aes(col=published, shape=as.factor(trait))) +
  geom_text_repel(aes(label=egene), size=2.5) +
  labs(x=expression(paste('-log'['10']*'(MA ', italic('p'), ')', sep='')),
    y=expression(paste('-log'['10']*'(eQTL ', italic('p'), ')', sep=''))
  ) +
  scale_colour_manual(values=c('blue', 'red'), labels=c('Previously unpublished', 'Published')) +
  guides(colour=guide_legend(title=''), shape=guide_legend(title='Trait')) +
  theme(legend.position='top')
dev.off()

overlaps <- NULL
#overlaps <- rbind(overlaps, data.frame(fi='rs9987289', fg='rs9987289', t2d=NA, stringsAsFactors=FALSE))
#overlaps <- rbind(overlaps, data.frame(fi=NA, fg='rs11708067', t2d='rs11708067', stringsAsFactors=FALSE))
overlaps <- rbind(overlaps, data.frame(fi=NA, fg='rs2191349', t2d='rs2215383', stringsAsFactors=FALSE))
#overlaps <- rbind(overlaps, data.frame(fi=NA, fg='rs11603334', t2d='rs11603334', stringsAsFactors=FALSE))

table_write <- fullrestable
table_write$eGene <- paste('\\emph{', table_write$eGene, '}', sep='')

write.table(table_write, col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t&\t', eol='\\\\\n\\hline\n', file='/well/got2d/apayne/thesis/chapter6/gwas_eqtl_coloc_5e-08_ld.txt.latexformat')

library(rtracklayer)
gencode <- readGFF('/well/got2d/rna-seq/resources/gencode.v19.annotation.gtf')
plottable$gene_type <- gencode[match(plottable$egene, gencode$gene_name),'gene_type']

#Compare STAG3L1 and PMS2P3 snps from Varshney and here
v1 <- 'rs1167800'
v2 <- 'rs1483121'
gt_matrix <- matrix(unlist(lapply(c(v1, v2), FUN=function(x) as.genotype.allele.count(round(gt[,x], 0), alleles=c('a','b')))), byrow=FALSE, ncol=2)
colnames(gt_matrix) <- c(v1, v2)
currgts <- data.frame(lapply(colnames(gt_matrix), FUN=function(x) genotype(gt_matrix[,x])))
names(currgts) <- colnames(gt_matrix)
LD(currgts)$r^2


##################
# Plot SV colocs #
##################
library(data.table)
library(parallel)
library(genetics)
library(ggplot2)

eqtl <- fread('/well/got2d/apayne/thesis/chapter6/fulleqtls.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
load('/well/got2d/apayne/thesis/chapter6/filtered.Rda')

for (currgwas in c('t2d', 'fg', 'fi')) {
  currtab <- subset(read.table('/well/got2d/apayne/thesis/chapter6/eqtl_gwas_coloc/results_after_ld/full_significant_results_ld_with_top_eqtl.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE), r2>=0.8&gwas==toupper(currgwas))
  currsums <- fread(paste('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/', currgwas, '.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

  for (currind in 1:nrow(currtab)) {
#  tophit <- currtab[match(max(log(currtab$p.value)*log(currtab$P)), log(currtab$p.value)*log(currtab$P)),]
    tophit <- currtab[currind,]
    currgene_name <- tophit$gene_name
    snpset <- subset(eqtl, gene_name==tophit$gene_name)
    snpset$gwasp <- currsums[match(snpset$SNP, currsums$ID),'P']

    maxeqtl <- max(-log10(snpset$"p-value"))
    mingwas <- min(log10(snpset$gwasp))

    snpset$scaledeqtl <- -log10(snpset$"p-value")/maxeqtl
    snpset$scaledgwas <- -log10(snpset$gwasp)/mingwas
    snpset$istop <- 'no'
    snpset[snpset$SNP==tophit$SNP,'istop'] <- 'yes'
    snpset$istop <- factor(snpset$istop, levels=c('yes', 'no'))

    currgt <- gt[,match(snpset$SNP, colnames(gt))]
    gt_matrix <- matrix(unlist(mclapply(1:ncol(currgt), FUN=function(x) as.genotype.allele.count(round(currgt[,x], 0), alleles=c('a','b')), mc.cores=32)), byrow=FALSE, ncol=ncol(currgt))
    colnames(gt_matrix) <- colnames(currgt)

    currgts <- data.frame(mclapply(colnames(gt_matrix), FUN=function(x) genotype(gt_matrix[,x]), mc.cores=32))
    names(currgts) <- colnames(gt_matrix)

    lds <- unlist(mclapply(names(currgts)[names(currgts)!=tophit$SNP], FUN=function(x) (LD(currgts[c(tophit$SNP, x)])$r^2)[1,2], mc.cores=32))
    lds <- c(1, lds)
    names(lds) <- c(tophit$SNP, names(currgts)[names(currgts)!=tophit$SNP])

    snpset$ldwithtop <- lds[match(snpset$SNP, names(lds))]
    snpset$snp_pos_mb <- snpset$snp_pos/1000000

    pdf(paste('/well/got2d/apayne/thesis/chapter6/eqtl_gwas_coloc/plots_after_ld_filtering/', currgwas, '_', currgene_name, '_', as.character(currind), '_eqtl_coloc_miami.pdf', sep=''), width=5.9, height=5.9)
    print(ggplot(snpset) +
    geom_point(aes(x=snp_pos_mb,
        y=scaledeqtl,
        colour=ldwithtop,
        shape=istop
      )
    ) +
    geom_point(
      aes(x=snp_pos_mb,
        y=scaledgwas,
        colour=ldwithtop,
        shape=istop
      )
    ) +
    geom_hline(yintercept=0,
      col='black',
      size=1.5
    ) +
    geom_segment(
      aes(x=snpset[1,'gene_start']/1000000,
        y=0,
        xend=snpset[1,'gene_end']/1000000,
        yend=0,
        size=1
      ),
      col='darkorange'
    ) +
    scale_colour_gradient(low='black',
      high='red',
      breaks=c(0, 0.5, 1),
      labels=c('0', '0.5', '1'),
      limits=c(0, 1)
    ) +
    scale_shape_manual(values=c(17, 16),
      guide=guide_legend(override.aes=list(
          shape=c(15, 17),
          colour=c('darkorange', 'red')
        ),
        order=1
      ),
      labels=c(bquote(~italic(.(tophit$gene_name))), expression(paste('Top MA SNP in region      ', bold(italic(R))^2, ' with top SNP:', sep='')))
    ) +
    guides(size=FALSE,
      colour=guide_colorbar(order=0)
    ) +
    scale_y_continuous(breaks=c(floor(abs(mingwas))/mingwas, 0, floor(maxeqtl)/maxeqtl),
      labels=c(as.character(floor(abs(mingwas))), '0', as.character(floor(maxeqtl)))
    ) +
    labs(x=paste('Genomic position on chromosome ', as.character(snpset[1,'CHR']), ' (Mb)', sep=''),
#      y=expression(paste(italic(paste(toupper(currgwas), ' MA', sep='')), '                    -log10(p)                    ', italic('eQTL'), sep='')),
       y=bquote(~italic(paste(.(toupper(currgwas)), ' MA', sep='')) ~ '                 -log'['10'] ~ '('*italic('p')*')                      ' ~ italic('eQTL')),
      shape='',
      colour=''
    ) +
    theme(legend.position='top', legend.spacing=unit(0, 'cm'))
    )
    dev.off()
    }
  }


#### LD FOR MADD
gt_sigs <- gt[,c('rs7944584', 'rs1483121')]
gt_matrix <- matrix(unlist(lapply(1:ncol(gt_sigs), FUN=function(x) as.genotype.allele.count(round(gt_sigs[,x], 0), alleles=c('a','b')))), byrow=FALSE, ncol=ncol(gt_sigs))
colnames(gt_matrix) <- colnames(gt_sigs)
currgts <- data.frame(lapply(colnames(gt_sigs), FUN=function(x) genotype(gt_matrix[,x])))
names(currgts) <- colnames(gt_sigs)
currld <- LD(currgts)$r^2

t2d_de_islets <- read.table('/well/got2d/apayne/islet_t2d_networks/de_by_t2d.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
t2d_de_islets[t2d_de_islets$gene_name=='MADD',]


########################
# RUN MR AND HET TESTS #
########################
library(parallel)
library(MendelianRandomization)
library(genetics)
load('/well/got2d/apayne/thesis/chapter6/lasso/filtered.Rda')
genesums <- read.table('/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
genecoefs <- read.table('/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
gt <- gt[,colnames(gt)%in%genecoefs$snp]
rm(covs, final_genes, modeller, snplocs, v)
gc()

mr_runs <- function(gwas, currgene) {
  snps <- subset(genecoefs, gene==currgene)$snp

  currgt <- t(t(gt[,match(snps, colnames(gt))]))
  gt_matrix <- matrix(unlist(lapply(1:ncol(currgt), FUN=function(x) as.genotype.allele.count(round(currgt[,x], 0), alleles=c('a','b')))), byrow=FALSE, ncol=ncol(currgt))
  colnames(gt_matrix) <- colnames(currgt)

  if (length(snps)>1) {
    currgts <- data.frame(lapply(snps, FUN=function(x) genotype(gt_matrix[,x])))
    names(currgts) <- snps
    lds <- LD(currgts)$r
    lds[lower.tri(lds)] <- t(lds)[lower.tri(t(lds))]
    diag(lds) <- rep(1, length(snps))
  }

  mr_ob_nocor <- mr_input(bx=subset(genecoefs, gene==currgene)$coefs, by=gwas[match(snps, gwas$ID),'EFFECT'], byse=gwas[match(snps, gwas$ID),'SE'], exposure=paste(currgene, '_expression', sep=''), outcome='GWAS', snps=snps)

  if (length(snps)>2) {
    mr_ob_cor <- mr_input(bx=subset(genecoefs, gene==currgene)$coefs, by=gwas[match(snps, gwas$ID),'EFFECT'], byse=gwas[match(snps, gwas$ID),'SE'], exposure=paste(currgene, '_expression', sep=''), outcome='GWAS', snps=snps, correlation=lds)
    mr_res_nocor_ivw <- tryCatch(mr_ivw(mr_ob_nocor),
      error=function(x) {
        return(list(Estimate=as.numeric(NA), StdError=as.numeric(NA), Heter.Stat=c(as.numeric(NA), as.numeric(NA))))
      },
      warning=function(x) {
        return(list(Estimate=as.numeric(NA), StdError=as.numeric(NA), Heter.Stat=c(as.numeric(NA), as.numeric(NA))))
      }
    )
    mr_res_nocor_egger <- tryCatch(mr_egger(mr_ob_nocor),
      error=function(x) {
        return(data.frame(Intercept=as.numeric(NA), StdError.Int=as.numeric(NA), stringsAsFactors=FALSE))
      },
      warning=function(x) {
        return(data.frame(Intercept=as.numeric(NA), StdError.Int=as.numeric(NA), stringsAsFactors=FALSE))
      }
    )
    mr_res_cor_ivw <- tryCatch(mr_ivw(mr_ob_cor),
      error=function(x) {
        return(list(Estimate=as.numeric(NA), StdError=as.numeric(NA), Heter.Stat=c(as.numeric(NA), as.numeric(NA))))
      },
      warning=function(x) {
        return(list(Estimate=as.numeric(NA), StdError=as.numeric(NA), Heter.Stat=c(as.numeric(NA), as.numeric(NA))))
      }
    )
    mr_res_cor_egger <- tryCatch(mr_egger(mr_ob_cor),
      error=function(x) {
        return(data.frame(Intercept=as.numeric(NA), StdError.Int=as.numeric(NA), stringsAsFactors=FALSE))
      },
      warning=function(x) {
        return(data.frame(Intercept=as.numeric(NA), StdError.Int=as.numeric(NA), stringsAsFactors=FALSE))
      }
    )
    rframenocor <- data.frame(
    gene=currgene,
    nsnps=length(snps),
    fullinds.r2=subset(genesums, gene==currgene)$fullinds.r2,
    fullinds.p=subset(genesums, gene==currgene)$fullinds.p,
    ivw_est=mr_res_nocor_ivw$Estimate,
    ivw_se=mr_res_nocor_ivw$StdError,
    het_stat=mr_res_nocor_ivw$Heter.Stat[1],
    het_p=mr_res_nocor_ivw$Heter.Stat[2],
    egger_int=mr_res_nocor_egger$Intercept,
    egger_int_se=mr_res_nocor_egger$StdError.Int,
    stringsAsFactors=FALSE)
    rframecor <- data.frame(
    gene=currgene,
    nsnps=length(snps),
    fullinds.r2=subset(genesums, gene==currgene)$fullinds.r2,
    fullinds.p=subset(genesums, gene==currgene)$fullinds.p,
    ivw_est=mr_res_cor_ivw$Estimate,
    ivw_se=mr_res_cor_ivw$StdError,
    het_stat=mr_res_cor_ivw$Heter.Stat[1],
    het_p=mr_res_cor_ivw$Heter.Stat[2],
    egger_int=mr_res_cor_egger$Intercept,
    egger_int_se=mr_res_cor_egger$StdError.Int,
    stringsAsFactors=FALSE)
  } else {
    if (length(snps)==1) {
      mr_ob_cor <- mr_ob_nocor
      ivw_res_nocor <- ivw_res_cor <- mr_ivw(mr_ob_nocor)
    } else {
      mr_ob_cor <- mr_input(bx=subset(genecoefs, gene==currgene)$coefs, by=gwas[match(snps, gwas$ID),'EFFECT'], byse=gwas[match(snps, gwas$ID),'SE'], exposure=paste(currgene, '_expression', sep=''), outcome='GWAS', snps=snps, correlation=lds)
      ivw_res_nocor <- mr_ivw(mr_ob_nocor)
      ivw_res_cor <- mr_ivw(mr_ob_cor)
    }
    rframenocor <- data.frame(
    gene=currgene,
    nsnps=length(snps),
    fullinds.r2=subset(genesums, gene==currgene)$fullinds.r2,
    fullinds.p=subset(genesums, gene==currgene)$fullinds.p,
    ivw_est=ivw_res_nocor$Estimate,
    ivw_se=ivw_res_nocor$StdError,
    het_stat=ivw_res_nocor$Heter.Stat[1],
    het_p=ivw_res_nocor$Heter.Stat[2],
    egger_int=NA,
    egger_int_se=NA,
    stringsAsFactors=FALSE)
    rframecor <- data.frame(
    gene=currgene,
    nsnps=length(snps),
    fullinds.r2=subset(genesums, gene==currgene)$fullinds.r2,
    fullinds.p=subset(genesums, gene==currgene)$fullinds.p,
    ivw_est=ivw_res_cor$Estimate,
    ivw_se=ivw_res_cor$StdError,
    het_stat=ivw_res_cor$Heter.Stat[1],
    het_p=ivw_res_cor$Heter.Stat[2],
    egger_int=NA,
    egger_int_se=NA,
    stringsAsFactors=FALSE)
  }

  return(list(mr_ob_nocor, mr_ob_cor, rframenocor, rframecor))
}

for (i in c('t2d', 'fi', 'fg')) {
  gwas <- read.table(paste('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/', i, '.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  gc()
  mr_results <- mclapply(genesums$gene, FUN=mr_runs, gwas=gwas, mc.cores=56)
  names(mr_results) <- genesums$gene
  save(mr_results, file=paste('/well/got2d/apayne/thesis/chapter6/mr_results/', i, '_mr_results.Rda', sep=''))
  mr_nocor <- do.call('rbind', mclapply(genesums$gene, FUN=function(x) mr_results[[x]][[3]], mc.cores=56))
  mr_cor <- do.call('rbind', mclapply(genesums$gene, FUN=function(x) mr_results[[x]][[4]], mc.cores=56))
  write.table(mr_nocor, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter6/mr_results/', i, '_no_correlation_included.txt', sep=''))
  write.table(mr_cor, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter6/mr_results/', i, '_correlation_included.txt', sep=''))
}


################
# RUN METAXCAN #
################
library(parallel)
load('/well/got2d/apayne/thesis/chapter6/lasso/filtered.Rda')
genesums <- read.table('/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
genecoefs <- read.table('/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
gt <- gt[,colnames(gt)%in%genecoefs$snp]
rm(covs, final_genes, modeller, snplocs, v)
gc()

metaxcan <- function(currgene, gwas) {
  if (currgene%in%gcgenes) gc()
  currcoefs <- subset(genecoefs, gene==currgene)
  covs <- as.matrix(var(gt[,match(currcoefs$snp, colnames(gt))]))

  w <- currcoefs$coefs
  b <- gwas[match(currcoefs$snp, gwas$ID),'EFFECT']
  bse <- gwas[match(currcoefs$snp, gwas$ID),'SE']
  sl <- sqrt(diag(covs))
  sg <- sqrt(t(w)%*%covs%*%t(t(w)))

  metaz <- sum(w*b*sl/(sg*bse))

  return(data.frame(gene=currgene,
    nsnps=genesums[genesums$gene==currgene,'nsnps'],
    fullinds.r2=genesums[genesums$gene==currgene,'fullinds.r2'],
    fullinds.p=genesums[genesums$gene==currgene,'fullinds.p'],
    metaxcan.z=metaz,
    metaxcan.p=2*pnorm(-abs(metaz)),
    stringsAsFactors=FALSE)
  )
}

gcgenes <- sample(genesums$gene, 1800)

for (currgwas in c('fi', 'fg')) {
  gwas <- read.table(paste('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/', currgwas, '.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE)
  metares <- do.call('rbind', mclapply(genesums$gene, FUN=metaxcan, gwas, mc.cores=56))
  write.table(metares, col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t', file=paste('/well/got2d/apayne/thesis/chapter6/metaxcan_results/', currgwas, '.metaxcan.res', sep=''))
  print(currgwas)
}


###########################
# ANALYSE MR AND METAXCAN #
###########################
library(rtracklayer)
library(ggplot2)
gencode <- readGFF('/well/got2d/rna-seq/resources/gencode.v19.annotation.gtf')
t2dres <- read.table('/well/got2d/apayne/thesis/chapter6/mr_results/t2d_no_correlation_included.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
fgres <- read.table('/well/got2d/apayne/thesis/chapter6/mr_results/fg_no_correlation_included.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
fires <- read.table('/well/got2d/apayne/thesis/chapter6/mr_results/fi_no_correlation_included.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

t2dres$ivw_p <- 2*pnorm(-abs(t2dres$ivw_est/t2dres$ivw_se))
t2dres$egger_p <- 2*pnorm(-abs(t2dres$egger_int/t2dres$egger_int_se))
t2dres$ivw_fdr <- p.adjust(t2dres$ivw_p, method='BH')
t2dres$egger_fdr <- p.adjust(t2dres$egger_p, method='BH')
fgres$ivw_p <- 2*pnorm(-abs(fgres$ivw_est/fgres$ivw_se))
fgres$egger_p <- 2*pnorm(-abs(fgres$egger_int/fgres$egger_int_se))
fgres$ivw_fdr <- p.adjust(fgres$ivw_p, method='BH')
fgres$egger_fdr <- p.adjust(fgres$egger_p, method='BH')
fires$ivw_p <- 2*pnorm(-abs(fires$ivw_est/fires$ivw_se))
fires$egger_p <- 2*pnorm(-abs(fires$egger_int/fires$egger_int_se))
fires$ivw_fdr <- p.adjust(fires$ivw_p, method='BH')
fires$egger_fdr <- p.adjust(fires$egger_p, method='BH')
t2dres$het_fdr <- p.adjust(t2dres$het_p, method='BH')
fgres$het_fdr <- p.adjust(fgres$het_p, method='BH')
fires$het_fdr <- p.adjust(fires$het_p, method='BH')

t2dres$model_fdr <- p.adjust(t2dres$fullinds.p, method='BH')
fgres$model_fdr <- p.adjust(fgres$fullinds.p, method='BH')
fires$model_fdr <- p.adjust(fires$fullinds.p, method='BH')

t2dres$gene_name <- gencode[match(t2dres$gene, gencode$gene_id),'gene_name']
fgres$gene_name <- gencode[match(fgres$gene, gencode$gene_id),'gene_name']
fires$gene_name <- gencode[match(fires$gene, gencode$gene_id),'gene_name']

t2dres$gwas <- 'T2D'
fgres$gwas <- 'FG'
fires$gwas <- 'FI'

write.table(t2dres, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/mr_results/t2d_nocor_results_with_info.txt')
write.table(fgres, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/mr_results/fg_nocor_results_with_info.txt')
write.table(fires, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/mr_results/fi_nocor_results_with_info.txt')

#t2dsig <- subset(t2dres, ivw_fdr<=0.05&(egger_fdr>=0.05|is.na(egger_fdr))&model_fdr<=0.05)
#fgsig <- subset(fgres, ivw_fdr<=0.05&(egger_fdr>=0.05|is.na(egger_fdr))&model_fdr<=0.05)
#fisig <- subset(fires, ivw_fdr<=0.05&(egger_fdr>=0.05|is.na(egger_fdr))&model_fdr<=0.05)

fullres <- rbind(t2dres, fgres, fires)
fullres$full_egger_fdr <- p.adjust(fullres$egger_p, method='BH')
fullres$full_ivw_fdr <- p.adjust(fullres$ivw_p, method='BH')
fullres$full_het_fdr <- p.adjust(fullres$het_p, method='BH')
fullres_sig <- subset(fullres, model_fdr<=0.05&(full_egger_fdr>=0.05|is.na(full_egger_fdr))&full_ivw_fdr<=0.05&(full_het_fdr>=0.05|is.na(full_het_fdr)))
write.table(fullres_sig, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/mr_results/overall_significant_nocor_allgwas_together.txt')

publishedloci <- read.table('/well/got2d/apayne/thesis/chapter6/published_eqtl_gwas_coloc.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
publishedloci <- publishedloci[,-6]
publishedloci <- publishedloci[rowSums(publishedloci[,-1])>0,]
publishedloci <- publishedloci[order(publishedloci$ottosson, publishedloci$taneera, publishedloci$fadista, publishedloci$vandebunt, publishedloci$varshney, decreasing=TRUE),]
publishedloci <- publishedloci[,c(1, 3, 4, 2, 5, 6)]
publishedloci[publishedloci==1] <- 'X'
publishedloci[publishedloci==0] <- ''

fullres_sig$published <- fullres_sig$gene_name%in%publishedloci$gene_name

pdf('/well/got2d/apayne/thesis/chapter6/mr_full_results_gwas_vs_model_p.pdf', width=5.9, height=5.9)
ggplot(fullres_sig, aes(x=fullinds.r2, y=-log10(ivw_p))) +
  geom_point(aes(col=published, shape=as.factor(gwas))) +
  geom_text_repel(data=subset(fullres_sig, -log10(ivw_p)>6|fullinds.r2>0.26), aes(label=gene_name), size=2.5) +
  labs(x=expression(paste('Final model ', italic('R')^2, sep='')),
    y=expression(paste('-log'['10']*'(IVW ', italic('p'), ')', sep=''))
  ) +
  scale_colour_manual(values=c('blue', 'red'), labels=c('Previously unpublished', 'Published')) +
  guides(colour=guide_legend(title=''), shape=guide_legend(title='Trait')) +
  theme(legend.position='top')
dev.off()

  labs(x=expression(paste('-log'['10']*'(meta-analysis ', italic('p'), ')', sep='')),
    y=expression(paste('-log'['10']*'(eQTL ', italic('p'), ')', sep=''))
  ) +
  scale_colour_manual(values=c('blue', 'red'), labels=c('Previously unpublished', 'Published')) +




fullres_sig$chr <- as.character(gencode[match(fullres_sig$gene, gencode$gene_id),'seqid'])
fullres_sig$start <- gencode[match(fullres_sig$gene, gencode$gene_id),'start']
fullres_sig$end <- gencode[match(fullres_sig$gene, gencode$gene_id),'end']
fullres_sig$writepos <- paste(substring(fullres_sig$chr, 4, nchar(fullres_sig$chr)), ':', fullres_sig$start, '-', fullres_sig$end, sep='')
out_table <- fullres_sig[,c('gene_name', 'writepos', 'nsnps', 'fullinds.r2', 'gwas', 'ivw_p')]
out_table$fullinds.r2 <- round(out_table$fullinds.r2, 2)
out_table$ivw_p <- signif(out_table$ivw_p, 3)
out_table$gene_name <- paste(out_table$gene_name, '}', sep='')
out_table <- out_table[order(out_table$ivw_p, decreasing=FALSE),]
write.table(out_table, sep='\t&\t', row.names=FALSE, col.names=FALSE, quote=FALSE, eol='\\\\\n\\hline\n\\emph{', file='/well/got2d/apayne/thesis/chapter6/full_mr_results_table.latex')


#######
####
#######
####
#LD stuff
eqtl_results <- read.table('/well/got2d/apayne/thesis/chapter6/eqtl_gwas_coloc/results_after_ld/full_significant_results_ld_with_top_eqtl.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
eqtl_results <- subset(eqtl_results, r2>=0.8)
sum(out_table$gene_name%in%eqtl_results$gene_name)

out_table$snp_gene_trait <- paste(out_table$writepos, out_table$gene_name, out_table$gwas)
match(subset(out_table, gene_name%in%eqtl_results$gene_name)$snp_gene_trait, out_table[order(out_table$ivw_p, decreasing=FALSE),'snp_gene_trait'])

#lncrnas
fullres_sig$gene_type <- gencode[match(fullres_sig$gene, gencode$gene_id),'gene_type']
genesums <- read.table('/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
lncs <- c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA")
genesums$gene_type <- gencode[match(genesums$gene, gencode$gene_id),'gene_type']
sum(genesums$gene_type%in%lncs)
nrow(genesums)
sum(fullres_sig$gene_type%in%lncs)
nrow(fullres_sig)
prop.test(matrix(c(sum(genesums$gene_type%in%lncs), nrow(genesums)-sum(genesums$gene_type%in%lncs), sum(fullres_sig$gene_type%in%lncs), nrow(fullres_sig)-sum(fullres_sig$gene_type%in%lncs)), byrow=TRUE, ncol=2))


########################
# PLOT SOME MR RESULTS #
########################
library(data.table)
library(parallel)
library(genetics)
library(ggplot2)

eqtl <- fread('/well/got2d/apayne/thesis/chapter6/fulleqtls.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
load('/well/got2d/apayne/thesis/chapter6/filtered.Rda')

for (currgwas in c('T2D', 'FI', 'FG')) {
  sigmr <- subset(read.table('/well/got2d/apayne/thesis/chapter6/mr_results/overall_significant_nocor_allgwas_together.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE), gwas==currgwas)
  currsums <- fread(paste('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/', tolower(currgwas), '.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  currcoefs <- read.table('/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

  for (j in sigmr$gene) {
    genecoefs <- subset(currcoefs, gene==j)
    tophit <- sigmr[sigmr$gene==j,]
    snpset <- subset(eqtl, gene==tophit$gene)
    snpset$gwasb <- currsums[match(snpset$SNP, currsums$ID),'EFFECT']/currsums[match(snpset$SNP, currsums$ID),'SE']
    snpset$gwaslogp <- -log10(currsums[match(snpset$SNP, currsums$ID),'P'])
    maxeqtl <- max(abs(snpset$beta))
    mingwas <- min(-(abs(snpset$gwaslogp)))
    snpset$scaledeqtl <- abs(snpset$beta)/maxeqtl
    snpset$scaledgwas <- abs(snpset$gwaslogp)/mingwas
    snpset$ismodel <- FALSE
    snpset[snpset$SNP%in%genecoefs$snp,'ismodel'] <- TRUE
    snpset$ismodel <- factor(snpset$ismodel, levels=c(TRUE, FALSE))

    snpset$samedirection <- sign(snpset$beta*snpset$gwasb)==sign(tophit$ivw_est)
    snpset$snp_pos_mb <- snpset$snp_pos/1000000

    pdf(paste('/well/got2d/apayne/thesis/chapter6/mr_results/plots_tidy/', currgwas, '_', snpset[1,'gene_name'], '_model_snps_eqtl_gwas.pdf', sep=''), width=5.9, height=3.9)
    print(ggplot(snpset) +
    geom_point(aes(x=snp_pos_mb,
        y=scaledeqtl,
        colour=samedirection,
        shape=ismodel,
        size=as.numeric(ismodel)
      )
    ) +
    geom_point(
      aes(x=snp_pos_mb,
        y=scaledgwas,
        colour=samedirection,
        shape=ismodel,
        size=as.numeric(ismodel)
      )
    ) +
    geom_hline(yintercept=0,
      col='black',
      size=1.5
    ) +
    geom_segment(
      aes(x=snpset[1,'gene_start']/1000000,
        y=0,
        xend=snpset[1,'gene_end']/1000000,
        yend=0,
        size=1
      ),
      col='darkorange'
    ) +
    scale_colour_manual(values=c('red', 'blue'),
      guide=guide_legend(override.aes=list(
          shape=c(16, 16),
          colour=c('blue', 'red')
        ),
        order=0,
        title='Directional consistency\n(eQTL and MA)'
      ),
      labels=c('Consistent', 'Inconsistent')
    ) +
    scale_shape_manual(values=c(17, 16),
      guide=guide_legend(override.aes=list(
          shape=c(15, 17),
          colour=c('darkorange', 'black'),
          size=c(1, 2)
        ),
        order=1,
        title='Labels'
      ),
      labels=c(bquote(italic(.(snpset[1,'gene_name']))), 'SNP in model')
    ) +
    scale_size(range=c(2, 1)) +
    guides(size=FALSE) +
    scale_y_continuous(breaks=c(-ceiling(mingwas)/mingwas, 0, 1),
      labels=c(as.character(abs(ceiling(mingwas))), '0', as.character(round(abs(maxeqtl), 2)))
    ) +
    labs(x=paste('Genomic position on chromosome ', as.character(snpset[1,'CHR']), ' (Mb)', sep=''),
      y=bquote('  -log'['10']*'('*.(currgwas)* ~ 'MA '*italic('p')*')    ' ~ '              ' ~ '|'*italic("\u03B2")['eQTL']*'|              ')
    )
    )
    dev.off()
  }

}


library(MendelianRandomization)
fullsigs <- read.table('/well/got2d/apayne/thesis/chapter6/mr_results/overall_significant_nocor_allgwas_together.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

for (i in c('t2d', 'fg', 'fi')) {
  load(paste('/well/got2d/apayne/thesis/chapter6/mr_results/', i, '_mr_results.Rda', sep=''))
  currsigs <- subset(fullsigs, gwas==toupper(i))

  for (j in currsigs$gene) {

    if (!is.na(mr_results[[j]][[3]]$egger_int)) {
      currname <- currsigs[match(j, currsigs$gene),'gene_name']
      currplot <- mr_plot(mr_results[[j]][[1]], line='egger', orientate=TRUE, interactive=FALSE, labels=FALSE, error=FALSE)
      currplot$labels$y <- paste('Genetic association with ', toupper(i), sep='')
      currplot$labels$x <- bquote(paste('Genetic association with' ~ italic(.(currname)) ~ 'expression', sep=''))
      pdf(paste('/well/got2d/apayne/thesis/chapter6/mr_results/egger_plots/', i, '_', currname, '_egger_plot.pdf', sep=''), width=5.9, height=3.9)
        print(plot(currplot))
      dev.off()
    }

  }

}


#########################
# pull metaxcan results #
#########################
library(rtracklayer)
gencode <- readGFF('/well/got2d/rna-seq/resources/gencode.v19.annotation.gtf')
gencode[,1] <- as.character(gencode[,1])
t2d <- read.table('/well/got2d/apayne/thesis/chapter6/metaxcan_results/t2d.metaxcan.res', sep='\t', header=TRUE, stringsAsFactors=FALSE)
fg <- read.table('/well/got2d/apayne/thesis/chapter6/metaxcan_results/fg.metaxcan.res', sep='\t', header=TRUE, stringsAsFactors=FALSE)
fi <- read.table('/well/got2d/apayne/thesis/chapter6/metaxcan_results/fi.metaxcan.res', sep='\t', header=TRUE, stringsAsFactors=FALSE)
t2d$gwas <- 'T2D'
fg$gwas <- 'FG'
fi$gwas <- 'FI'

fullresults <- rbind(t2d, fg, fi)
fullresults$model_fdr <- p.adjust(fullresults$fullinds.r2, method='BH')
fullresults$metaxcan_fdr <- p.adjust(fullresults$metaxcan.p, method='BH')
fullresults$gene_name <- gencode[match(fullresults$gene, gencode$gene_id),'gene_name']
fullresults$chr <- gencode[match(fullresults$gene, gencode$gene_id),'seqid']
fullresults$start <- gencode[match(fullresults$gene, gencode$gene_id),'start']
fullresults$end <- gencode[match(fullresults$gene, gencode$gene_id),'end']
fullresults$chr <- substring(fullresults$chr, 4, nchar(fullresults$chr))
fullresults$chrpos <- paste(fullresults$chr, ':', fullresults$start, '-', fullresults$end, sep='')
fullresults_sig <- subset(fullresults, model_fdr>0.05&metaxcan_fdr<=0.05)
fullresults_sig <- fullresults_sig[order(fullresults_sig$metaxcan.p, decreasing=FALSE),]
write.table(fullresults_sig, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter6/metaxcan_results/full_results_sig.txt')

fullwrite <- fullresults_sig[,c('gene_name', 'chrpos', 'nsnps', 'fullinds.r2', 'gwas', 'metaxcan.p')]
fullwrite[,1] <- paste(fullwrite[,1], '}', sep='')
fullwrite$fullinds.r2 <- round(fullwrite$fullinds.r2, 2)
fullwrite$metaxcan.p <- signif(fullwrite$metaxcan.p, 3)

fullresults_sig$gene_type <- gencode[match(fullresults_sig$gene, gencode$gene_id),'gene_type']
genesums <- read.table('/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
lncs <- c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA")
genesums$gene_type <- gencode[match(genesums$gene, gencode$gene_id),'gene_type']
sum(genesums$gene_type%in%lncs)
nrow(genesums)
sum(fullresults_sig$gene_type%in%lncs)
nrow(fullresults_sig)
prop.test(matrix(c(1575, 11801-1575, 6, 52-6), byrow=TRUE, ncol=2))


write.table(fullwrite, col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t&\t', eol='\\\\\n\\hline\n\\emph{', file='/well/got2d/apayne/thesis/chapter6/metaxcan_results.latex')

#######
fullresults <- read.table('/well/got2d/apayne/thesis/chapter6/metaxcan_results/full_results_sig.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

publishedloci <- read.table('/well/got2d/apayne/thesis/chapter6/published_eqtl_gwas_coloc.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
publishedloci <- publishedloci[,-6]
publishedloci <- publishedloci[rowSums(publishedloci[,-1])>0,]
publishedloci <- publishedloci[order(publishedloci$ottosson, publishedloci$taneera, publishedloci$fadista, publishedloci$vandebunt, publishedloci$varshney, decreasing=TRUE),]
publishedloci <- publishedloci[,c(1, 3, 4, 2, 5, 6)]
publishedloci[publishedloci==1] <- 'X'
publishedloci[publishedloci==0] <- ''

fullresults$published <- fullresults$gene_name%in%publishedloci$gene_name

pdf('/well/got2d/apayne/thesis/chapter6/metaxcan_full_results_gwas_vs_model_p.pdf', width=5.9, height=5.9)
ggplot(fullresults, aes(x=fullinds.r2, y=-log10(metaxcan.p))) +
  geom_point(aes(col=published, shape=as.factor(gwas))) +
  geom_text_repel(data=subset(fullresults, fullinds.r2>=0.15|-log10(metaxcan.p)>9), aes(label=gene_name), size=2.5) +
  labs(x=expression(paste('Final model ', italic('R')^2, sep='')),
    y=expression(paste('-log'['10'], '(S-PrediXcan ', italic('p'), ')', sep=''))
  ) +
  scale_colour_manual(values=c('blue', 'red'), labels=c('Previously unpublished', 'Published')) +
  guides(colour=guide_legend(title=''), shape=guide_legend(title='Trait')) +
  theme(legend.position='top')
dev.off()


ggplot(fullres_sig, aes(x=fullinds.r2, y=-log10(ivw_p))) +
  geom_point(aes(col=published, shape=as.factor(gwas))) +
  geom_text_repel(data=subset(fullres_sig, -log10(ivw_p)>6|fullinds.r2>0.26), aes(label=gene_name), size=2.5) +
  labs(x=expression(paste('Final model ', italic('R')^2, sep='')),
    y=expression(paste('-log'['10']*'(IVW ', italic('p'), ')', sep=''))
  ) +
  scale_colour_manual(values=c('blue', 'red'), labels=c('Previously unpublished', 'Published')) +
  guides(colour=guide_legend(title=''), shape=guide_legend(title='Trait')) +
  theme(legend.position='top')






#######################################
# COMPARE THREE METHODS AND PUBLISHED #
#######################################
library(rtracklayer)
gencode <- readGFF('/well/got2d/rna-seq/resources/gencode.v19.annotation.gtf')
gencode[,1] <- as.character(gencode[,1])

publishedloci <- read.table('/well/got2d/apayne/thesis/chapter6/published_eqtl_gwas_coloc.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
publishedloci <- publishedloci[,-6]
publishedloci <- publishedloci[rowSums(publishedloci[,-1])>0,]
publishedloci <- publishedloci[order(publishedloci$ottosson, publishedloci$taneera, publishedloci$fadista, publishedloci$vandebunt, publishedloci$varshney, decreasing=TRUE),]
publishedloci <- publishedloci[,c(1, 3, 4, 2, 5, 6)]

eqtlres <- read.table('/well/got2d/apayne/thesis/chapter6/eqtl_gwas_coloc/results_after_ld/full_significant_results_ld_with_top_eqtl.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
eqtlres <- subset(eqtlres, r2>=0.8)
mrres <- read.table('/well/got2d/apayne/thesis/chapter6/mr_results/overall_significant_nocor_allgwas_together.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
metares <- read.table('/well/got2d/apayne/thesis/chapter6/metaxcan_results/full_results_sig.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

sum(metares$gene_name%in%publishedloci$gene_name)
sum(metares$gene_name%in%eqtlres$gene_name)
sum(mrres$gene_name%in%eqtlres$gene_name)

unique_per_study <- rbind(data.frame(gene=unique(eqtlres$gene_name), method=rep('eqtl', length(unique(eqtlres$gene_name))), stringsAsFactors=FALSE),
  data.frame(gene=unique(mrres$gene_name), method=rep('mr', length(unique(mrres$gene_name))), stringsAsFactors=FALSE),
  data.frame(gene=unique(metares$gene_name), method=rep('meta', length(unique(metares$gene_name))), stringsAsFactors=FALSE)
)

unique_per_study_table <- table(unique_per_study)

final_table <- cbind(data.frame(gene=rownames(unique_per_study_table), stringsAsFactors=FALSE), as.data.frame.matrix(unique_per_study_table))
final_table_ge2 <- final_table[rowSums(final_table[,-1])>=2,]
final_table_ge2 <- final_table_ge2[order(final_table_ge2$eqtl, final_table_ge2$mr, final_table_ge2$meta, decreasing=TRUE),]
final_table_ge2 <- final_table_ge2[,c('gene', 'eqtl', 'mr', 'meta')]
final_table_ge2[final_table_ge2==1] <- 'X'
final_table_ge2[final_table_ge2==0] <- ''
#write.table(final_table_ge2, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t', file='/well/got2d/apayne/thesis/chapter6/genes_in_more_than_one_method_ld.txt')
#final_table_ge2$gene <- paste(final_table_ge2$gene, '}', sep='')
#write.table(final_table_ge2, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t&\t', eol='\\\\\n\\hline\n\\emph{', file='/well/got2d/apayne/thesis/chapter6/genes_in_more_than_one_method_ld.latex')
final_table_ge2$gene_type <- gencode[match(final_table_ge2$gene, gencode$gene_name),'gene_type']
lncs <- c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA")
sum(final_table_ge2$gene_type%in%lncs)
final_table_ge2[final_table_ge2$gene_type%in%lncs,]


###################################
# GO enrichment of final gene set #
###################################
library("org.Hs.eg.db")
library(GOstats)
library(edgeR)
library(rtracklayer)
library(data.table)
gencode <- readGFF('/well/got2d/GTEx_v7/reference_files/gencode.v19.genes.v7.patched_contigs.gtf')
gencode <- gencode[!gencode$seqid%in%c('MT'),]
load('/well/got2d/apayne/thesis/chapter6/filtered.Rda')
final_genes <- read.table('/well/got2d/apayne/thesis/chapter6/genes_in_more_than_one_method_ld.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

envirogenes <- rownames(v)
envirogenes <- unlist(strsplit(envirogenes, '.', fixed=TRUE))[seq(1, 2*nrow(v), 2)]
environmap <- select(org.Hs.eg.db, keys=envirogenes, columns=c('ENSEMBL', 'ENTREZID'), keytype='ENSEMBL')
universeentrez <- environmap[match(envirogenes, environmap$ENSEMBL),'ENTREZID']
universeentrez <- unique(universeentrez[!is.na(universeentrez)])

siggenes <- gencode[match(final_genes[,1], gencode$gene_name),'gene_id']
#RP11-158I13.2 not available for some reason
siggenes <- siggenes[!is.na(siggenes)]
siggenes <- c(siggenes, 'ENSG00000259080.1')
siggenes <- unlist(strsplit(siggenes, '.', fixed=TRUE))[seq(1, 2*nrow(final_genes), 2)]
sigmap <- select(org.Hs.eg.db, keys=siggenes, columns=c('ENSEMBL', 'ENTREZID'), keytype='ENSEMBL')
sigentrez <- sigmap[match(siggenes, sigmap$ENSEMBL),'ENTREZID']
sigentrez <- unique(sigentrez[!is.na(sigentrez)])

params <- new('GOHyperGParams', geneIds=sigentrez, universeGeneIds=universeentrez, annotation='org.Hs.eg.db', ontology='BP', pvalueCutoff=1, conditional=FALSE, testDirection='over')
hgover <- summary(hyperGTest(params))
hgover$fdr <- p.adjust(hgover$Pvalue)


##########################
# COMPARE WITH CHAPTER 4 #
##########################
library(rtracklayer)
gencode <- readGFF('/well/got2d/rna-seq/resources/gencode.v19.annotation.gtf')
load('/well/got2d/apayne/thesis/chapter6/filtered.Rda')
publishedloci <- read.table('/well/got2d/apayne/thesis/chapter6/published_eqtl_gwas_coloc.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
publishedloci <- publishedloci[,-6]
publishedloci <- publishedloci[rowSums(publishedloci[,-1])>0,]
publishedloci <- publishedloci[order(publishedloci$taneera, publishedloci$fadista, publishedloci$vandebunt, publishedloci$varshney, decreasing=TRUE),]
publishedloci <- publishedloci[,c(1, 3, 4, 2, 5)]

eqtlres <- read.table('/well/got2d/apayne/thesis/chapter6/eqtl_gwas_coloc/full_significant_results.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
mrres <- read.table('/well/got2d/apayne/thesis/chapter6/mr_results/overall_significant_nocor_allgwas_together.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
metares <- read.table('/well/got2d/apayne/thesis/chapter6/metaxcan_results/full_results_sig.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

sum(metares$gene_name%in%publishedloci$gene_name)
sum(metares$gene_name%in%eqtlres$gene_name)
sum(mrres$gene_name%in%eqtlres$gene_name)

unique_per_study <- rbind(data.frame(gene=unique(eqtlres$gene_name), method=rep('eqtl', length(unique(eqtlres$gene_name))), stringsAsFactors=FALSE),
  data.frame(gene=unique(mrres$gene_name), method=rep('mr', length(unique(mrres$gene_name))), stringsAsFactors=FALSE),
  data.frame(gene=unique(metares$gene_name), method=rep('meta', length(unique(metares$gene_name))), stringsAsFactors=FALSE)
)

unique_per_study_table <- table(unique_per_study)

final_table <- cbind(data.frame(gene=rownames(unique_per_study_table), stringsAsFactors=FALSE), as.data.frame.matrix(unique_per_study_table))
final_table_ge2 <- final_table[rowSums(final_table[,-1])>=2,]
final_table_ge2 <- final_table_ge2[order(final_table_ge2$eqtl, final_table_ge2$mr, final_table_ge2$meta, decreasing=TRUE),]
final_table_ge2 <- final_table_ge2[,c('gene', 'eqtl', 'mr', 'meta')]

de_results <- read.table('/well/got2d/apayne/islet_t2d_networks/de_by_t2d.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
islet_sigmod <- read.table('/well/got2d/apayne/thesis/chapter4/islet_work/significant_modules_05_bh/islet.lncRNA_de_modules.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
islet_sigmod$gene_name <- gencode[match(islet_sigmod$gene_id, gencode$gene_id),'gene_name']

de_sig <- subset(de_results, adj.P.Val<=0.05)
subset(final_table_ge2, gene%in%de_sig$gene_name)
subset(final_table_ge2, gene%in%islet_sigmod$gene_name)

mr_de <- subset(final_table_ge2, gene%in%de_sig$gene_name)
mr_de$gene_type <- gencode[match(mr_de$gene, gencode$gene_name),'gene_type']

####HSD17B12 and RP11
g1 <- gencode[match('HSD17B12', gencode$gene_name),'gene_id']
g2 <- gencode[match('RP11-613D13.5', gencode$gene_name),'gene_id']
cor.test(v[g1,], v[g2,])
de_sig[match(mr_de$gene, de_sig$gene_name),]

genecoefs <- read.table('/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
subset(genecoefs, gene==g1)
subset(genecoefs, gene==g2)

load('/well/got2d/apayne/thesis/chapter6/mr_results/t2d_mr_results.Rda')
mr_results[[g1]]$betaX * mr_results[[g1]]$betaY

metax <- read.table('/well/got2d/apayne/thesis/chapter6/metaxcan_results/t2d.metaxcan.res', sep='\t', header=TRUE, stringsAsFactors=FALSE)
subset(metax, gene==g1)

#PROX1 and PROX1-AS1
prox1 <- gencode[match('PROX1', gencode$gene_name),'gene_id']
prox1as <- gencode[match('PROX1-AS1', gencode$gene_name),'gene_id']
cor.test(v[prox1,], v[prox1as,])

#SPC25
library(data.table)
eqtl <- fread('/well/got2d/apayne/thesis/chapter6/fulleqtls.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
head(subset(eqtl, gene==gencode[match('SPC25', gencode$gene_name),'gene_id']))
fg <- read.table('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/fg.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
fg[fg$ID==subset(eqtl, gene==gencode[match('SPC25', gencode$gene_name),'gene_id'])[1,'SNP'],]
subset(eqtl, gene==gencode[match('SPC25', gencode$gene_name),'gene_id']&SNP=='rs560887')
genecoefs <- read.table('/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
subset(genecoefs, gene==gencode[match('SPC25', gencode$gene_name),'gene_id'])

#ACP2
head(subset(eqtl, gene==gencode[match('ACP2', gencode$gene_name),'gene_id']))
fg <- read.table('/well/got2d/apayne/thesis/chapter6/gwas_datasets/cleaned/fg.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
fg[fg$ID%in%subset(genecoefs, gene==gencode[match('ACP2', gencode$gene_name),'gene_id'])[,'snp'],]
subset(genecoefs, gene==gencode[match('ACP2', gencode$gene_name),'gene_id'])

gt_matrix <- matrix(unlist(lapply(match(c(currtest$SNP, topeqtl$SNP), colnames(gt)), FUN=function(x) as.genotype.allele.count(round(gt[,x], 0), alleles=c('a', 'b')))), byrow=FALSE, ncol=2)
colnames(gt_matrix) <- c(currtest$SNP, topeqtl$SNP)
currgts <- data.frame(lapply(colnames(gt_matrix), FUN=function(x) genotype(gt_matrix[,x])))
names(currgts) <- colnames(gt_matrix)

  return(LD(currgts)$r^2)

'rs12222581'



#CTB-12O2.1
library(data.table)
currgene <- gencode[match('CTB-12O2.1', gencode$gene_name),'gene_id']
eqtl <- fread('/well/got2d/apayne/thesis/chapter6/fulleqtls.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
ctbeqtl <- subset(eqtl, gene==currgene)
genesums <- read.table('/well/got2d/apayne/thesis/chapter6/lasso/finalmodels_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
subset(genesums, gene==currgene)



###########################
# HSD17B12 tests
library(sva)
library(edgeR)
library(limma)
library(data.table)
library(rtracklayer)

load('/well/got2d/apayne/thesis/chapter6/filtered.Rda')

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



subt2d <- t2d[colnames(v)%in%rownames(gt)]
subcovs <- svars[colnames(v)%in%rownames(gt),]
subv <- v[,colnames(v)%in%rownames(gt)]
subgt <- gt[match(colnames(subv), rownames(gt)),]

hsd <- gencode[match('HSD17B12', gencode$gene_name),'gene_id']

subresid <- resid(lm(subv[hsd,] ~ subcovs))
subt2d <- as.integer(subt2d=='Yes')



model <- glm(subt2d ~ subresid, family=binomial(link='logit'))
#increased expression -> increased t2d risk

eqtl <- lm(subresid ~ subgt[,'rs1061810'])
#Alt allele -> reduced expression

gwas <- glm(subt2d ~ subgt[,'rs1061810'], family=binomial(link='logit'))
#Alt allele -> reduced T2D risk

snplocs[snplocs$ID=="rs1061810",]
fulldiagram <- fread('/well/got2d/apayne/thesis/chapter6/gwas_datasets/METAANALYSIS_DIAGRAM_SE1.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

fulldiagram[fulldiagram[,1]=='11:43877934',]
subgtrs <- subgt[,'rs1061810']
subgtrsround <- round(subgtrs, 0)

summary(lm(subresid[subt2d==1] ~ subgt[subt2d==1,'rs1061810']))
boxplot(subv[hsd,subt2d==1] ~ as.factor(subgtrsround[subt2d==1]))






