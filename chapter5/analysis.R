library(sva)
library(edgeR)
library(limma)
library(data.table)
library(rtracklayer)

##################################
# INPUT FOR MATRIX EQTL RAW PREP #
##################################

gencode <- readGFF('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz')
phenotype <- read.table('/well/got2d/mvdbunt/InsPIRE/current_analysis/Cov_AllFilesIslets.csv', header=TRUE, stringsAsFactors=FALSE, sep=",")
expr <- fread('/well/got2d/apayne/InsPIRE_data/islet.overlapsremoved.gene.reads.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
load('/well/got2d/apayne/InsPIRE_data/genotypes_maf05.Rda')
snplocs <- fread('/well/got2d/apayne/InsPIRE_data/snplocs_full.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
snplocs <- snplocs[match(colnames(genotypes), snplocs$ID),]
gc()

rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.matrix(expr)
expr <- expr[,match(rownames(genotypes), colnames(expr))]

sex <- phenotype[match(colnames(expr), phenotype$NewSampleID),'Gender']

genereads <- expr[rowSums(expr>6)>=10,]
gc()
dge <- DGEList(counts=genereads)
dge <- dge[rowSums(cpm(dge)>1, na.rm=TRUE)>=10,]
dge$counts <- dge$counts + 1

design <- model.matrix(~ as.factor(sex))
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=F)$E

svars <- sva(dat=v, mod=design)$sv
covs <- cbind(as.data.frame(svars), data.frame(sex=(sex=='F')+1))
colnames(covs) <- c(paste('SVA', 1:ncol(svars), sep=''), 'SEX')
rownames(covs) <- colnames(v)

gc()
save.image(file='/well/got2d/apayne/thesis/chapter5/basic_filtering_full_workspace.Rda')

#################################
# WRITING INPUT FOR MATRIX EQTL #
#################################

#NOTE THAT I FILTER OUT SNPS WITH ANY SINGLE NA IN ANY INDIVIDUAL

covs <- as.data.frame(t(as.matrix(covs)))
covs <- cbind(data.frame(id=rownames(covs), stringsAsFactors=FALSE), covs)
write.table(covs, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t', file='/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/Covariates.txt')

gencode_genes <- subset(gencode, type=='gene')
gencode_genes$seqid <- as.character(gencode_genes$seqid)

for (i in as.character(1:22)) {
  gc()
  currgencode <- subset(gencode_genes, seqid==paste('chr', i, sep='')&gene_id%in%rownames(v))
  currv <- v[match(currgencode$gene_id, rownames(v)),]
  currsnplocs <- subset(snplocs, CHR==i)
  currgt <- genotypes[,match(currsnplocs$ID, colnames(genotypes))]

  geneloc <- currgencode[,c('gene_id', 'seqid', 'start', 'end')]
  colnames(geneloc) <- c('geneid', 'chr', 's1', 's2')
  geneloc$chr <- substring(geneloc$chr, 4, nchar(geneloc$chr))

  snpsloc <- currsnplocs[,c('ID', 'CHR', 'POS')]
  colnames(snpsloc) <- c('snp', 'chr', 'pos')

  GE <- cbind(data.frame(id=rownames(currv), stringsAsFactors=FALSE), as.data.frame(currv))

  SNP <- cbind(data.frame(id=colnames(currgt), stringsAsFactors=FALSE), as.data.frame(t(currgt)))
  SNP <- SNP[rowSums(is.na(SNP))==0,]

  snpsloc <- snpsloc[match(SNP$id, snpsloc$snp),]

  system2('mkdir', paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, sep=''))
  fwrite(geneloc, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, '/geneloc.txt', sep=''))
  fwrite(snpsloc, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, '/snpsloc.txt', sep=''))
  fwrite(GE, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, '/GE.txt', sep=''))
  fwrite(SNP, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, '/SNP.txt', sep=''))

  print(paste('Finished chromosome ', i, sep=''))
}

for (i in as.character(1:22)) {
  GE <- fread(paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, '/GE.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  covs <- fread('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/Covariates.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  GE_resid <- matrix(unlist(lapply(1:nrow(GE), FUN=function(x) resid(lm(unlist(GE[x,-1]) ~ t(as.matrix(covs[,-1])))))), ncol=ncol(GE)-1, byrow=TRUE)
  colnames(GE_resid) <- colnames(GE)[-1]
  GE_resid <- cbind(data.frame(id=GE[,1], stringsAsFactors=FALSE), as.data.frame(GE_resid))
  fwrite(GE_resid, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, '/GE_resid.txt', sep=''))
}


###############################
# WRITING INPUT FOR MV MODELS #
###############################
# NOTE THAT THIS FOLLOWS DIRECT FROM FIRST CHUNK OF CODE - FILTERING, SV, ETC
gencode[,1] <- as.character(gencode[,1])
gencode[,2] <- as.character(gencode[,2])
gencode[,2] <- as.character(gencode[,3])

final_genes <- subset(gencode, type=='gene')[match(rownames(expr), subset(gencode, type=='gene')$gene_id),]
gt <- genotypes[,colSums(is.na(genotypes))==0]
gc()

snplocs <- snplocs[match(colnames(gt), snplocs$ID),]
final_genes <- subset(final_genes, seqid%in%paste('chr', unique(snplocs$CHR), sep=''))
final_genes <- subset(final_genes, gene_id%in%rownames(v))
v <- v[rownames(v)%in%final_genes$gene_id,]
v <- v[match(final_genes$gene_id, rownames(v)),]

save(covs, v, final_genes, gt, snplocs, file='/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/filtered.Rda')


#######################
# RUNNING MATRIX EQTL #
#######################
library(MatrixEQTL)

for (i in as.character(1:22)) {
  base.dir <- '/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/'
  useModel <- modelLINEAR
  SNP_file_name <- paste(base.dir, 'chr', i, '/SNP.txt', sep='')
  snps_location_file_name <- paste(base.dir, 'chr', i, '/snpsloc.txt', sep='')
#  expression_file_name <- paste(base.dir, 'chr', i, '/GE.txt', sep='')
  expression_file_name <- paste(base.dir, 'chr', i, '/GE_resid.txt', sep='')
  gene_location_file_name <- paste(base.dir, 'chr', i, '/geneloc.txt', sep='')
#  covariates_file_name <- paste(base.dir, 'Covariates.txt', sep='')
#  output_file_name_cis <- paste(base.dir, 'chr', i, '/Matrix_EQTL_cis.txt', sep='')
  output_file_name_cis <- paste(base.dir, 'chr', i, '/Matrix_EQTL_cis_covariate_residuals.txt', sep='')
  output_file_name_tra <- NULL
  pvOutputThreshold_cis <- 1
  pvOutputThreshold_tra <- 0
  errorCovariance <- numeric()
  cisDist <- 1e6
  snps = SlicedData$new()
  snps$fileDelimiter = "\t"
  snps$fileOmitCharacters = "NA"
  snps$fileSkipRows = 1
  snps$fileSkipColumns = 1
  snps$fileSliceSize = 500000
  snps$LoadFile(SNP_file_name)
  gene = SlicedData$new()
  gene$fileDelimiter = "\t"
  gene$fileOmitCharacters = "NA"
  gene$fileSkipRows = 1
  gene$fileSkipColumns = 1
  gene$fileSliceSize = 500000
  gene$LoadFile(expression_file_name)
#  cvrt = SlicedData$new()
#  cvrt$fileDelimiter = "\t"
#  cvrt$fileOmitCharacters = "NA"
#  cvrt$fileSkipRows = 1
#  cvrt$fileSkipColumns = 1
#  if(length(covariates_file_name)>0) {
#    cvrt$LoadFile(covariates_file_name)
#  }
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
#    cvrt = cvrt,
    output_file_name = output_file_name_tra,
    pvOutputThreshold = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)
    rm(list=setdiff(ls(), 'i'))
    gc()
}

topmodeller <- function(i) {
  fitmod <- lm(unlist(GE[GE[,1]==i,-1]) ~ unlist(SNP[SNP[,1]==tops[tops$gene==i,'SNP'],-1]))
  return(list(data.frame(gene=i, snp=tops[tops$gene==i,'SNP'], coefs=summary(fitmod)$coefficients[2,1], stringsAsFactors=FALSE),
    data.frame(gene=i, nsnps=1, fullinds.r2=summary(fitmod)$r.squared, fullinds.p=pf(summary(fitmod)$fstatistic[1], summary(fitmod)$fstatistic[2], summary(fitmod)$fstatistic[3], lower.tail=F))
  ))
}


fullmatrixresults <- NULL
for (i in as.character(1:22)) {
  gc()
  matrixresults <- fread(paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, '/Matrix_EQTL_cis_covariate_residuals.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  matrixresults$CHR <- i
  fullmatrixresults <- rbind(fullmatrixresults, matrixresults)
  print(i)
}

fullmatrixresults$snp_pos <- snplocs[match(fullmatrixresults$SNP, snplocs$ID),'POS']
fullmatrixresults$gene_start <- final_genes[match(fullmatrixresults$gene, final_genes$gene_id),'start']
fullmatrixresults$gene_end <- final_genes[match(fullmatrixresults$gene, final_genes$gene_id),'end']
fullmatrixresults$gene_name <- final_genes[match(fullmatrixresults$gene, final_genes$gene_id),'gene_name']
fullmatrixresults$FDR_overall <- p.adjust(fullmatrixresults$"p-value", method='BH')
fullmatrixresults$strand <- final_genes[match(fullmatrixresults$gene, final_genes$gene_id),'strand']
fullmatrixresults$TSS <- ((fullmatrixresults$strand=='+')*1*fullmatrixresults$gene_start) + ((fullmatrixresults$strand=='-')*1*fullmatrixresults$gene_end)
fullmatrixresults$TES <- ((fullmatrixresults$strand=='-')*1*fullmatrixresults$gene_start) + ((fullmatrixresults$strand=='+')*1*fullmatrixresults$gene_end)

fwrite(fullmatrixresults, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter5/matrix_eqtl/fullmatrixresults_residualised_expr_input.txt')


for (i in as.character(1:22)) {
  GE <- fread(paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, '/GE_resid.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  SNP <- fread(paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, '/SNP.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  geneloc <- fread(paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, '/geneloc.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  snpsloc <- fread(paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, '/snpsloc.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  matrixresults <- fread(paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, '/Matrix_EQTL_cis_covariate_residuals.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  matrixresults$snploc <- snpsloc[match(matrixresults$SNP, snpsloc$snp),'pos']
  matrixresults$start <- geneloc[match(matrixresults$gene, geneloc$geneid),'s1']
  matrixresults$end <- geneloc[match(matrixresults$gene, geneloc$geneid),'s2']
  tops <- matrixresults[match(unique(matrixresults$gene), matrixresults$gene),]

  outframes <- mclapply(tops$gene, FUN=topmodeller, mc.cores=16)
  coefframe <- rbindlist(mclapply(1:length(outframes), FUN=function(x) outframes[[x]][[1]], mc.cores=16))
  genesumframe <- rbindlist(mclapply(1:length(outframes), FUN=function(x) outframes[[x]][[2]], mc.cores=16))
  write.table(coefframe, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, '/top_eqtl_snps_and_coefficients.txt', sep=''))
  write.table(genesumframe, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, '/top_eqtl_gene_model_summaries.txt', sep=''))
  print(paste('Chromosome ', i, ' completed.', sep=''))
}

fullcoefs <- NULL
fullgenesums <- NULL

for (i in as.character(1:22)) {
  currcoefs <- fread(paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, '/top_eqtl_snps_and_coefficients.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  currgenesums <- fread(paste('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/raw_input/chr', i, '/top_eqtl_gene_model_summaries.txt', sep=''), sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
  fullcoefs <- rbind(fullcoefs, currcoefs)
  fullgenesums <- rbind(fullgenesums, currgenesums)
}

write.table(fullcoefs, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter5/matrix_eqtl/fullcoefs.txt')
write.table(fullgenesums, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter5/matrix_eqtl/fullgenesums.txt')


#################
# RUN MV MODELS #
#################
library(glmnet)
library(parallel)

load('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/filtered.Rda')
snplocs <- snplocs[,-c(4, 5)]
gc()

final_genes[,1] <- substring(final_genes[,1], 4, nchar(final_genes[,1]))
final_genes[,1] <- as.integer(final_genes[,1])
gc()

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

    cvfit <- cv.glmnet(testgenotypes, testexpr, nfolds=10, type.measure='mse', alpha=1, keep=TRUE)
#    save(cvfit, file=paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_first_time/full_lassos/', i, '.lasso', sep=''))
    save(cvfit, file=paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_second_time/full_lassos/', i, '.lasso', sep=''))
  }

}

mclapply(final_genes$gene_id, FUN=modeller, mc.cores=22)

#lassolist <- list.files('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_first_time/full_lassos/')
lassolist <- list.files('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_second_time/full_lassos/')
#load(paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_first_time/full_lassos/', lassolist[1], sep=''))
load(paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_second_time/full_lassos/', lassolist[1], sep=''))
lasso_model <- list()
lasso_model[[lassolist[1]]] <- cvfit
counter <- 1

for (i in lassolist[2:length(lassolist)]) {
#  load(paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_first_time/full_lassos/', i, sep=''))
  load(paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_second_time/full_lassos/', i, sep=''))
  lasso_model[[i]] <- cvfit
  counter <- counter + 1
  print(counter)
}


coefs_1 <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_first_time/full_lasso_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
sums_1 <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_first_time/full_lasso_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

coefs_200 <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/full_lasso_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
sums_200 <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/full_lasso_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

overlaps_1 <- sums_1[sums_1$gene%in%sums_200$gene,]
overlaps_200 <- sums_200[match(overlaps_1$gene, sums_200$gene),]
names(lasso_model) <- substring(names(lasso_model), 1, nchar(names(lasso_model))-6)

#save(lasso_model, file='/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_first_time/full_lasso_all_genes.Rda')
save(lasso_model, file='/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_second_time/full_lasso_all_genes.Rda')


#####################
# RUN LASSO 200 ITS #
#####################
load('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/filtered.Rda')
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

    while (counter<200) {
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

save(covs, v, final_genes, gt, snplocs, modeller, file='/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/filtered.Rda')

#basedir <- "/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/" #Base dii of analysis
basedir <- "/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/second_run/" #Base dii of analysis
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

load('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/filtered.Rda')
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
}

for (j in as.character(1:22)) {
#  submitdir <- paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/chr', j, '/submission_scripts/', sep='') #Location of submission files
  submitdir <- paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/second_run/chr', j, '/submission_scripts/', sep='') #Location of submission files
  files <- list.files(submitdir)
  for (file in files) system2('qsub', paste(submitdir, file, sep=""))
}


full_lasso_model <- list()

for (i in as.character(1:22)) {
#  lassolist <- list.files(paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/chr', i, '/final_data_files', sep=''))
  lassolist <- list.files(paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/second_run/chr', i, '/final_data_files', sep=''))

  for (j in lassolist) {
#    load(paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/chr', i, '/final_data_files/', j, sep=''))
    load(paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/second_run/chr', i, '/final_data_files/', j, sep=''))
    full_lasso_model <- c(full_lasso_model, lasso_model)
  }

}

lasso_model <- full_lasso_model
#save(lasso_model, file='/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/full_lasso_all_genes.Rda')
save(lasso_model, file='/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/second_run/full_lasso_all_genes.Rda')


#########################################
# Process full lasso models single runs #
#########################################
library(glmnet)
library(parallel)
library(data.table)

load('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/filtered.Rda')
rm(final_genes, gt, snplocs)
gc()
covs <- as.matrix(covs)

modeller <- function(i) {
  gc()
  coefs <- as.matrix(coef(lasso_model[[i]], s=lasso_model[[i]]$lambda.min))
  coefs <- coefs[coefs[,1]!=0,]
  testexpr <- resid(lm(v[i,] ~ covs))

  if (length(coefs)%in%c(0,1)) {
    return(NA)
  } else {
    snps <- names(coefs)[-1]
    res <- summary(lm(testexpr ~ lasso_model[[i]]$fit.preval[,match(lasso_model[[i]]$lambda.min, lasso_model[[i]]$lambda)]))
    r2 <- res$r.squared
    qval <- pf(res$fstatistic[1], res$fstatistic[2], res$fstatistic[3], lower.tail=F)

    return(list(data.frame(gene=rep(i, length(snps)), snp=snps, coefs=coefs[-1], stringsAsFactors=FALSE),
      data.frame(gene=i, nsnps=length(snps), preval.r2=r2, preval.p=qval)
    ))
  }

}

for (runnum in c('first', 'second')) {
  load(paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_', runnum, '_time/full_lasso_all_genes.Rda', sep=''))
  outframes <- mclapply(names(lasso_model), modeller, mc.cores=48)
  navector <- unlist(lapply(1:length(outframes), FUN=function(x) !is.list(outframes[[x]])))
  goodinds <- (1:length(outframes))[!navector]
  cleanframes <- outframes[goodinds]
  coefframe <- data.frame(rbindlist(mclapply(1:length(cleanframes), FUN=function(x) cleanframes[[x]][[1]], mc.cores=32)))
  genesumframe <- data.frame(rbindlist(mclapply(1:length(cleanframes), FUN=function(x) cleanframes[[x]][[2]], mc.cores=32)))
  write.table(coefframe, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_', runnum, '_time/full_lasso_snps_and_coefficients.txt', sep=''))
  write.table(genesumframe, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_', runnum, '_time/full_lasso_gene_model_summaries_temp.txt', sep=''))
}

r1_coefs <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_first_time/full_lasso_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
r2_coefs <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_second_time/full_lasso_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

r1_sums <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_first_time/full_lasso_gene_model_summaries_temp.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
r2_sums <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_second_time/full_lasso_gene_model_summaries_temp.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

load('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/filtered.Rda')
rm(final_genes, snplocs)
gt <- gt[,colnames(gt)%in%c(r1_coefs$snp, r2_coefs$snp)]
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
  qval <- pf(res$fstatistic[1], res$fstatistic[2], res$fstatistic[3], lower.tail=F)
  return(cbind(whichsums[match(i, whichsums$gene),], data.frame(fullinds.r2=r2, fullinds.p=qval, stringsAsFactors=FALSE)))
}

r1_finalr2 <- rbindlist(mclapply(r1_sums$gene, FUN=r2_calc, r1_coefs, r1_sums, mc.cores=56))
r2_finalr2 <- rbindlist(mclapply(r2_sums$gene, FUN=r2_calc, r2_coefs, r2_sums, mc.cores=56))

write.table(r1_finalr2, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_first_time/full_lasso_gene_model_summaries.txt')
write.table(r2_finalr2, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_second_time/full_lasso_gene_model_summaries.txt')

#system2('rm', c('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_first_time/full_lasso_gene_model_summaries_temp.txt', '/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_second_time/full_lasso_gene_model_summaries_temp.txt'))


#####################################
# PROCESS DATA FROM 200 ITERATIONS  #
#####################################
library(glmnet)
library(parallel)
library(data.table)

load('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/filtered.Rda')
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

for (runnum in c('first', 'second')) {
  gc()
  load(paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/', runnum, '_run/full_lasso_all_genes.Rda', sep=''))
  outframes <- mclapply(names(lasso_model), modeller, mc.cores=48)
  navector <- unlist(lapply(1:length(outframes), FUN=function(x) !is.list(outframes[[x]])))
  goodinds <- (1:length(outframes))[!navector]
  cleanframes <- outframes[goodinds]
  coefframe <- data.frame(rbindlist(mclapply(1:length(cleanframes), FUN=function(x) cleanframes[[x]][[1]], mc.cores=24)))
  genesumframe <- data.frame(rbindlist(mclapply(1:length(cleanframes), FUN=function(x) cleanframes[[x]][[2]], mc.cores=24)))
  write.table(coefframe, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/', runnum, '_run/full_lasso_snps_and_coefficients.txt', sep=''))
  write.table(genesumframe, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/', runnum, '_run/full_lasso_gene_model_summaries_temp.txt', sep=''))
}


library(data.table)
library(parallel)

r1_coefs <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/full_lasso_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
r2_coefs <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/second_run/full_lasso_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

r1_sums <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/full_lasso_gene_model_summaries_temp.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
r2_sums <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/second_run/full_lasso_gene_model_summaries_temp.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

load('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/filtered.Rda')
rm(final_genes, snplocs)
gt <- gt[,colnames(gt)%in%c(r1_coefs$snp, r2_coefs$snp)]
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

r1_finalr2 <- rbindlist(mclapply(r1_sums$gene, FUN=r2_calc, r1_coefs, r1_sums, mc.cores=48))
r2_finalr2 <- rbindlist(mclapply(r2_sums$gene, FUN=r2_calc, r2_coefs, r2_sums, mc.cores=48))

write.table(r1_finalr2, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/full_lasso_gene_model_summaries.txt')
write.table(r2_finalr2, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/second_run/full_lasso_gene_model_summaries.txt')


####################
# FILTER MV MODELS #
####################
coef_table <- read.table('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/full_lasso_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
load('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/filtered.Rda')
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

  while (counter2<25) {
    counter2 <- counter2 + 1
    cvfit <- cv.glmnet(testgenotypes, testexpr, nfolds=10, type.measure='mse', alpha=0, nlambda=100, keep=TRUE)
    MSEs <- cbind(MSEs, cvfit$cvm[match(ls, cvfit$lambda)])
  }

  rownames(MSEs) <- ls
  lambda.min <- as.numeric(names(which.min(rowMeans(MSEs, na.rm=TRUE))))

  write.table(c('a', 'a'), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE, file=paste('/well/got2d/apayne/thesis/chapter5/temp_iteration_counter/', i, sep=''))

  return(list(cvfit, lambda.min, currsnps))
}

its <- unique(coef_table$gene)
names(its) <- its
gc()

finalmodels <- mclapply(its, FUN=modeller, mc.cores=56)
save(finalmodels, file='/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/finalmodels.Rda')

#################################################
# Extract coef and calculate r2 filtered models #
#################################################
library(glmnet)
library(parallel)
library(data.table)

load('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/finalmodels.Rda')
load('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/filtered.Rda')
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

outframes <- mclapply(names(finalmodels), modeller, mc.cores=56)
navector <- unlist(lapply(1:length(outframes), FUN=function(x) !is.list(outframes[[x]])))
goodinds <- (1:length(outframes))[!navector]
cleanframes <- outframes[goodinds]
coefframe <- data.frame(rbindlist(mclapply(1:length(cleanframes), FUN=function(x) cleanframes[[x]][[1]], mc.cores=24)))
genesumframe <- data.frame(rbindlist(mclapply(1:length(cleanframes), FUN=function(x) cleanframes[[x]][[2]], mc.cores=24)))
write.table(coefframe, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/finalmodels_snps_and_coefficients.txt')
write.table(genesumframe, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/finalmodels_gene_model_summaries_temp.txt')


coefs <- read.table('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/finalmodels_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
gsums <- read.table('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/finalmodels_gene_model_summaries_temp.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
load('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/filtered.Rda')
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

finalr2 <- rbindlist(mclapply(gsums$gene, FUN=r2_calc, coefs, gsums, mc.cores=56))

write.table(finalr2, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/finalmodels_gene_model_summaries.txt')


#######################
# LD OLD LASSO VS NEW #
#######################
library(parallel)
library(genetics)
lasso_coefs <- read.table('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/full_lasso_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
final_coefs <- read.table('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/finalmodels_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
load('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/filtered.Rda')

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
ld_list_lasso <- mclapply(its, FUN=ld_checker_lasso, mc.cores=56)
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
ld_list_final <- mclapply(its, FUN=ld_checker_final, mc.cores=63)
names(ld_list_final) <- its

ldsums_lasso <- matrix(0, ncol=20, nrow=length(ld_list_lasso))
rownames(ldsums_lasso) <- names(ld_list_lasso)
colnames(ldsums_lasso) <- as.character(seq(0.05, 1, 0.05))

for (i in 1:length(ld_list_lasso)) {
  ldsums_lasso[names(ld_list_lasso)[i],] <- unlist(lapply(seq(0.05, 1, 0.05), FUN=function(x) sum(ld_list_lasso[[i]]>=x, na.rm=TRUE)))
}

ldwrite_lasso <- cbind(rownames(ldsums_lasso), ldsums_lasso)
colnames(ldwrite_lasso) <- c('gene_id', colnames(ldsums_lasso))
write.table(ldwrite_lasso, quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t', file='/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/ld_r2_summaries_lasso_models.txt')

ldsums_final <- matrix(0, ncol=20, nrow=length(ld_list_final))
rownames(ldsums_final) <- names(ld_list_final)
colnames(ldsums_final) <- as.character(seq(0.05, 1, 0.05))

for (i in 1:length(ld_list_final)) {
  ldsums_final[names(ld_list_final)[i],] <- unlist(lapply(seq(0.05, 1, 0.05), FUN=function(x) sum(ld_list_final[[i]]>=x, na.rm=TRUE)))
}

ldwrite_final <- cbind(rownames(ldsums_final), ldsums_final)
colnames(ldwrite_final) <- c('gene_id', colnames(ldsums_final))
write.table(ldwrite_final, quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t', file='/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/ld_r2_summaries_final_models.txt')



###################################
###################################
### ACTUAL THESIS ANALYSIS WORK ###
###################################
###################################
library(data.table)
library(ggplot2)
library(gridExtra)
matrix_results <- fread('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/fullmatrixresults_residualised_expr_input.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
#colnames(matrix_results) <- c('SNP', 'gene', 'beta', 't-stat', 'p-value', 'FDR', 'CHR', 'snp_pos', 'gene_start', 'gene_end', 'gene_name', 'FDR_overall', 'strand', 'TSS')
#fwrite(matrix_results, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter5/matrix_eqtl/fullmatrixresults_residualised_expr_input.txt')
load('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/filtered.Rda')
gc()

lncs <- subset(final_genes, gene_type%in%c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA"))

totalgenes <- unique(matrix_results$gene)
egenes <- unique(subset(matrix_results, FDR_overall<=0.05)$gene)
length(unique(subset(matrix_results, FDR_overall<=0.01)$gene)) #compare to inspire
snps_per_gene <- table(matrix_results$gene)
gc()

dim(matrix_results)
median(snps_per_gene)
min(snps_per_gene)
max(snps_per_gene)

matrix_coefs <- read.table('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/fullcoefs.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
matrix_sums <- read.table('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/fullgenesums.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

matrix_results_sig <- subset(matrix_results, FDR_overall<=0.05)
matrix_results_sig_lncs <- subset(matrix_results_sig, gene%in%lncs$gene_id)
topsnps <- matrix_results_sig[match(unique(matrix_results_sig$gene), matrix_results_sig$gene),]
topsnps_lncs <- subset(topsnps, gene%in%lncs$gene_id)

boxplots_table <- NULL

for (i in c(0.001, 0.005, 0.01, 0.05, 0.1)) {
  currtab <- subset(matrix_results, FDR_overall<=i)
  currtab <- currtab[match(unique(currtab$gene), currtab$gene),]
  currtab$thresh <- i
  boxplots_table <- rbind(boxplots_table, currtab)
}

#rdist <- runif(length(tssdist), -1000000, 1000000)
#tesdist <- topsnps$TES-topsnps$snp_pos
#topsnps$genelen <- topsnps$gene_end-topsnps$gene_start

pdf('/well/got2d/apayne/thesis/chapter5/eqtl_tss_dist.pdf', width=5.9, height=3.5)
ggplot(subset(boxplots_table, thresh==0.05), aes((snp_pos-TSS)/1000000)) +
  geom_density() +
  labs(x="SNP distance from TSS (Mb)",
    y='Density'
  )
#p2 <- ggplot(topsnps, aes((snp_pos-TES)/1000000)) +
#  geom_density() +
#  labs(x="SNP distance from TES (1000000's bp)",
#    y='Density'
#  )
#grid.arrange(p1, p2, layout_matrix=matrix(c(1, 2), byrow=TRUE, nrow=1))
dev.off()

pdf('/well/got2d/apayne/thesis/chapter5/eqtl_tss_dist_boxplots.pdf', width=5.9, height=3.5)
ggplot(boxplots_table, aes(x=as.factor(thresh), y=abs(snp_pos-TSS)/1000000)) +
  geom_boxplot() +
  labs(x="FDR threshold",
    y="TSS distance (Mb)"
  )
dev.off()


matrix_sums_sig <- subset(matrix_sums, gene%in%topsnps$gene)
matrix_sums_sig_lncs <- subset(matrix_sums_sig, gene%in%lncs$gene_id)
matrix_sums_sig_nolncs <- subset(matrix_sums_sig, !gene%in%lncs$gene_id)

ks.test(matrix_sums_sig$fullinds.r2, matrix_sums_sig_lncs$fullinds.r2)
ks.test(matrix_sums_sig_nolncs$fullinds.r2, matrix_sums_sig_lncs$fullinds.r2)
ks.test(matrix_sums_sig$fullinds.r2, matrix_sums_sig_nolncs$fullinds.r2)

min(matrix_sums_sig$fullinds.r2); max(matrix_sums_sig$fullinds.r2); median(matrix_sums_sig$fullinds.r2)
min(matrix_sums_sig_lncs$fullinds.r2); max(matrix_sums_sig_lncs$fullinds.r2); median(matrix_sums_sig_lncs$fullinds.r2)
min(matrix_sums_sig_nolncs$fullinds.r2); max(matrix_sums_sig_nolncs$fullinds.r2); median(matrix_sums_sig_nolncs$fullinds.r2)



######################################################
# COMPARE FIRST AND SECOND RUN FOR SINGLE INTERATION #
######################################################
library(data.table)
library(ggplot2)

coefs_1 <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_first_time/full_lasso_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
coefs_2 <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_second_time/full_lasso_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

sums_1 <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_first_time/full_lasso_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
sums_2 <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_1_it_second_time/full_lasso_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

overlapgenes <- sums_1$gene[sums_1$gene%in%sums_2$gene]
uniquegenes_1 <- sums_1$gene[!sums_1$gene%in%overlapgenes]
uniquegenes_2 <- sums_2$gene[!sums_2$gene%in%overlapgenes]

compareframe <- data.frame(gene=overlapgenes, stringsAsFactors=FALSE)
compareframe$r1_nsnps <- sums_1[match(overlapgenes, sums_1$gene),'nsnps']
compareframe$r2_nsnps <- sums_2[match(overlapgenes, sums_2$gene),'nsnps']
compareframe$r1_r2 <- sums_1[match(overlapgenes, sums_1$gene),'fullinds.r2']
compareframe$r2_r2 <- sums_2[match(overlapgenes, sums_2$gene),'fullinds.r2']

pdf('/well/got2d/apayne/thesis/chapter5/single_runs_mv_r2.pdf', width=5.9, height=5.9)
ggplot(compareframe, aes(r1_r2, r2_r2)) +
  geom_point() +
  labs(x='First iteration R-squared',
    y='Second iteration R-squared'
  )  
dev.off()

identical_values <- unlist(mclapply(overlapgenes, FUN=function(x) all.equal(subset(coefs_1, gene==x)$snp, subset(coefs_2, gene==x)$snp), mc.cores=48))
sum(identical_values=='TRUE')

cor(compareframe$r1_r2, compareframe$r2_r2)


###################################
# Compare iterations between 200s #
###################################
library(data.table)
library(ggplot2)

coefs_1 <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/full_lasso_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
coefs_2 <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/second_run/full_lasso_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

sums_1 <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/full_lasso_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
sums_2 <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/second_run/full_lasso_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

overlapgenes <- sums_1$gene[sums_1$gene%in%sums_2$gene]
uniquegenes_1 <- sums_1$gene[!sums_1$gene%in%overlapgenes]
uniquegenes_2 <- sums_2$gene[!sums_2$gene%in%overlapgenes]

compareframe <- data.frame(gene=overlapgenes, stringsAsFactors=FALSE)
compareframe$r1_nsnps <- sums_1[match(overlapgenes, sums_1$gene),'nsnps']
compareframe$r2_nsnps <- sums_2[match(overlapgenes, sums_2$gene),'nsnps']
compareframe$r1_r2 <- sums_1[match(overlapgenes, sums_1$gene),'fullinds.r2']
compareframe$r2_r2 <- sums_2[match(overlapgenes, sums_2$gene),'fullinds.r2']

pdf('/well/got2d/apayne/thesis/chapter5/200_runs_mv_r2.pdf', width=5.9, height=5.9)
ggplot(compareframe, aes(r1_r2, r2_r2)) +
  geom_point() +
  labs(x='First iteration R-squared',
    y='Second iteration R-squared'
  )
dev.off()

identical_values <- unlist(mclapply(overlapgenes, FUN=function(x) all.equal(subset(coefs_1, gene==x)$snp, subset(coefs_2, gene==x)$snp), mc.cores=48))
sum(identical_values=='TRUE')

cor(compareframe$r1_r2, compareframe$r2_r2)


###########################
# Compare 200 its to eQTL #
###########################
library(data.table)
library(ggplot2)
library(parallel)

coefs_1 <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/full_lasso_snps_and_coefficients.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
coefs_2 <- fread('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/fullcoefs.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

sums_1 <- fread('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/full_lasso_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
sums_2 <- fread('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/fullgenesums.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

coefs_2 <- coefs_2[match(unique(coefs_1$gene), coefs_2$gene),]
sums_2 <- sums_2[match(unique(coefs_1$gene), sums_2$gene),]

nrow(sums_1)/19503
nrow(subset(sums_1, gene%in%lncs$gene_id))/2503

median(subset(sums_1, gene%in%lncs$gene_id)$nsnps)
median(subset(sums_1, !gene%in%lncs$gene_id)$nsnps)

fulleqtls <- fread('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/fullmatrixresults_residualised_expr_input.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

topeqtlpicked <- unlist(mclapply(sums_2$gene, FUN=function(x) subset(coefs_2, gene==x)$snp%in%subset(coefs_1, gene==x)$snp, mc.cores=48))
sum(topeqtlpicked)

notopeqtls <- subset(fulleqtls, gene%in%unique(sums_2$gene)[!topeqtlpicked])
eqtlranks_percentile <- unlist(mclapply(unique(notopeqtls$gene), FUN=function(x) (min(match(subset(coefs_1, gene==x)$snp, subset(notopeqtls, gene==x)$SNP)))/nrow(subset(notopeqtls, gene==x)), mc.cores=48))
max(eqtlranks_percentile)

notopgenes <- unique(notopeqtls$gene)
rank(sums_1$fullinds.r2)[match(notopgenes, sums_1$gene)]

sum(sums_1$fullinds.r2>=sums_2$fullinds.r2)
max(subset(sums_2, fullinds.r2>sums_1$fullinds.r2)$fullinds.r2/subset(sums_1, fullinds.r2<sums_2$fullinds.r2)$fullinds.r2)

finalsums <- read.table('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/finalmodels_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
finalsums <- finalsums[match(unique(coefs_1$gene), finalsums$gene),]

plotframe <- data.frame(gene=sums_1$gene, eqtlr2=sums_2$fullinds.r2, lassor2=sums_1$fullinds.r2, finalr2=finalsums$fullinds.r2, stringsAsFactors=FALSE)

pdf('/well/got2d/apayne/thesis/chapter5/compare_lasso_and_final_to_top_eqtl.pdf', width=5.9, height=7.9)
p1 <- ggplot(plotframe, aes(eqtlr2, lassor2)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, col='red') +
  labs(x='Top eQTL R-squared',
    y='LASSO model R-squared'
  )

p2 <- ggplot(plotframe, aes(eqtlr2, finalr2)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, col='red') +
  labs(x='Top eQTL R-squared',
    y='Final model R-squared'
  )

grid.arrange(p1, p2, layout_matrix=matrix(c(1, 2), nrow=2))
dev.off()

median(plotframe$lassor2/plotframe$eqtlr2)


######################################
# LD for old vs final models 200 its #
######################################
library(ggplot2)
library(gridExtra)
ldsums_lasso <- read.table('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/ld_r2_summaries_lasso_models.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
ldsums_final <- read.table('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/ld_r2_summaries_final_models.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
rownames(ldsums_lasso) <- ldsums_lasso[,1]
rownames(ldsums_final) <- ldsums_final[,1]
ldsums_lasso <- ldsums_lasso[,-1]
ldsums_final <- ldsums_final[,-1]

lasso_summaries <- read.table('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/full_lasso_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
final_summaries <- read.table('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/finalmodels_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

multisnps <- lasso_summaries$nsnps>1
ldsums_lasso_multi <- ldsums_lasso[multisnps,]

sum(ldsums_lasso_multi[,'X0.8']>=1)

median((lasso_summaries$nsnps-final_summaries$nsnps)/lasso_summaries$nsnps)
min((lasso_summaries$fullinds.r2-final_summaries$fullinds.r2)/lasso_summaries$fullinds.r2)

ldsums_final_multi <- ldsums_final[multisnps,]
sum(ldsums_final_multi[,'X0.8']>=1)

compareframe <- data.frame(gene=lasso_summaries$gene, nsnps.lasso=lasso_summaries$nsnps, nsnps.final=final_summaries$nsnps, r2.lasso=lasso_summaries$fullinds.r2, r2.final=final_summaries$fullinds.r2, stringsAsFactors=FALSE)

pdf('/well/got2d/apayne/thesis/chapter5/nsnps_r2_lasso_final.pdf', width=3.6, height=7)
p1 <- ggplot(compareframe) +
  geom_point(aes(nsnps.lasso, nsnps.final)) +
  xlim(0, max(lasso_summaries$nsnps)) +
  ylim(0, max(lasso_summaries$nsnps)) +
  labs(x='Number of SNPs selected by LASSO',
    y='Number of SNPs remaining after filtering'
  )

p2 <- ggplot(compareframe) +
  geom_point(aes(r2.lasso, r2.final)) +
  xlim(0, 1) +
  ylim(0, 1) +
  labs(x='Model R-squared from full LASSO',
    y='Model R-squared after filtering and ridge'
  )

grid.arrange(p1, p2, layout_matrix=matrix(c(1, 2), nrow=2))
dev.off()


ldsums_lasso <- ldsums_lasso[match(rownames(ldsums_final), rownames(ldsums_lasso)),]
lasso_summaries <- lasso_summaries[match(rownames(ldsums_final), lasso_summaries$gene),]




lasso_summaries <- read.table('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/full_lasso_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
final_summaries <- read.table('/well/got2d/apayne/thesis/chapter5/mv_models/full_models_all_complete_snps/run_200_it/first_run/finalmodels_gene_model_summaries.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
eqtl_summaries <- read.table('/well/got2d/apayne/thesis/chapter5/matrix_eqtl/fullgenesums.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
eqtl_summaries <- eqtl_summaries[match(final_summaries$gene, eqtl_summaries$gene),]

final_summaries_lncs <- subset(final_summaries, gene%in%lncs$gene_id)
final_summaries_nolncs <- subset(final_summaries, !gene%in%lncs$gene_id)

min(final_summaries$fullinds.r2); max(final_summaries$fullinds.r2); median(final_summaries$fullinds.r2)
min(final_summaries_lncs$fullinds.r2); max(final_summaries_lncs$fullinds.r2); median(final_summaries_lncs$fullinds.r2)
min(final_summaries_nolncs$fullinds.r2); max(final_summaries_nolncs$fullinds.r2); median(final_summaries_nolncs$fullinds.r2)

