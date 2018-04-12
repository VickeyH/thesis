library(mclapply)

filelist <- list.files('/well/got2d/rna-seq/data/riboMinus_Oxford_islets/merged/bam')
filelist <- filelist[substring(filelist, nchar(filelist)-3, nchar(filelist))!='.bai']
filelist <- filelist[filelist!='H531.bam']
filelist <- substring(filelist, 1, nchar(filelist)-4)

stringtie <- function(files) {#for (file in filelist) {
  system2('mkdir', paste('/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/', files, sep=''))
  system2('/home/apayne/stringtie-1.3.3b/stringtie', c(
    paste('/well/got2d/rna-seq/data/riboMinus_Oxford_islets/merged/bam/', files, '.bam', sep=''),
    '-o', paste('/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/', files, '/out', sep=''),
    '-p', '8',
    '-G', '/well/got2d/apayne/ribominus_full/gencode.v19.annotation.gtf',
    '--rf',
    '-A', paste('/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/', files, '/gene_abund', sep=''),
    '-C', paste('/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/', files, '/full_covered_reference', sep=''))
  )
}

mclapply(filelist, FUN=stringtie, mc.cores=6)


system2('mkdir', '/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/full')
system2('/home/apayne/stringtie-1.3.3b/stringtie', c(
  '--merge',
  '-o', '/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/full/final.gtf',
  '-G', '/well/got2d/apayne/ribominus_full/gencode.v19.annotation.gtf',
  '-m', '200',
  paste('/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/', filelist, '/out', sep=''))
)

system2('/home/apayne/gffcompare/gffcompare', c(
  '-r', '/well/got2d/apayne/ribominus_full/gencode.v19.annotation.gtf',
  '-o', '/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/full/gffcomp',
  '/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/full/final.gtf')
)




###################
# Examine new GFF #
###################
library(rtracklayer)
library(zoo)
library(data.table)

newgff <- readGFF('/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/full/gffcomp.annotated.gtf')
gencode <- readGFF('/well/got2d/apayne/ribominus_full/gencode.v19.annotation.gtf')
newgff[,1] <- as.character(newgff[,1])
newgff[,2] <- as.character(newgff[,2])
newgff[,3] <- as.character(newgff[,3])
gencode[,1] <- as.character(gencode[,1])
gencode[,2] <- as.character(gencode[,2])
gencode[,3] <- as.character(gencode[,3])

newgff <- subset(newgff, seqid%in%gencode$seqid)
newgff[substring(newgff$transcript_id, 1, 1)=='E','ref_gene_id'] <- gencode[match(newgff[substring(newgff$transcript_id, 1, 1)=='E','transcript_id'], gencode$transcript_id),'gene_id']
newgff[substring(newgff$transcript_id, 1, 1)=='E','gene_name'] <- gencode[match(newgff[substring(newgff$transcript_id, 1, 1)=='E','transcript_id'], gencode$transcript_id),'gene_name']
newgff[is.na(newgff$ref_gene_id),'ref_gene_id'] <- newgff[is.na(newgff$ref_gene_id),'gene_id']
newgff$class_code <- na.locf(newgff$class_code)
newgff$gene_id <- newgff$ref_gene_id

gencode_genes <- subset(gencode, type=='gene')

finalgff <- subset(newgff, class_code%in%c('u', '='))

distchecker <- function(i) {
  currsub <- subset(finalgff, gene_id==i)

  if (substring(i, 1, 1)=='E') return(currsub)

  chrom <- currsub$seqid[1]
  cstart <- min(currsub$start)
  cend <- max(currsub$end)

  kb10 <- subset(gencode_genes, seqid==chrom&((cstart<(start-10000)&cend>(start-10000))|(cstart>(start-10000)&cstart<(end+10000))))

  if (nrow(kb10)==0) return(currsub) else return(NULL)
}

filteredgff <- rbindlist(mclapply(unique(finalgff$gene_id), FUN=distchecker, mc.cores=48))
filteredgffdf <- as.data.frame(filteredgff)

outputgff <- c(
  '##description: combined gtf from /well/got2d/rna-seq/resources/gencode.v19.annotation.gtf and stringtie assemblies from 55 islet samples',
  '##provider: GENCODE and stringtie',
  '##contact: apayne@well.ox.ac.uk',
  '##format: gtf',
  '##date: 06/12/2017'
)

gff_strings <- paste(
  filteredgffdf[,1], '\t',
  filteredgffdf[,2], '\t',
  filteredgffdf[,3], '\t',
  filteredgffdf[,4], '\t',
  filteredgffdf[,5], '\t',
  filteredgffdf[,6], '\t',
  filteredgffdf[,7], '\t',
  filteredgffdf[,8], '\t',
  'transcript_id', ' ', filteredgffdf[,9], '; ',
  'gene_id', ' ', filteredgffdf[,10], '; ',
  'xloc', ' ', filteredgffdf[,11], '; ',
  'class_code', ' ', filteredgffdf[,12], '; ',
  'tss_id', ' ', filteredgffdf[,13], '; ',
  'exon_number', ' ', filteredgffdf[,14], '; ',
  'gene_name', ' ', filteredgffdf[,15], '; ',
  'ref_gene_id', ' ', filteredgffdf[,16], '; ',
  'cmp_ref', ' ', filteredgffdf[,17], '; ',
  'contained_in', ' ', filteredgffdf[,18], '; ',
  'cmp_ref_gene', ' ', filteredgffdf[,19], ';',
  sep=''
)

finaloutput <- c(outputgff, gff_strings)
write.table(finaloutput, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter3/novel_genes/final_filtered_known_and_intergenic_only.gtf')

sum(substring(unique(filteredgff$gene_id), 1, 1)=='M')
sum(table(subset(filteredgff, substring(gene_id, 1, 1)=='M'&type=='exon')$gene_id)==1)
sum(table(subset(filteredgff, substring(gene_id, 1, 1)=='M'&type=='exon')$gene_id)>1)

############################
# Re-quantify with new GFF #
############################
library(parallel)
library(data.table)
library(rtracklayer)
filelist <- list.files('/well/got2d/rna-seq/data/riboMinus_Oxford_islets/merged/bam')
filelist <- filelist[substring(filelist, nchar(filelist)-3, nchar(filelist))!='.bai']
filelist <- filelist[filelist!='H531.bam']
filelist <- substring(filelist, 1, nchar(filelist)-4)

starmap <- function(file) {
  system2('/well/got2d/rna-seq/dependencies/bin/featureCounts', c(
    '-p',
    '-t', 'exon',
    '-g', 'ref_gene_id',
    '-s', '2',
    '-T', '2',
    '-B',
    '-C',
    '-a', '/well/got2d/apayne/thesis/chapter3/novel_genes/final_filtered_known_and_intergenic_only.gtf',
    '-o', paste('/well/got2d/apayne/thesis/chapter3/novel_genes/counts/', file, '.gene.counts', sep=''),
    paste('/well/got2d/rna-seq/data/riboMinus_Oxford_islets/merged/bam/', file, '.bam', sep='')
  ))
}

mclapply(filelist, FUN=starmap, mc.cores=27)

gencode <- readGFF('/well/got2d/apayne/thesis/chapter3/novel_genes/final_filtered_known_and_intergenic_only.gtf')
gencode[gencode$gene_name=='NA','gene_name'] <- gencode[gencode$gene_name=='NA','ref_gene_id']

fullcounts <- fread(paste('/well/got2d/apayne/thesis/chapter3/novel_genes/counts/', filelist[1], '.gene.counts', sep=''), sep='\t', skip=1, header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

fullcounts <- cbind(fullcounts[,1], gencode[match(fullcounts[,1], gencode$ref_gene_id),'gene_name'])

for (i in filelist) {
  fullcounts <- cbind(fullcounts, fread(paste('/well/got2d/apayne/thesis/chapter3/novel_genes/counts/', i, '.gene.counts', sep=''), sep='\t', skip=1, header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)[,7])
  print(i)
}

colnames(fullcounts) <- c('gene_id', 'gene_name', filelist)
write.table(fullcounts, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t', file='/well/got2d/apayne/thesis/chapter3/novel_genes/final_counts_using_known_and_intergenic_gtf.gene.counts.tsv')


#######################
# ANALYSE NEW REGIONS #
#######################
library(parallel)
library(data.table)
library(rtracklayer)
library(edgeR)
library(matrixStats)
gencode <- readGFF('/well/got2d/apayne/thesis/chapter3/novel_genes/final_filtered_known_and_intergenic_only.gtf')
gencode[gencode$gene_name=='NA','gene_name'] <- gencode[gencode$gene_name=='NA','ref_gene_id']
gencode_old <- readGFF('/well/got2d/apayne/ribominus_full/gencode.v19.annotation.gtf')
lncrna <- c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense",  "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "macro_lncRNA", "bidirectional_promoter_lncrna")

fullcounts <- fread('/well/got2d/apayne/thesis/chapter3/novel_genes/final_counts_using_known_and_intergenic_gtf.gene.counts.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
rownames(fullcounts) <- fullcounts[,1]
fullcounts <- fullcounts[,-c(1, 2)]
fullcounts <- as.matrix(fullcounts)
dge <- DGEList(counts=fullcounts)
dge <- dge[rowSums(cpm(dge)>1, na.rm=TRUE)>=5,]
dge$counts <- dge$counts + 1
dge <- calcNormFactors(dge)
v <- voom(dge, plot=F)$E

oldcounts <- fread('/well/got2d/rna-seq/data/riboMinus_Oxford_islets/09.05.2017.riboMinus_Oxford_islets.gene.counts.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
rownames(oldcounts) <- oldcounts[,1]
oldcounts <- oldcounts[,-c(1, 2)]
oldcounts <- oldcounts[,match(colnames(v), colnames(oldcounts))]
dgeo <- DGEList(counts=oldcounts)
dgeo <- dgeo[rowSums(cpm(dgeo)>1, na.rm=TRUE)>=5,]
dgeo$counts <- dgeo$counts + 1
dgeo <- calcNormFactors(dgeo)
vo <- voom(dgeo, plot=F)$E

sum(rownames(vo)%in%rownames(v))
mediano <- rowMedians(vo)

sum(substring(rownames(v), 1, 1)=='M')
table(gencode[match(rownames(v)[substring(rownames(v), 1, 1)=='M'], gencode$gene_id),'seqid'])

mediann <- rowMedians(v)
names(mediann) <- rownames(v)
round(range(rank(mediann)[substring(names(mediann), 1, 1)=='M'])/length(mediann), 4)
sum(rank(mediann)[substring(names(mediann), 1, 1)=='M']<=nrow(v)/2)
median(rank(mediann)[substring(names(mediann), 1, 1)=='M'])
median(rank(mediann)[substring(names(mediann), 1, 1)=='M'])/nrow(v)

mediann <- rowMedians(v)
names(mediann) <- rownames(v)
denses <- data.frame(dense=rank(mediann)[substring(rownames(v), 1, 1)=='M'], stringsAsFactors=FALSE)
rownames(denses) <- rownames(v)[substring(rownames(v), 1, 1)=='M']

sum(denses[,1]<=nrow(v)/2)

lncdenses <- data.frame(dense=rank(mediann)[names(mediann)%in%subset(gencode_old, gene_type=='lincRNA')$gene_id], stringsAsFactors=FALSE)
rownames(lncdenses) <- names(mediann)[names(mediann)%in%subset(gencode_old, gene_type=='lincRNA')$gene_id]
range(lncdenses[,1])
sum(lncdenses[,1]<nrow(v)/2)
median(lncdenses[,1])

library(ggplot2)
library(gridExtra)

ks.test(denses[,1], lncdenses[,1])

pdf('/well/got2d/apayne/thesis/chapter3/ning_lnc_ranks.pdf', width=5.9, height=7.9)
p1 <- ggplot(denses, aes(dense)) +
  geom_density() +
  labs(x='NIER expression rank',
    y='Density'
  ) +
  coord_cartesian(xlim=c(0, nrow(v)))

p2 <- ggplot(lncdenses, aes(dense)) +
  geom_density() +
  labs(x='lincRNA expression rank',
    y='Density'
  ) +
  coord_cartesian(xlim=c(0, nrow(v)))

grid.arrange(p1, p2, ncol=1)
dev.off()

lens <- unlist(lapply(rownames(denses), FUN=function(x) max(subset(gencode, gene_id==x)$end)-min(subset(gencode, gene_id==x)$start)))
names(lens) <- rownames(denses)


##################################################
# TEST NIER CORRELATION WITH ALL ANNOTATED GENES #
##################################################
niers <- rownames(v)[substring(rownames(v), 1, 1)=='M']
v_ann <- v[substring(rownames(v), 1, 1)=='E',]
v_niers <- v[substring(rownames(v), 1, 1)=='M',]

corchecker <- function(currnier) {
  rbindlist(lapply(rownames(v_ann), FUN=function(x) {
    currcor <- cor.test(v_niers[currnier,], v_ann[x,])
    return(data.frame(nier=currnier, gene_id=x, correl=currcor$estimate, pval=currcor$p.value, stringsAsFactors=FALSE))
  }))
}

correlations <- rbindlist(mclapply(niers, FUN=corchecker, mc.cores=62))
correlations$p.adj <- p.adjust(correlations$pval, method='BH')
correlations$gene_name <- gencode[match(correlations$gene_id, gencode$gene_id),'gene_name']

correlations_sig <- subset(correlations, correl>=0.85)


























#fullcounts7sl <- fullcounts[!rownames(fullcounts)%in%gencode[gencode$gene_name%in%c('RN7SL2', 'RN7SL1', 'RN7SK', 'RPPH1'),'gene_id'],]
#sum(fullcounts7sl[substring(rownames(fullcounts7sl), 1, 1)=='M',])/sum(as.numeric(fullcounts7sl))
#sum(fullcounts7sl[rownames(fullcounts7sl)%in%subset(gencode_old, gene_type%in%lncrna)$gene_id,])/sum(as.numeric(fullcounts7sl))


#gencode_novex <- subset(gencode, gene_id%in%rownames(v)[substring(rownames(v), 1, 1)=='M'])
#table(gencode_novex[match(unique(gencode_novex$gene_id), gencode_novex$gene_id),'seqid'])

#topning <- rownames(denses)[denses[,1]==max(denses[,1])]
#subset(gencode, gene_id==topning)
#v[topning,]

#cor.test(v['MSTRG.2041',], v['ENSG00000231752.1',])

#genecors <- unlist(lapply(1:nrow(v), FUN=function(x) cor(v['MSTRG.2041',], v[x,])))
#names(genecors) <- rownames(v)

#rank(mediann)[names(mediann)=='ENSG00000231752.1']


#lens <- unlist(lapply(rownames(denses), FUN=function(x) max(subset(gencode, gene_id==x)$end)-min(subset(gencode, gene_id==x)$start)))
#names(lens) <- rownames(denses)




































#################################
# Re-quantify poly with new GFF #
#################################
library(parallel)
library(data.table)
library(rtracklayer)
filelist <- list.files('/well/got2d/rna-seq/data/Oxford_islets/stranded_libraries/merged/bam')
filelist <- filelist[substring(filelist, nchar(filelist)-3, nchar(filelist))!='.bai']
filelist <- substring(filelist, 1, nchar(filelist)-4)
#system2('mkdir', '/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/polya_counts')

#starmap <- function(file) {
#  system2('/home/apayne/subread-1.5.1-source/bin/featureCounts', c(
#    '-p',
#    '-t', 'exon',
#    '-g', 'ref_gene_id',
#    '-s', '2',
#    '-T', '2',
#    '-B',
#    '-C',
#    '-a', '/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/full/final_filtered_known_and_intergenic_only.gtf',
#    '-o', paste('/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/polya_counts/', file, '.gene.counts', sep=''),
#    paste('/well/got2d/rna-seq/data/Oxford_islets/stranded_libraries/merged/bam/', file, '.bam', sep='')
#  ))
#}

#mclapply(filelist, FUN=starmap, mc.cores=28)

gencode <- readGFF('/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/full/final_filtered_known_and_intergenic_only.gtf')
gencode[gencode$gene_name=='NA','gene_name'] <- gencode[gencode$gene_name=='NA','ref_gene_id']

fullcounts <- fread(paste('/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/', filelist[1], '.gene.counts', sep=''), sep='\t', skip=1, header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

fullcounts <- cbind(fullcounts[,1], gencode[match(fullcounts[,1], gencode$ref_gene_id),'gene_name'])

for (i in filelist) {
  fullcounts <- cbind(fullcounts, fread(paste('/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/polya_counts/', i, '.gene.counts', sep=''), sep='\t', skip=1, header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)[,7])
  print(i)
}

colnames(fullcounts) <- c('gene_id', 'gene_name', filelist)
write.table(fullcounts, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t', file='/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/final_polya_counts_using_known_and_intergenic_gtf.gene.counts.tsv')
































#####################################
# Some basic analysis on new counts #
#####################################
library(data.table)
library(edgeR)
library(rtracklayer)
library(sva)

gencode <- readGFF('/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/full/final_filtered_known_and_intergenic_only.gtf')
gencode[gencode$gene_name=='NA','gene_name'] <- gencode[gencode$gene_name=='NA','ref_gene_id']

riboexpr <- fread('/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/final_counts_using_known_and_intergenic_gtf.gene.counts.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
rownames(riboexpr) <- riboexpr[,1]
riboexpr <- riboexpr[,-c(1, 2)]
riboexpr <- as.matrix(riboexpr)

polyexpr <- fread('/well/got2d/apayne/ribominus_full/stringtie_full_gencode_reference/final_polya_counts_using_known_and_intergenic_gtf.gene.counts.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
rownames(polyexpr) <- polyexpr[,1]
polyexpr <- polyexpr[,-c(1, 2)]
polyexpr <- as.matrix(polyexpr)


riboexpr <- riboexpr[,colnames(riboexpr)%in%colnames(polyexpr)]
polyexpr <- polyexpr[,match(colnames(riboexpr), colnames(polyexpr))]

ribodge <- DGEList(counts=riboexpr)
ribodge <- ribodge[rowSums(cpm(ribodge)>1, na.rm=TRUE)>=5,]
ribodge$counts <- ribodge$counts + 1
ribodge <- calcNormFactors(ribodge)
ribov <- voom(ribodge, plot=F)$E
polydge <- DGEList(counts=polyexpr)
polydge <- polydge[rowSums(cpm(polydge)>1, na.rm=TRUE)>=5,]
polydge$counts <- polydge$counts + 1
polydge <- calcNormFactors(polydge)
polyv <- voom(polydge, plot=F)$E

design <- model.matrix(~ as.factor(substring(colnames(ribov), 1, 1)))

ribosvars <- sva(dat=ribov, mod=design)$sv
polysvars <- sva(dat=polyv, mod=design)$sv

ribonew <- matrix(0, nrow=nrow(ribov), ncol=ncol(ribov))
rownames(ribonew) <- rownames(ribov)
colnames(ribonew) <- colnames(ribov)
polynew <- matrix(0, nrow=nrow(polyv), ncol=ncol(polyv))
rownames(polynew) <- rownames(polyv)
colnames(polynew) <- colnames(polyv)

for (k in 1:nrow(ribov)) ribonew[k,] <- resid(lm(ribov[k,] ~ ribosvars + as.factor(substring(colnames(ribov), 1, 1))))
for (k in 1:nrow(polyv)) polynew[k,] <- resid(lm(polyv[k,] ~ polysvars + as.factor(substring(colnames(polyv), 1, 1))))


ribosub <- ribonew[rownames(ribonew)%in%rownames(polynew),]
polysub <- polynew[match(rownames(ribosub), rownames(polynew)),]



plot(density(unlist(lapply(1:nrow(ribosub), FUN=function(x) cor(ribosub[x,], polysub[x,])))))


ribo_nov_unique <- rownames(ribonew)[substring(rownames(ribonew), 1, 1)=='M'&!rownames(ribonew)%in%rownames(polynew)]
nov_both <- rownames(ribonew)[rownames(ribonew)%in%rownames(polynew)&substring(rownames(ribonew), 1, 1)=='M']
poly_nov_unique <- rownames(polynew)[substring(rownames(polynew), 1, 1)=='M'&!rownames(polynew)%in%rownames(ribonew)]


t2d_loci <- read.table('/well/got2d/apayne/T2D_GWAS_loci_30082017.AMah.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)
t2d_loci$CHR <- paste('chr', t2d_loci$CHR, sep='')
t2d_loci <- t2d_loci[,c('SNP', 'CHR', 'POS_37', 'GENE_LABEL', 'Population')]

matches <- NULL

for (i in c(ribo_nov_unique, nov_both)) {
  currgencode <- gencode[match(i, gencode$gene_name),]
  currloci <- subset(t2d_loci, CHR==currgencode$seqid)
  currloci <- subset(currloci, POS_37>=(currgencode$start-10000)&POS_37<(currgencode$end+10000))

  if(nrow(currloci)!=0) {
    currloci$gene_id <- i
    currloci$start <- currgencode$start
    currloci$end <- currgencode$end
    matches <- rbind(matches, currloci)
  }

}

matches$gene_ann_start <- gencode[match(matches$GENE_LABEL, gencode$gene_name),'start']
matches$gene_ann_end <- gencode[match(matches$GENE_LABEL, gencode$gene_name),'end']
matches$gene_ann_id <- gencode[match(matches$GENE_LABEL, gencode$gene_name),'ref_gene_id']
matches$cor <- unlist(lapply(1:nrow(matches), FUN=function(x) 
  if (matches[x,'gene_ann_id']%in%rownames(ribonew)) {
    return(cor(ribonew[matches[x,'gene_id'],], ribonew[matches[x,'gene_ann_id'],]))
  } else {
    return(NA)
  } ))










library(data.table)
library(parallel)

filelist <- list.files('/well/got2d/apayne/GTEx_v7/metabolite_paper/lasso_cluster_files/')
filelist <- filelist[filelist!='templates']

metablist <- read.table('/well/got2d/apayne/GTEx_v7/metabolite_paper/metaboliteMap.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

qval_table <- NULL

for (tiss in filelist) {
  currgwas <- fread(paste('/well/got2d/apayne/GTEx_v7/metabolite_paper/lasso_cluster_files/', tiss, '/final_tables_manual_p/M00054.txt', sep=''), header=TRUE, sep=',', stringsAsFactors=FALSE, data.table=FALSE)
  currgwas$Tissue <- tiss
  qval_table <- rbind(qval_table, currgwas)
}

qval_table$qval_studywide <- p.adjust(qval_table$pred_perf_pval, method='BH')


a <- function(tiss) {
  currtiss <- NULL
  currqvals <- subset(qval_table, Tissue==tiss)

  for (i in metablist[,1]) {
    currgwas <- fread(paste('/well/got2d/apayne/GTEx_v7/metabolite_paper/lasso_cluster_files/', tiss, '/final_tables_manual_p/', i, '.txt', sep=''), header=TRUE, sep=',', stringsAsFactors=FALSE, data.table=FALSE)
    currgwas$metabolite <- i
    currgwas$metabolite_name <- metablist[match(i, metablist[,1]),2]
    currgwas$FDR_by_gwas <- p.adjust(currgwas$pvalue, method='BH')
    currgwas$Tissue <- rep(tiss, nrow(currgwas))
    currgwas$qval_studywide <- currqvals[match(currgwas$gene, currqvals$gene),'qval_studywide']
    currtiss <- rbind(currtiss, currgwas)
  }

  currtiss$FDR_by_tissue <- p.adjust(currtiss$pvalue, method='BH')
  return(currtiss)
}

finallist <- do.call('rbind', mclapply(filelist, FUN=a, mc.cores=43))

finallist$FDR_studywide <- p.adjust(finallist$pvalue, method='BH')

write.table(finallist, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/GTEx_v7/metabolite_paper/metaxcan_full_results.txt')

filterlist <- subset(finallist, FDR_studywide<=0.01&qval_studywide<=0.01)

write.table(filterlist, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/GTEx_v7/metabolite_paper/metaxcan_full_results_filtered_for_model_fit_and_pval_at_0.01.txt')

















ribonov <- rownames(ribodge$counts)[substring(rownames(ribodge$counts), 1, 1)=='M']
polynov <- rownames(polydge$counts)[substring(rownames(polydge$counts), 1, 1)=='M']

novoverlap <- ribonov[ribonov%in%polynov]

cors <- c()
for (i in novoverlap) cors <- c(cors, cor(c(ribodge[i,]$counts), c(polydge[i,]$counts)))
names(cors) <- novoverlap

cors <- cors[order(cors, decreasing=TRUE)]

gencode[gencode$gene_name=='MSTRG.552',]

subset(gencode, start>21500000&end<22500000)




subset(gencode, start>22000000&end<22500000&seqid=='chr1'&type=='transcript')[,c(4, 5, 15)]








dge <- DGEList(counts=expr)
dge <- dge[rowSums(cpm(dge)>1, na.rm=TRUE)>=5,]

novexpressed <- rownames(dge$counts)[substring(rownames(dge$counts), 1, 1)=='M']
novexpressed_counts <- dge$counts[novexpressed,]
novexpressed_counts <- novexpressed_counts[order(rowSums(novexpressed_counts), decreasing=TRUE),]







polynovexpressed <- rownames(polydge$counts)[substring(rownames(polydge$counts), 1, 1)=='M']
polynovexpressed_counts <- polydge$counts[polynovexpressed,]
polynovexpressed_counts <- polynovexpressed_counts[order(rowSums(polynovexpressed_counts), decreasing=TRUE),]



ribo10 <- rownames(novexpressed_counts)[1:10]
poly10 <- rownames(polynovexpressed_counts)[1:10]










gencode[match(rownames(novexpressed_counts)[1:5], gencode$ref_gene_id),]

cors <- lapply(rownames(dge), FUN=function(i) cor(dge['MSTRG.2041',], dge[i,]))
names(cors) <- rownames(dge)
cors <- unlist(cors)
cors <- cors[names(cors)!='MSTRG.2041']

MSTRG.17944

gencode[gencode$gene_name%in%c('MSTRG.2041', 'MSTRG.17944'),]

M2041start <- min(gencode[gencode$gene_name=='MSTRG.2041',4])
M2041end <- max(gencode[gencode$gene_name=='MSTRG.2041',5])
M17944start <- min(gencode[gencode$gene_name=='MSTRG.17944',4])
M17944end <- max(gencode[gencode$gene_name=='MSTRG.17944',5])




########################
# LOOK AT NOTCH2 LOCUS #
########################
load('/well/got2d/apayne/InsPIRE_data/genotypes_maf05.Rda')
phenotype <- read.table('/well/got2d/mvdbunt/InsPIRE/current_analysis/Cov_AllFilesIslets.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)
rownames(genotypes) <- phenotype[match(rownames(genotypes), phenotype$NewSampleID),'GenotypeID']

dge$counts <- dge$counts + 0.5
dge <- calcNormFactors(dge)
dge <- voom(dge)$E

qtlexpr <- dge['MSTRG.15421',colnames(dge)%in%rownames(genotypes)]
qtlgeno <- genotypes[match(names(qtlexpr), rownames(genotypes)),]
qtlgeno <- t(na.omit(t(qtlgeno)))

pvals <- mclapply(colnames(qtlgeno)[1:500000], FUN=function(i) summary(lm(qtlexpr ~ qtlgeno[,i]))$coefficients[2,4], mc.cores=4)

summary(lm(qtlexpr ~ qtlgeno[,'rs665268']))

cors <- lapply(rownames(dge), FUN=function(i) cor(dge['MSTRG.15421',], dge[i,]))
names(cors) <- rownames(dge)
cors <- unlist(cors)
cors <- cors[names(cors)!='MSTRG.15421']




summary(lm(qtlexpr ~ qtlgeno[,'rs2641348']))











dge$counts <- dge$counts + 0.5
dge <- calcNormFactors(dge)
dge <- voom(dge)$E

novexpr <- expr[substring(rownames(expr), 1, 1)=='M',]





















