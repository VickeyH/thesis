library(parallel)
library(data.table)
library(edgeR)
library(rtracklayer)

gencode <- readGFF('/well/got2d/apayne/thesis/chapter3/novel_genes/final_filtered_known_and_intergenic_only.gtf')
gencode[gencode$gene_name=='NA','gene_name'] <- gencode[gencode$gene_name=='NA','ref_gene_id']

fullcounts <- fread('/well/got2d/apayne/thesis/chapter3/novel_genes/final_counts_using_known_and_intergenic_gtf.gene.counts.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
rownames(fullcounts) <- fullcounts[,1]
fullcounts <- fullcounts[,-c(1, 2)]
fullcounts <- as.matrix(fullcounts)
dge <- DGEList(counts=fullcounts)
dge <- dge[rowSums(cpm(dge)>1, na.rm=TRUE)>=5,]
dge$counts <- dge$counts + 1
dge <- calcNormFactors(dge)
v <- voom(dge, plot=F)$E

gencode_novel <- read.table('/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_moran/ap_novel_lincs.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

ak_lnc <- read.table('/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/mmc4.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)

ap_generegion <- do.call('rbind', mclapply(unique(gencode_novel$gene_name), FUN=function(x) {
  data.frame(
    gene_id=x,
    chr=gencode_novel[match(x, gencode_novel$gene_id),'seqid'],
    start=min(subset(gencode_novel, gene_id==x)$start),
    end=max(subset(gencode_novel, gene_id==x)$end),
    stringsAsFactors=FALSE
  )}, mc.cores=32))


ap_overlapping_ak <- NULL
ak_overlapping_ap <- NULL

for (i in ap_generegion$gene_id) {
  istart <- subset(ap_generegion, gene_id==i)$start-10000
  iend <- subset(ap_generegion, gene_id==i)$end+10000
  ichr <- subset(ap_generegion, gene_id==i)$chr
  akchr <- subset(ak_lnc, chr==ichr)
  aksub <- subset(akchr, (start<istart&end>istart)|(start>istart&start<iend))
  if (nrow(aksub)>0) {
    ap_overlapping_ak <- rbind(ap_overlapping_ak, ap_generegion[ap_generegion$gene_id==i,])
    ak_overlapping_ap <- rbind(ak_overlapping_ap, aksub)
  }
  print(i)
}

write.table(ap_overlapping_ak, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/ap_novel_genes_overlapping_ak_novel_genes.txt')
write.table(ak_overlapping_ap, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/ak_novel_genes_overlapping_ap_novel_genes.txt')

ak_gff <- readGFF('/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/final_filtered_known_and_intergenic_only_akerman.gtf')
ak_filter_genes <- unique(ak_gff[substring(ak_gff$gene_id, 1, 1)=='H','gene_id'])
sum(ak_overlapping_ap$LNC_name%in%ak_filter_genes)

head(ak_overlapping_ap[!ak_overlapping_ap$LNC_name%in%ak_filter_genes,], 2)


################################
# RUN MO THROUGH SAME PIPELINE #
################################
library(parallel)
ak <- read.table('/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/mmc4.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)
ak <- ak[is.na(ak$Ensembl.gene.name),]

ak_final_transcripts <- data.frame(
  seqid=ak$chr,
  source=ak$LncRNA.transcript.type,
  type=rep('transcript', nrow(ak)),
  start=ak$start,
  end=ak$end,
  score=rep(NA, nrow(ak)),
  strand=ak$strand,
  phase=rep(NA, nrow(ak)),
  transcript_id=ak$LNC_name,
  gene_id=ak$LNC_name,
  xloc=rep(NA, nrow(ak)),
  class_code=rep('=', nrow(ak)),
  tss_id=rep(NA, nrow(ak)),
  exon_number=rep(NA, nrow(ak)),
  gene_name=ak$LNC_name,
  ref_gene_id=ak$LNC_name,
  cmp_ref=rep(NA, nrow(ak)),
  contained_in=rep(NA, nrow(ak)),
  cmp_ref_gene=rep(NA, nrow(ak)),
  stringsAsFactors=FALSE
)

exoner <- function(i) {
  exons <- ak[i,'blockCount']

  if (exons==1) {
    return(data.frame(
      seqid=ak[i,'chr'],
      source=rep(ak[i,'LncRNA.transcript.type'], exons),
      type=rep('exon', exons),
      start=ak[i,'start'],
      end=ak[i,'end'],
      score=rep(NA, exons),
      strand=ak[i,'strand'],
      phase=rep(NA, exons),
      transcript_id=ak[i,'LNC_name'],
      gene_id=ak[i,'LNC_name'],
      xloc=rep(NA, exons),
      class_code=rep('=', exons),
      tss_id=rep(NA, exons),
      exon_number=1,
      gene_name=ak[i,'LNC_name'],
      ref_gene_id=ak[i,'LNC_name'],
      cmp_ref=rep(NA, exons),
      contained_in=rep(NA, exons),
      cmp_ref_gene=rep(NA, exons),
      stringsAsFactors=FALSE
    ))
  } else {
    sizes <- as.integer(unlist(strsplit(ak[i,'blockSize'], ',', fixed=TRUE)))
    starts <- as.integer(unlist(strsplit(ak[i,'blockStart'], ',', fixed=TRUE)))+ak[i,'start']
    ends <- starts + sizes
    return(data.frame(
      seqid=rep(ak[i,'chr'], exons),
      source=rep(ak[i,'LncRNA.transcript.type'], exons),
      type=rep('exon', exons),
      start=starts,
      end=ends,
      score=rep(NA, exons),
      strand=rep(ak[i,'strand'], exons),
      phase=rep(NA, exons),
      transcript_id=rep(ak[i,'LNC_name'], exons),
      gene_id=rep(ak[i,'LNC_name'], exons),
      xloc=rep(NA, exons),
      class_code=rep('=', exons),
      tss_id=rep(NA, exons),
      exon_number=1:exons,
      gene_name=rep(ak[i,'LNC_name'], exons),
      ref_gene_id=rep(ak[i,'LNC_name'], exons),
      cmp_ref=rep(NA, exons),
      contained_in=rep(NA, exons),
      cmp_ref_gene=rep(NA, exons),
      stringsAsFactors=FALSE
    ))
  }
}

ak_final_exons <- do.call('rbind', lapply(1:nrow(ak), FUN=exoner))

write.table(ak_final_transcripts, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/akerman_transcripts.txt')
write.table(ak_final_exons, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/akerman_exons.txt')

######################################
# REPLACE STRINGTIE EXONS WITH MORAN #
######################################
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
finalgff <- finalgff[substring(finalgff$gene_id, 1, 1)=='E',]
aktrans <- read.table('/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/akerman_transcripts.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
akexons <- read.table('/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/akerman_exons.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
finalgff <- rbind(finalgff, aktrans, akexons)

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
write.table(finaloutput, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE, file='/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/final_filtered_known_and_intergenic_only_akerman.gtf')

sum(substring(unique(filteredgff$gene_id), 1, 1)=='H')

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
    '-T', '4',
    '-B',
    '-C',
    '-a', '/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/final_filtered_known_and_intergenic_only_akerman.gtf',
    '-o', paste('/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/counts/', file, '.gene.counts', sep=''),
    paste('/well/got2d/rna-seq/data/riboMinus_Oxford_islets/merged/bam/', file, '.bam', sep='')
  ))
}

mclapply(filelist, FUN=starmap, mc.cores=27)

gencode <- readGFF('/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/final_filtered_known_and_intergenic_only_akerman.gtf')
gencode[gencode$gene_name=='NA','gene_name'] <- gencode[gencode$gene_name=='NA','ref_gene_id']

fullcounts <- fread(paste('/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/counts/', filelist[1], '.gene.counts', sep=''), sep='\t', skip=1, header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

fullcounts <- cbind(fullcounts[,1], gencode[match(fullcounts[,1], gencode$ref_gene_id),'gene_name'])

for (i in filelist) {
  fullcounts <- cbind(fullcounts, fread(paste('/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/counts/', i, '.gene.counts', sep=''), sep='\t', skip=1, header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)[,7])
  print(i)
}

colnames(fullcounts) <- c('gene_id', 'gene_name', filelist)
write.table(fullcounts, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t', file='/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/final_counts_using_known_and_intergenic_gtf_akerman.gene.counts.tsv')


#############
# COMPARING #
#############
library(parallel)
library(data.table)
library(rtracklayer)
library(edgeR)
library(matrixStats)
gencodeak <- readGFF('/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/final_filtered_known_and_intergenic_only_akerman.gtf')
gencodeak[gencodeak$gene_name=='NA','gene_name'] <- gencodeak[gencodeak$gene_name=='NA','ref_gene_id']
gencodeap <- readGFF('/well/got2d/apayne/thesis/chapter3/novel_genes/final_filtered_known_and_intergenic_only.gtf')
gencodeap[gencodeap$gene_name=='NA','gene_name'] <- gencodeap[gencodeap$gene_name=='NA','ref_gene_id']
gencode_old <- readGFF('/well/got2d/apayne/ribominus_full/gencode.v19.annotation.gtf')
lncrna <- c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense",  "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "macro_lncRNA", "bidirectional_promoter_lncrna")

gencodeaknov <- gencodeak[substring(gencodeak$gene_id, 1, 1)=='H',]
length(unique(gencodeaknov$gene_id))
table(gencodeaknov[match(unique(gencodeaknov$gene_id), gencodeaknov$gene_id),'source'])

akerman <- fread('/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/final_counts_using_known_and_intergenic_gtf_akerman.gene.counts.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
rownames(akerman) <- akerman[,1]
akerman <- akerman[,-c(1, 2)]
akerman <- as.matrix(akerman)
dge <- DGEList(counts=akerman)
dge <- dge[rowSums(cpm(dge)>1, na.rm=TRUE)>=5,]
dge$counts <- dge$counts + 1
dge <- calcNormFactors(dge)
vak <- voom(dge, plot=F)$E

ap <- fread('/well/got2d/apayne/thesis/chapter3/novel_genes/final_counts_using_known_and_intergenic_gtf.gene.counts.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
rownames(ap) <- ap[,1]
ap <- ap[,-c(1, 2)]
ap <- as.matrix(ap)
dge <- DGEList(counts=ap)
dge <- dge[rowSums(cpm(dge)>1, na.rm=TRUE)>=5,]
dge$counts <- dge$counts + 1
dge <- calcNormFactors(dge)
vap <- voom(dge, plot=F)$E

oldcounts <- fread('/well/got2d/rna-seq/data/riboMinus_Oxford_islets/09.05.2017.riboMinus_Oxford_islets.gene.counts.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
rownames(oldcounts) <- oldcounts[,1]
oldcounts <- oldcounts[,-c(1, 2)]
oldcounts <- oldcounts[,match(colnames(v), colnames(oldcounts))]
dgeo <- DGEList(counts=oldcounts)
dgeo <- dgeo[rowSums(cpm(dgeo)>1, na.rm=TRUE)>=5,]
dgeo$counts <- dgeo$counts + 1
dgeo <- calcNormFactors(dgeo)
vold <- voom(dgeo, plot=F)$E

sum(rownames(vold)%in%rownames(vap))
sum(rownames(vold)%in%rownames(vak))
medianold <- rowMedians(vo)
medianak <- rowMedians(vak)
medianap <- rowMedians(vap)

sum(substring(rownames(vap), 1, 1)=='M')
sum(substring(rownames(vak), 1, 1)=='H')

table(gencodeak[match(rownames(vak)[substring(rownames(vak), 1, 1)=='H'], gencodeak$gene_id),'source'])

ap_genes <- rownames(vap)[substring(rownames(vap), 1, 1)=='M']
ak_genes <- rownames(vak)[substring(rownames(vak), 1, 1)=='H']

table(gencodeap[match(rownames(vap)[substring(rownames(vap), 1, 1)=='M'], gencodeap$gene_id),'seqid'])
table(gencodeak[match(rownames(vak)[substring(rownames(vak), 1, 1)=='H'], gencodeak$gene_id),'seqid'])

#ap_overlapping_ak <- read.table('/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/ap_novel_genes_overlapping_ak_novel_genes.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)
#ak_overlapping_ap <- read.table('/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/ak_novel_genes_overlapping_ap_novel_genes.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)

#sum(ak_genes%in%ak_overlapping_ap$LNC_name)

ap_generegion <- data.frame(gene_id=ap_genes, stringsAsFactors=FALSE)
ap_generegion$chr <- gencodeap[match(ap_generegion$gene_id, gencodeap$gene_id),'seqid']
ap_generegion$start <- lapply(ap_generegion$gene_id, FUN=function(x) min(subset(gencodeap, gene_id==x)$start))
ap_generegion$end <- lapply(ap_generegion$gene_id, FUN=function(x) max(subset(gencodeap, gene_id==x)$end))
ak_lnc <- read.table('/well/got2d/apayne/thesis/chapter3/novel_genes/compare_to_akerman/mmc4.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)
ak_lnc <- subset(ak_lnc, LNC_name%in%rownames(vak))

ap_overlapping_ak <- NULL
ak_overlapping_ap <- NULL

for (i in ap_generegion$gene_id) {
  istart <- as.integer(subset(ap_generegion, gene_id==i)$start)
  iend <- as.integer(subset(ap_generegion, gene_id==i)$end)
  ichr <- subset(ap_generegion, gene_id==i)$chr
  akchr <- subset(ak_lnc, chr==ichr)
  aksub <- subset(akchr, (start<=istart&end>=istart)|(start>=istart&start<=iend))
  if (nrow(aksub)>0) {
    ap_overlapping_ak <- rbind(ap_overlapping_ak, ap_generegion[ap_generegion$gene_id==i,])
    ak_overlapping_ap <- rbind(ak_overlapping_ap, aksub)
  }
  print(i)
}

ak_overlapping_ap <- ak_overlapping_ap[match(unique(ak_overlapping_ap$LNC_name), ak_overlapping_ap$LNC_name),]
ap_overlapping_ak <- ap_overlapping_ak[match(unique(ap_overlapping_ak$gene_id), ap_overlapping_ak$gene_id),]


sum(ak_genes%in%ak_overlapping_ap$LNC_name)
table(gencodeak[match(ak_genes[ak_genes%in%ak_overlapping_ap$LNC_name], gencodeak$gene_id),'source'])


ak_overlapping_lnc_full <- ak_lnc[match(ak_genes[ak_genes%in%ak_overlapping_ap$LNC_name], ak_lnc$LNC_name),]

> ak_overlapping_lnc_full[,1:12]
> ap_overlapping_ak














ap_10k_ak <- NULL
ak_10k_ap <- NULL

for (i in ap_generegion$gene_id) {
  istart <- as.integer(subset(ap_generegion, gene_id==i)$start)-10000
  iend <- as.integer(subset(ap_generegion, gene_id==i)$end)+10000
  ichr <- subset(ap_generegion, gene_id==i)$chr
  akchr <- subset(ak_lnc, chr==ichr)
  aksub <- subset(akchr, (start<=istart&end>=istart)|(start>=istart&start<=iend))
  if (nrow(aksub)>0) {
    ap_10k_ak <- rbind(ap_10k_ak, ap_generegion[ap_generegion$gene_id==i,])
    ak_10k_ap <- rbind(ak_10k_ap, aksub)
  }
  print(i)
}

sum(ak_genes%in%ak_10k_ap$LNC_name)
table(gencodeak[match(ak_genes[ak_genes%in%ak_10k_ap$LNC_name], gencodeak$gene_id),'source'])

akgeneranks <- rank(rank(rowMedians(vak[ak_genes,])))
names(akgeneranks) <- ak_genes
akgeneranks[ak_genes[ak_genes%in%ak_10k_ap$LNC_name]]





mediann <- rowMedians(v)
names(mediann) <- rownames(v)
range(rank(mediann)[substring(names(mediann), 1, 1)=='M'])
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


