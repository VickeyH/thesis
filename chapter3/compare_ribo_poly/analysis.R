##################################
# Picard for read type summaries #
##################################
library(data.table)
library(parallel)

system2('mkdir', '/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/picard_summaries')
system2('mkdir', '/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/picard_summaries/ribo')
system2('mkdir', '/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/picard_summaries/poly_stranded')
system2('mkdir', '/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/picard_summaries/poly_unstranded')

flagstatsr <- list.files('/well/got2d/rna-seq/data/riboMinus_Oxford_islets/merged/bam')
flagstatsr <- flagstatsr[substring(flagstatsr, nchar(flagstatsr), nchar(flagstatsr))!='i']
flagstatsps <- list.files('/well/got2d/rna-seq/data/Oxford_islets/stranded_libraries/merged/bam')
flagstatsps <- flagstatsps[substring(flagstatsps, nchar(flagstatsps), nchar(flagstatsps))!='i']
flagstatspu <- list.files('/well/got2d/rna-seq/data/Oxford_islets/unstranded_libraries/merged/bam')
flagstatspu <- flagstatspu[substring(flagstatspu, nchar(flagstatspu), nchar(flagstatspu))!='i']

flagstatsr <- flagstatsr[flagstatsr%in%flagstatsps|flagstatsr%in%flagstatspu]

picarder <- function(i) system2('java', c(
  '-jar', '/users/mccarthy/apayne/picard/picard.jar',
  'CollectRnaSeqMetrics',
  paste('I=/well/got2d/rna-seq/data/riboMinus_Oxford_islets/merged/bam/', i, sep=''),
  paste('O=/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/picard_summaries/ribo/', substring(i, 1, nchar(i)-4), '.rna_metrics', sep=''),
  'REF_FLAT=/well/got2d/apayne/ribominus_full/refflat_from_original_gtf.txt',
  'RIBOSOMAL_INTERVALS=/well/got2d/apayne/ribominus_full/rRNA.interval_list',
  'STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND')
)

mclapply(flagstatsr, FUN=picarder, mc.cores=11)


flagstatsps <- flagstatsps[flagstatsps%in%flagstatsr]

picarder <- function(i) system2('java', c(
  '-jar', '/users/mccarthy/apayne/picard/picard.jar',
  'CollectRnaSeqMetrics',
  paste('I=/well/got2d/rna-seq/data/Oxford_islets/stranded_libraries/merged/bam/', i, sep=''),
  paste('O=/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/picard_summaries/poly_stranded/', substring(i, 1, nchar(i)-4), '.rna_metrics', sep=''),
  'REF_FLAT=/well/got2d/apayne/ribominus_full/refflat_from_original_gtf.txt',
  'RIBOSOMAL_INTERVALS=/well/got2d/apayne/ribominus_full/rRNA.interval_list',
  'STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND')
)

mclapply(flagstatsps, FUN=picarder, mc.cores=14)


flagstatspu <- flagstatspu[flagstatspu%in%flagstatsr]

picarder <- function(i) system2('java', c(
  '-jar', '/users/mccarthy/apayne/picard/picard.jar',
  'CollectRnaSeqMetrics',
  paste('I=/well/got2d/rna-seq/data/Oxford_islets/unstranded_libraries/merged/bam/', i, sep=''),
  paste('O=/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/picard_summaries/poly_unstranded/', substring(i, 1, nchar(i)-4), '.rna_metrics', sep=''),
  'REF_FLAT=/well/got2d/apayne/ribominus_full/refflat_from_original_gtf.txt',
  'RIBOSOMAL_INTERVALS=/well/got2d/apayne/ribominus_full/rRNA.interval_list',
  'STRAND_SPECIFICITY=NONE')
)

mclapply(flagstatspu, FUN=picarder, mc.cores=3)


#######################
# PSEUDOMONAS MAPPING #
#######################
system2('/well/got2d/rna-seq/dependencies/bin/STAR', c(
'--runThreadN', '32',
'--genomeDir', '/goddn/2/users/apayne/ribo_deplete_islet/pseudomonas_checking/genomedir',
'--genomeLoad', 'NoSharedMemory',
'--genomeSAindexNbases', '5',
'--readFilesIn', '/goddn/2/users/apayne/ribo_deplete_islet/pe1/WTCHG_218_1.fastq', '/goddn/2/users/apayne/ribo_deplete_islet/pe2/WTCHG_218_2.fastq',
'--outFileNamePrefix', '/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/h531_pseudo/WTCHG_218_pseudomonas.',
'--outMultimapperOrder', 'Random',
'--outSAMtype', 'BAM', 'Unsorted',
'--outSAMattributes', 'NH', 'HI', 'AS', 'NM', 'MD',
'--outSAMunmapped', 'Within',
'--outSAMmapqUnique', '50',
'--outSAMheaderHD', '@HD', 'VN@1.4', 'SO:unsorted',
'--outSAMmultNmax', '1',
'--outFilterType', 'BySJout',
'--outFilterMultimapNmax', '20',
'--outFilterMismatchNmax', '999',
'--outFilterMismatchNoverLmax', '0.04',
'--alignIntronMin', '20',
'--alignIntronMax', '1000000',
'--alignMatesGapMax', '1000000',
'--alignSJoverhangMin', '8',
'--alignSJDBoverhangMin', '1',
'--sjdbScore', '1',
'--twopassMode', 'Basic'
))

system2('samtools', c(
  'flagstat', '/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/h531_pseudo/WTCHG_218_pseudomonas.Aligned.out.bam',
  '>', '/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/h531_pseudo/WTCHG_218_pseudomonas.Aligned.out.bam.flagstat'
))


system2('samtools', c(
  'flagstat', '/well/got2d/rna-seq/data/Oxford_islets/unstranded_libraries/merged/bam/H531.bam',
  '>', '/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/h531_pseudo/polya_sample/H531.bam.flagstat'
))


##################
# COMPILE RESULTS #
###################
sample_x <- 'H531'
flagstatsr <- list.files('/well/got2d/rna-seq/data/riboMinus_Oxford_islets/merged/bam')
flagstatsr <- flagstatsr[substring(flagstatsr, nchar(flagstatsr), nchar(flagstatsr))!='i']
flagstatsr <- flagstatsr[flagstatsr!=paste(sample_x, '.bam', sep='')]
flagstatsps <- list.files('/well/got2d/rna-seq/data/Oxford_islets/stranded_libraries/merged/bam')
flagstatsps <- flagstatsps[substring(flagstatsps, nchar(flagstatsps), nchar(flagstatsps))!='i']
flagstatsr <- flagstatsr[flagstatsr%in%flagstatsps]
flagstatsps <- flagstatsps[flagstatsps%in%flagstatsr]

ribopics <- list.files('/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/picard_summaries/ribo')
polyspics <- list.files('/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/picard_summaries/poly_stranded')

riboinfo <- matrix(0, nrow=length(ribopics), ncol=30)
rownames(riboinfo) <- substring(ribopics, 1, nchar(ribopics)-12)
polysinfo <- matrix(0, nrow=length(polyspics), ncol=30)
rownames(polysinfo) <- substring(polyspics, 1, nchar(polyspics)-12)
colnames(riboinfo) <- colnames(polysinfo) <- scan('/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/picard_summaries/ribo/H264.rna_metrics', sep='\t', skip=6, nlines=1, what='character')

for (i in ribopics) riboinfo[substring(i, 1, nchar(i)-12),] <- unlist(read.table(paste('/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/picard_summaries/ribo/', i, sep=''), sep='\t', header=TRUE, skip=6, nrows=1, stringsAsFactors=FALSE))
for (i in polyspics) polysinfo[substring(i, 1, nchar(i)-12),] <- unlist(read.table(paste('/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/picard_summaries/poly_stranded/', i, sep=''), sep='\t', header=TRUE, skip=6, nrows=1, stringsAsFactors=FALSE))

polyinfo <- polysinfo

riboinfo <- riboinfo[rownames(riboinfo)%in%rownames(polyinfo),]
riboinfo <- riboinfo[rownames(riboinfo)!=sample_x,]
polyinfo <- polyinfo[match(rownames(riboinfo), rownames(polyinfo)),]

ribo_key_stats <- as.data.frame(riboinfo[,c(1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20)])
poly_key_stats <- as.data.frame(polyinfo[,c(1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20)])
colnames(ribo_key_stats) <- colnames(poly_key_stats) <- c('total', 'aligned', 'rRNA', 'coding', 'utr', 'intron', 'intergene', 'rRNAp', 'codingp', 'utrp', 'intronp', 'intergenep')

ribo_key_stats$mapp <- ribo_key_stats[,2]/ribo_key_stats[,1]
poly_key_stats$mapp <- poly_key_stats[,2]/poly_key_stats[,1]
ribo_key_stats$exon <- ribo_key_stats[,4]+ribo_key_stats[,5]
poly_key_stats$exon <- poly_key_stats[,4]+poly_key_stats[,5]
ribo_key_stats$exonp <- ribo_key_stats[,14]/ribo_key_stats[,2]
poly_key_stats$exonp <- poly_key_stats[,14]/poly_key_stats[,2]

plottable <- data.frame(Feature=rep(c('Exonic', 'Intronic', 'Intergenic'), each=2), Method=rep(c('RD', 'PA'), 3), Proportion=c(
  median(ribo_key_stats$exonp),
  median(poly_key_stats$exonp),
  median(ribo_key_stats$intronp),
  median(poly_key_stats$intronp),
  median(ribo_key_stats$intergenep),
  median(poly_key_stats$intergenep)
))

plottable$Feature <- factor(plottable$Feature, levels=c('Exonic', 'Intronic', 'Intergenic'))
plottable$Method <- factor(plottable$Method, levels=c('RD', 'PA'))

median(poly_key_stats[,1]/200)
median(ribo_key_stats[,1]/150)
min(poly_key_stats[,1]/200)
min(ribo_key_stats[,1]/150)
max(poly_key_stats[,1]/200)
max(ribo_key_stats[,1]/150)

sum(poly_key_stats[,2])/(1000000000*200)
sum(poly_key_stats[,1])/(1000000000*200)
sum(ribo_key_stats[,2])/(1000000000*150)
sum(ribo_key_stats[,1])/(1000000000*150)

prop.test(matrix(c(sum(poly_key_stats[,2]), sum(poly_key_stats[,1])-sum(poly_key_stats[,2]), sum(ribo_key_stats[,2]), sum(ribo_key_stats[,1])-sum(ribo_key_stats[,2])), byrow=TRUE, nrow=2))

sum(poly_key_stats[,2])/sum(poly_key_stats[,1])
sum(ribo_key_stats[,2])/sum(ribo_key_stats[,1])

sum(poly_key_stats[,14])/(1000000000*200)
sum(poly_key_stats[,2])/(1000000000*200)
sum(ribo_key_stats[,14])/(1000000000*150)
sum(ribo_key_stats[,2])/(1000000000*150)
sum(poly_key_stats[,14])/sum(poly_key_stats[,2])
sum(ribo_key_stats[,14])/sum(ribo_key_stats[,2])


sum(poly_key_stats[,6])/(1000000000*200)
sum(poly_key_stats[,2])/(1000000000*200)
sum(ribo_key_stats[,6])/(1000000000*150)
sum(ribo_key_stats[,2])/(1000000000*150)
sum(poly_key_stats[,6])/sum(poly_key_stats[,2])
sum(ribo_key_stats[,6])/sum(ribo_key_stats[,2])


sum(poly_key_stats[,7])/(1000000000*200)
sum(ribo_key_stats[,7])/(1000000000*150)
sum(poly_key_stats[,7])/sum(poly_key_stats[,2])
sum(ribo_key_stats[,7])/sum(ribo_key_stats[,2])


library(ggplot2)

pdf('/well/got2d/apayne/thesis/chapter3/rdpa_features.pdf', width=5.9, height=3.9)
ggplot(plottable, aes(Feature, Proportion)) +
  geom_bar(aes(fill=Method), position='dodge', stat='identity') +
  scale_fill_manual(values=c('darkorange', 'darkblue')) +
  theme(legend.position='top')
dev.off()

###################
# effective reads #
###################
median(ribo_key_stats$exon/ribo_key_stats$total) * .809 * (1 - .238) #.809 exon mapping rate, .238 rate of rn7s/rphh1
median(poly_key_stats$exon/poly_key_stats$total) * .860 # exon mapping rate

riboout <- data.frame(reads=round(ribo_key_stats[,1]/75000000, 1), stringsAsFactors=FALSE)
polyout <- data.frame(reads=round(poly_key_stats[,1]/100000000, 1), stringsAsFactors=FALSE)
riboout$mapped <- paste(round(ribo_key_stats[,2]/75000000, 1), ' (', round(ribo_key_stats[,13]*100, 1), '%)', sep='')
polyout$mapped <- paste(round(poly_key_stats[,2]/100000000, 1), ' (', round(poly_key_stats[,13]*100, 1), '%)', sep='')
riboout$rrna <- ribo_key_stats[,3]/75
polyout$rrna <- poly_key_stats[,3]/100
riboout$exon <- round(ribo_key_stats[,15]*100, 1)
polyout$exon <- round(poly_key_stats[,15]*100, 1)
riboout$intron <- round(ribo_key_stats[,11]*100, 1)
polyout$intron <- round(poly_key_stats[,11]*100, 1)
riboout$inter <- round(ribo_key_stats[,12]*100, 1)
polyout$inter <- round(poly_key_stats[,12]*100, 1)

finalout <- NULL

for (i in 1:nrow(riboout)) {
  finalout <- rbind(finalout, riboout[i,])
  finalout <- rbind(finalout, polyout[i,])
}

finalout <- cbind(data.frame(sample=rep(paste('Sample ', 1:nrow(riboinfo), sep=''), each=2), meth=rep(c('RiboZero', 'Poly-A'), nrow(riboinfo)), stringsAsFactors=FALSE), finalout)

write.table(finalout, col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t&\t', eol='\\\\\n\\hline\n', file='/well/got2d/apayne/thesis/chapter3/ribopolysummary.txt')


max(finalout$rrna)
max(c(ribo_key_stats$rRNAp, poly_key_stats$rRNAp))

sultan_table <- matrix(c(101028616, 93699257, 103480788, 96072826, 100766522, 92611619, 109645220, 100836127, 109221096, 103587793, 113936740, 107997806, 118950948, 110394536, 112056778, 103887722), ncol=2, byrow=TRUE)
rownames(sultan_table) <- c('poly1', 'poly2', 'poly3', 'poly4', 'ribo1', 'ribo2', 'ribo3', 'ribo4')
sultan_table <- cbind(sultan_table, sultan_table[,2]/sultan_table[,1])
median(sultan_table[1:4,3])
median(sultan_table[5:8,3])

sum(ribo_key_stats$mapp>poly_key_stats$mapp)

median(riboout$exon)
median(polyout$exon)
sum(riboout$exon>polyout$exon)

median(riboout$intron)
median(polyout$intron)
sum(riboout$intron>polyout$intron)

median(riboout$inter)
median(polyout$inter)
sum(riboout$inter>polyout$inter)


###################
# MAPPING METRICS #
###################
ribomapping <- c()
polysmapping <- c()

for (i in substring(flagstatsr, 1, nchar(flagstatsr)-4)) {
  ribomapping <- c(ribomapping, unlist(read.table(paste('/well/got2d/rna-seq/data/riboMinus_Oxford_islets/merged/counts/', i, '.gene.counts.summary', sep=''), skip=1, nrows=1, header=FALSE, colClasses=c('NULL', 'integer'))))
}
for (i in substring(flagstatsps, 1, nchar(flagstatsps)-4)) {
  polysmapping <- c(polysmapping, unlist(read.table(paste('/well/got2d/rna-seq/data/Oxford_islets/stranded_libraries/merged/counts/', i, '.gene.counts.summary', sep=''), skip=1, nrows=1, header=FALSE, colClasses=c('NULL', 'integer'))))
}

polyinfo <- cbind(polyinfo, matrix(polysmapping[match(substring(flagstatsps, 1, nchar(flagstatsps)-4), rownames(polyinfo))]*100*2, ncol=1))
riboinfo <- cbind(riboinfo, matrix(ribomapping[match(substring(flagstatsr, 1, nchar(flagstatsr)-4), rownames(riboinfo))]*75*2, ncol=1))
colnames(polyinfo)[ncol(polyinfo)] <- 'featurecounts_mapped_bases'
colnames(riboinfo)[ncol(riboinfo)] <- 'featurecounts_mapped_bases'

median((riboinfo[,ncol(riboinfo)])/(riboinfo[,4]+riboinfo[,5]))
median((polyinfo[,ncol(polyinfo)])/(polyinfo[,4]+polyinfo[,5]))
sum(polyinfo[,ncol(polyinfo)])/(1000000000*200)
sum(polyinfo[,4]+polyinfo[,5])/(1000000000*200)
sum(polyinfo[,ncol(polyinfo)])/sum(polyinfo[,4]+polyinfo[,5])
sum(riboinfo[,ncol(riboinfo)])/(1000000000*150)
sum(riboinfo[,4]+riboinfo[,5])/(1000000000*150)
sum(riboinfo[,ncol(riboinfo)])/sum(riboinfo[,4]+riboinfo[,5])


sum((riboinfo[,ncol(riboinfo)])/(riboinfo[,4]+riboinfo[,5])<(polyinfo[,ncol(polyinfo)])/(polyinfo[,4]+polyinfo[,5]))

################################
# Post-quantification analysis #
################################
library(data.table)

finalexpr <- fread('/well/got2d/rna-seq/data/riboMinus_Oxford_islets/09.05.2017.riboMinus_Oxford_islets.gene.counts.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)

rownames(finalexpr) <- finalexpr[,1]
finalexpr <- finalexpr[,-c(1, 2)]
finalexpr <- as.matrix(finalexpr)

oldexpr <- fread('/well/got2d/rna-seq/data/Oxford_islets/stranded_libraries/17.03.2017.stranded_libraries.gene.counts.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
rownames(oldexpr) <- oldexpr[,1]
oldexpr <- oldexpr[,-c(1, 2)]
oldexpr <- as.matrix(oldexpr)

riboexpr <- finalexpr[,colnames(finalexpr)%in%colnames(oldexpr)]
riboexpr <- riboexpr[,colnames(riboexpr)!=sample_x]
polyexpr <- oldexpr[,match(colnames(riboexpr), colnames(oldexpr))]
riboexpr <- riboexpr[match(rownames(polyexpr), rownames(riboexpr)),]

library(rtracklayer)
gencode <- readGFF('/well/got2d/rna-seq/resources/gencode.v19.annotation.gtf')
gencode[,1] <- as.character(gencode[,1])
gencode[,2] <- as.character(gencode[,2])
gencode[,3] <- as.character(gencode[,3])
gencode <- subset(gencode, type=='gene')

lncs <- subset(gencode, gene_type%in%c("processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA", "processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "bidirectional_promoter_lncrna", "macro_lncRNA"))

pcs <- subset(gencode, gene_type=="protein_coding")


ribosums <- rowSums(riboexpr)
names(ribosums) <- rownames(riboexpr)
polysums <- rowSums(polyexpr)
names(polysums) <- rownames(polyexpr)

for (i in names(head(sort(ribosums, decreasing=TRUE), 10))) print(paste(ribosums[i]*100/sum(riboexpr), ' % for ', gencode[match(i, gencode$gene_id),'gene_name'], sep=''))

for (i in names(head(sort(polysums, decreasing=TRUE), 10))) print(paste(polysums[i]*100/sum(polyexpr), ' % for ', gencode[match(i, gencode$gene_id),'gene_name'], sep=''))

for (i in names(head(sort(ribosums, decreasing=TRUE), 10))) print(paste(polysums[i]*100/sum(polyexpr), ' % for ', gencode[match(i, gencode$gene_id),'gene_name'], sep=''))


overnames <- c('RN7SL2', 'RN7SL1', 'RN7SK', 'RPPH1', 'INS', 'GCG')
gids <- gencode[match(overnames, gencode$gene_name),'gene_id']

plotframe <- data.frame(gname=rep(overnames, each=(2*ncol(riboexpr))), sid=rep(colnames(riboexpr), 2*length(overnames)), exprprop=rep(0, 6*2*ncol(riboexpr)), meth=rep(rep(c('RD', 'PA'), each=ncol(riboexpr)), length(overnames)), stringsAsFactors=FALSE)

for (i in overnames) {
  for (j in colnames(riboexpr)) {
      plotframe[plotframe$gname==i&plotframe$sid==j&plotframe$meth=='RD','exprprop'] <- riboexpr[gids[match(i, overnames)],j]/sum(riboexpr[,j])
      plotframe[plotframe$gname==i&plotframe$sid==j&plotframe$meth=='PA','exprprop'] <- polyexpr[gids[match(i, overnames)],j]/sum(polyexpr[,j])
  }
}


library(ggplot2)
plotframe$gname <- factor(plotframe$gname, levels=overnames)
plotframe$meth <- factor(plotframe$meth, levels=c('RD', 'PA'))

pdf('/well/got2d/apayne/thesis/chapter3/top_rd_rna_props.pdf', width=5.9, height=3.9)
ggplot(plotframe, aes(x=gname, y=exprprop, fill=meth)) +
  geom_boxplot() +
  labs(x='', y='Proportion of quantified reads', fill='Method') +
  scale_fill_manual(values=c('darkblue', 'darkorange')) +
  theme(legend.position='top')
dev.off()



#library(matrixStats)
#library(parallel)

#ribocols <- colSums(riboexpr)
#polycols <- colSums(polyexpr)

#ribo_sample_props <- do.call('rbind', mclapply(1:nrow(riboexpr), FUN=function(x) riboexpr[x,]/ribocols, mc.cores=48))
#poly_sample_props <- do.call('rbind', mclapply(1:nrow(polyexpr), FUN=function(x) polyexpr[x,]/polycols, mc.cores=48))

#ribomedians <- rowMedians(ribo_sample_props)
#names(ribomedians) <- rownames(riboexpr)
#polymedians <- rowMedians(poly_sample_props)
#names(polymedians) <- rownames(polyexpr)

#for (i in names(head(sort(ribomedians, decreasing=TRUE), 10))) print(paste(ribomedians[i]*100, ' % for ', gencode[match(i, gencode$gene_id),'gene_name'], sep=''))

#for (i in names(head(sort(polymedians, decreasing=TRUE), 10))) print(paste(polymedians[i]*100, ' % for ', gencode[match(i, gencode$gene_id),'gene_name'], sep=''))


##################################################################
# Now removing top four ribo transcripts for further comparisons #
##################################################################

riboexpr <- riboexpr[-match(names(head(sort(ribosums, decreasing=TRUE)))[1:4], rownames(riboexpr)),]
polyexpr <- polyexpr[-match(names(head(sort(ribosums, decreasing=TRUE)))[1:4], rownames(polyexpr)),]

library(edgeR)

polydge <- DGEList(counts=polyexpr)
polydge <- polydge[rowSums(cpm(polydge)>1, na.rm=TRUE)>=5,]
polydge$counts <- polydge$counts + 1
polydge <- calcNormFactors(polydge)
polyv <- voom(polydge)$E

ribodge <- DGEList(counts=riboexpr)
ribodge <- ribodge[rowSums(cpm(ribodge)>1, na.rm=TRUE)>=5,]
ribodge$counts <- ribodge$counts + 1
ribodge <- calcNormFactors(ribodge)
ribov <- voom(ribodge)$E


ribogenes <- data.frame(gene_id=rownames(ribov), stringsAsFactors=FALSE)
polygenes <- data.frame(gene_id=rownames(polyv), stringsAsFactors=FALSE)
ribogenes <- cbind(ribogenes, data.frame(gene_type=gencode[match(ribogenes$gene_id, gencode$gene_id),'gene_type'], gene_type_simple=rep('other', nrow(ribogenes)), stringsAsFactors=FALSE))
polygenes <- cbind(polygenes, data.frame(gene_type=gencode[match(polygenes$gene_id, gencode$gene_id),'gene_type'], gene_type_simple=rep('other', nrow(polygenes)), stringsAsFactors=FALSE))
ribogenes[ribogenes$gene_id%in%lncs$gene_id,'gene_type_simple'] <- 'lncRNA'
polygenes[polygenes$gene_id%in%lncs$gene_id,'gene_type_simple'] <- 'lncRNA'
ribogenes[ribogenes$gene_id%in%pcs$gene_id,'gene_type_simple'] <- 'protein_coding'
polygenes[polygenes$gene_id%in%pcs$gene_id,'gene_type_simple'] <- 'protein_coding'
ribogenes$totalexpr <- rowSums(ribov[match(ribogenes$gene_id, rownames(ribov)),])
polygenes$totalexpr <- rowSums(polyv[match(polygenes$gene_id, rownames(polyv)),])
ribogenes$method <- rep('RiboMinus', nrow(ribogenes))
polygenes$method <- rep('Poly(a)', nrow(polygenes))


ribooverlaps <- ribogenes[ribogenes$gene_id%in%polygenes$gene_id,]
polyoverlaps <- polygenes[match(ribooverlaps$gene_id, polygenes$gene_id),]
plotoverlaps <- cbind(ribooverlaps, data.frame(polyexpr=polyoverlaps$totalexpr, stringsAsFactors=FALSE))

cor(ribooverlaps$totalexpr, polyoverlaps$totalexpr)
cor(subset(ribooverlaps, gene_type_simple=='protein_coding')$totalexpr, subset(polyoverlaps, gene_type_simple=='protein_coding')$totalexpr)
cor(subset(ribooverlaps, gene_type_simple=='lncRNA')$totalexpr, subset(polyoverlaps, gene_type_simple=='lncRNA')$totalexpr)

colnames(plotoverlaps)[3] <- 'Type'
plotoverlaps[plotoverlaps$Type=='protein_coding','Type'] <- 'Coding'
plotoverlaps$Type <- factor(plotoverlaps$Type, levels=c('Coding', 'lncRNA'))

pdf('/well/got2d/apayne/thesis/chapter3/ribopoly_overlap_counts.pdf', width=5.9, height=5.9)
ggplot(subset(plotoverlaps, Type%in%c('Coding', 'lncRNA')), aes(x=polyexpr, y=totalexpr, col=Type)) +
  geom_point() +
  scale_colour_manual(values=c('darkblue', 'darkorange')) +
  theme(legend.position='top') +
  labs(x='PA total normalised expression',
    y='RD total normalised expression'
  )
dev.off()


ribo_unique <- ribogenes[!ribogenes$gene_id%in%polygenes$gene_id,]
poly_unique <- polygenes[!polygenes$gene_id%in%ribogenes$gene_id,]

table(ribo_unique$gene_type_simple)
table(poly_unique$gene_type_simple)

sultan_reads <- read.table('/well/got2d/apayne/thesis/chapter3/compare_ribo_poly/sultan_lncrna_readcounts.csv', sep=',', header=TRUE, colClasses=c('character', rep('NULL', 5), 'integer', rep('numeric', 6)))

library(edgeR)
sum(sultan_reads$g%in%unlist(strsplit(gencode$gene_id, '.', fixed=TRUE))[seq(1, 2*nrow(gencode), 2)])

pa1_lncs <- subset(sultan_reads, RPKM_PolyA_Qia_RNA>=0.5)$g
pa2_lncs <- subset(sultan_reads, RPKM_PolyA_Tri_RNA>=0.5)$g
rd1_lncs <- subset(sultan_reads, RPKM_RiboZ_Qia_RNA>=0.5)$g
rd2_lncs <- subset(sultan_reads, RPKM_RiboZ_Tri_RNA>=0.5)$g

length(pa1_lncs[pa1_lncs%in%pa2_lncs&pa1_lncs%in%rd1_lncs&pa1_lncs%in%rd2_lncs])

pa_lncs <- unique(c(pa1_lncs, pa2_lncs))
rd_lncs <- unique(c(rd1_lncs, rd2_lncs))
pa_only_lncs <- pa_lncs[!pa_lncs%in%rd_lncs]
rd_only_lncs <- rd_lncs[!rd_lncs%in%pa_lncs]

my_rd_unique <- subset(ribo_unique, gene_type_simple=='lncRNA')$gene_id
my_pa_unique <- subset(poly_unique, gene_type_simple=='lncRNA')$gene_id
my_rd_unique <- unlist(strsplit(my_rd_unique, '.', fixed=TRUE))[seq(1, 2*length(my_rd_unique))]
my_pa_unique <- unlist(strsplit(my_pa_unique, '.', fixed=TRUE))[seq(1, 2*length(my_pa_unique))]

my_pa_full <- unlist(strsplit(polygenes$gene_id, '.', fixed=TRUE))[seq(1, 2*nrow(polygenes))]
my_rd_full <- unlist(strsplit(ribogenes$gene_id, '.', fixed=TRUE))[seq(1, 2*nrow(ribogenes))]

sum(pa_only_lncs%in%my_pa_full&pa_only_lncs%in%my_rd_full)
sum(pa_only_lncs%in%my_pa_unique)
sum(rd_only_lncs%in%my_pa_full&rd_only_lncs%in%my_rd_full)
sum(rd_only_lncs%in%my_rd_unique)

#############################
# Discussion mapping rates
#############################
riboquant <- sum(riboexpr)/(sum(riboinfo[,1])/150)
polyquant <- sum(polyexpr)/(sum(polyinfo[,1])/200)
polyquant/riboquant

ribopc <- riboexpr[rownames(riboexpr)%in%subset(gencode, gene_type=='protein_coding')$gene_id,]
polypc <- polyexpr[rownames(polyexpr)%in%subset(gencode, gene_type=='protein_coding')$gene_id,]
riboqpc <- sum(ribopc)/(sum(riboinfo[,1])/150)
polyqpc <- sum(polypc)/(sum(polyinfo[,1]/200))
polyqpc/riboqpc

















#ribolncrnas <- subset(ribo_unique, gene_type_simple=='lncRNA')$gene_id

#ribolncplot <- rowSums(cpm(DGEList(counts=riboexpr))[ribolncrnas,])
#polylncplot <- rowSums(cpm(DGEList(counts=polyexpr))[ribolncrnas,])

#plotlncs <- data.frame(gene_id=ribolncrnas, riboexpr=ribolncplot, polyexpr=polylncplot)

#pdf('/well/got2d/apayne/thesis/chapter3/ribopoly_ribo_unique_lncrna.pdf', width=5.9, height=5.9)
#ggplot(plotlncs, aes(x=polyexpr, y=riboexpr)) +
#  geom_point() +
#  labs(x='PA total read counts',
#    y='RD total read counts'
#  )# +
#  coord_cartesian(xlim=c(-1, 16), ylim=c(-1, 16))
#dev.off()
