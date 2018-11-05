# process miseq targeted sequencing of crispr edits
# modifierd 11/30/2016 to handle Olga's and Connie's multi-multiplexed libraries (multiple sites per index)

#find /media/jbloom/d1/CRISPR_base_editor/targeted_sequencing/CRISPR_Pilot2-33767777/  -type f -iname '*.fastq.gz' -exec mv {} /media/jbloom/d1/CRISPR_base_editor/targeted_sequencing/fastq/ \;
# then used thunar in linux for batch rename to parse off text between name given in csv key file and name from miseq

library(Biostrings)
library(seqinr)
library(seqLogo)
library(ShortRead)
library("BSgenome.Scerevisiae.UCSC.sacCer3")
sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3

#----------------------------------------------------------------------------------------------------------------------
# change these variables as necessary------------------------------------------------------------
# put miseq data in fastq/ directory ... rename appropriately (e.g. OS_pool1_R1.fastq.gz)
git.dir='/Users/connie/Dropbox/mutation_rate/20180705miseq/'

#Load accessory functions
source(paste0(git.dir, 'process_targeted_sequencing_accessory_functions.R'))

fastq.dir=paste0(git.dir, 'fastq/')
out.dir=paste0(git.dir, 'out_downsample/') 
dir.create(out.dir)

trimmomatic.call='java -jar /Users/connie/Documents/local-tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 4 ' # trailing space is necessary'
# adapters are here (check these)
adapter.file=paste0(git.dir, 'Nextera_primers_PE.fa')
# key file
sample.key.file=paste0(git.dir, 'site_info.csv')
bwa.call='/Users/connie/Documents/local-tools/bwa-master/bwa mem -t 4'
#pear overlap 5bp (for previous "out" folder, the overlap required is 10bp)
pear.params='-v 5 -m 400 -j 4 -y 8G -q 20'
#-----------------------------------------------------------------------------------------------------------


# build lookup table for samples ---------------------------------------------------------------------------
samples=buildLookupTable(sample.key.file, sacCer3)


# ! Preprocessing steps  ---------------------------------------------------------------------------------

# run pear ----------------------------------------------------------------------------------------------
runPear(samples, pear.params)

# run Trimmomatic
runTrimmomatic(samples, trimmomatic.call, adapter.file)

# organize input by library
samples.by.library=split(samples, samples$library)

# read alignment
alignReads(samples.by.library, out.dir)

# read alignments into R 
sample.alignments=parseAlignments(samples.by.library, samples, out.dir)

# if you want to skip preprocessing, start here  
#saveRDS(sample.alignments, file=paste0(out.dir, 'sample.alignments.RDS'))

# ! -------------------------------------pre-processing End ---------------------------------------------- 

sample.alignments=readRDS(paste0(out.dir, 'sample.alignments.RDS'))

makeSummaryPDF(out.dir, sample.alignments) 

makeAlignmentTxtOut(out.dir, sample.alignments)

consensusMatrices=makeSequenceLogo(out.dir, samples, sample.alignments)
#saveRDS(consensusMatrices, file=paste0(out.dir, 'consensusMatrices.RDS'))
consensusMatrices=readRDS(paste0(out.dir,'consensusMatrices.RDS'))

#make sequence logos
for(libr in names(sample.alignments) ){
  expected.samples=paste(samples$library, samples$samp, samples$chrom, samples$start, samples$end, sep=':')
  for(samp in names(sample.alignments[[libr]]) ) {
    seq.match=match(samp, expected.samples)
    png(file=paste0(out.dir, samp, '_barplot.png'), width=6000, height=500)
    x=consensusMatrices[[libr]][[samp]]
    y=t(t(x)/(colSums(x)))
    bp.coords=barplot(y, col=c('red', 'green', 'yellow', 'blue'),xaxt='n', main=samp, xlab='pos', ylab='frequency', density=25)
    #points(bp.coords, rep(1,ncol(y)), type='n', ylim=c(0,1), xaxt='n', main=samp, xlab='pos', ylab='frequency')
    #ypos=apply(y,2,cumsum)
   # pam.pos=1+(samples[seq.match,'pam']-samples[seq.match, 'start'])
    
    text(bp.coords, y[1,], rownames(y)[1])#, cex=y[1,])
    text(bp.coords, y[2,], rownames(y)[2])#, cex=y[2,])
    text(bp.coords, y[3,], rownames(y)[3])#, cex=y[3,])
    text(bp.coords, y[4,], rownames(y)[4])#, cex=y[4,])
    axis(3, at=bp.coords, labels=colnames(y)) #, col.lab=ifelse(1:ncol(y) %in% pam.pos, 'red', 'black'))
    axis(1, at=bp.coords, labels=1:ncol(y)) #colnames(y))
   # abline(v=pam.pos, col='white', lwd=10)
    dev.off()
    #readline()
  }
}
