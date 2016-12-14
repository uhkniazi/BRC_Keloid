# Name: spikein_bam_check.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 28/11/2016
# Desc: Check number of reads aligned to spikeins

## set variables and source libraries
gcswd = getwd()
setwd('Keloid_main/')
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CBamQuality/master/CBamQuality.R'
download(url, 'CBamQuality.R')

# load the required packages
source('CBamQuality.R')
# delete the file after source
unlink('CBamQuality.R')

csPath = '/run/user/1000/gvfs/sftp:host=10.202.64.28,user=k1625253/users/k1625253/brc_scratch/Data/ProjectsData/BRC_Keloid/Aligned/S021/K1-1st_S021_R1_.fastq.gz_q10_sort_rd.bam'
# load bam file
bf = BamFile(csPath)

# get names of sequences
sn = seqinfo(bf)
gr = as(sn, 'GRanges')
seqlengths(bf)

ob = CBamScaffold(bf, 'ERCC-00003')

param = ScanBamParam(what=scanBamWhat(), which=gr['ERCC-00003'])
oGA = readGAlignments(bf, param = param)

# get the dna string set object with sequences
ds = mcols(oGA)$seq
ds = DNAStringSet(c(oRefSeq['ERCC-00003'], ds))
# read reference sequences
oRefSeq = readDNAStringSet('Data_external/ERCC92.fa')

oPWA = pairwiseAlignment(pattern = ds, subject = oRefSeq['ERCC-00003'], type='local')

# NAME: f_oDNAStringSetConvertPWAMatrix
# ARGS: pwa.matrix = a pairwise alignment object in matrix form, generated via a call like
#       as.matrix(pairwiseAlignment(patten, subject))
# DESC: converts the sequences in the matrix to a DNAStringSet object which can be exported to a FASTA
#       file to be viewed in a viewer of choice e.g. jalview
# RETS: a DNAStringSet object with all sequences from the pwa matrix
f_oDNAStringSetConvertPWAMatrix= function(pwa.matrix){
  require(Biostrings)
  # sequence to return
  seq.export = DNAStringSet()
  # extract each sequence from the matrix
  for (i in 1:nrow(pwa.matrix)){
    s = DNAStringSet(pwa.matrix[i,])
    s = DNAStringSet(unlist(s))
    seq.export = append(seq.export, s)
  }
  # set names for sequences
  names(seq.export) = rownames(pwa.matrix)
  return(seq.export)
}

seq = f_oDNAStringSetConvertPWAMatrix(as.matrix(oPWA))

writeXStringSet(seq, filepath = 'Results/ERCC-00003.fa')

### check counting agasint bam directly
library(GenomicAlignments)
library(rtracklayer)

# load the ERCC gtf with spikein information
oGRercc = import('~/Data/MetaData/ERCC92.gtf')

# reformat metadata column
f = oGRercc$gene_id
df = DataFrame(exon_id=oGRercc$transcript_id, exon_name=NA)
mcols(oGRercc) = df
oGRLercc = split(oGRercc, f)

mCounts.ercc = assays(summarizeOverlaps(oGRLercc, bf, ignore.strand = F, singleEnd=F))$counts