#### Motifbreak for all SNPs
library(BSgenome)
library(motifbreakR)
library(GenomicFeatures)
library(GenomicRanges)
library(dplyr)
library(data.table)
library(magrittr)
library(IRanges)
library(rtracklayer)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)

set.seed(20241021)
options(width=200,scipen=999)

project.dir = '/vast/scratch/users/vespasiani.d/gwas-ai-snps'
out.dir = paste(project.dir,'out',sep='/')
data.dir = paste(project.dir,'data',sep='/')
dir.create(project.dir,showWarnings=F)
setwd(project.dir)
dir.create(out.dir,showWarnings=F)
dir.create(data.dir,showWarnings=F)

##----------------------------------------------------------------
## Read input GWAS hits associated with AI
gwas.hits = fread(paste(data.dir,'AI_SNPs_unique_set_with_coords.txt',sep='/'))%>%setorderv('SNP')

##----------------------------------------------------------------
## Data wrangling to meet motifbreakR input file format requirements (check Details section of the snps.from.file function)
## This is a 4-column tab-separated bed file containing: seqnames, start, end, name (a variant ID being specified as seqnames:variantPosition:ref:alt ... see below)
## Each variant needs to have its own entry in the file, so if a variant has multiple alt alleles, these need to have each own row in the file (ie, 1 variant-ref-alt combination per row)

## Another important point in file formatting for motifbreakR is to make sure to have the variant coordinates specified correctly. 
## The start position of each variant needs to be in the same coordinates of those reported by NCBI. 
## For the set of SNPs contained in the AI_SNPs_unique_set_with_coords.txt file, these need to be -1 shifted.
## For indels, the file formatting requires a bit extra manipulation, so check out the notes on the relative sections to see how I have manipulated them

hits.motifbreakR <- copy(gwas.hits)[, ref := tstrsplit(allele,'/',fixed=T)[[1]]] ## extract the ref allele from the allele column in the file
hits.motifbreakR$alts <- stringr::str_remove(hits.motifbreakR$allele, paste(hits.motifbreakR$ref,'/',sep='')) ## remove ref allele from allele column 
hits.motifbreakR$numbAlt <- stringr::str_count(hits.motifbreakR$alts, "/") ## count how many alternative alleles are reported for each variant
hits.motifbreakR$numbAlt <- ifelse(hits.motifbreakR$numbAlt==0,hits.motifbreakR$numbAlt,hits.motifbreakR$numbAlt+1) ## this is just an internal command to then duplicate the rows based on the number of alt alleles each variant has
hits.motifbreakR <- hits.motifbreakR[,cbind(.SD,dup=1:numbAlt),by="SNP"] ## duplicate each variant by the number of alt alleles it has
hits.motifbreakR <- lapply(split(hits.motifbreakR,by='SNP'),function(s){
    if(sum(s$dup)==1){
        s$dup = 0
        s = unique(s)
        s$alt <- s$alts
    }else{
        alleles <-  data.table(
            alt = sapply(strsplit(s$alts, '/'), `[`, s$dup)[,1]
            )
        s$alt <- alleles$alt
    }
    s <- s[,c('numbAlt','dup','alts'):=NULL]
    return(s)
})%>%rbindlist()

##----------------------------------------------------------------
## Separate indels from SNPs and process variant type separately
snps <- copy(hits.motifbreakR)[nchar(ref)==1][nchar(alt)==1][ ! ref %like% '-'][ ! alt %like% '-'] ## n = 5516
indels <- copy(hits.motifbreakR)[ ! SNP %in% snps$SNP] ## n = 139

##----------------------------------------------------------------
## Manipulate SNPs
## For SNPs, the only requirement is to have the end column reporting the position of the actual variant.
## this means setting the end column as the current start and then subtracting 1 from the current start column
snps.motifbreakR <- copy(hits.motifbreakR)[nchar(ref)==1][nchar(alt)==1][ ! ref %like% '-'][ ! alt %like% '-']
snps.motifbreakR <- snps.motifbreakR[
    ,seqnames := paste('chr',chr,sep='')
    ][
        ,start := chrom_start-1
        ][
            ,end:=chrom_start
            ][
                ,name := paste(seqnames,end,ref,alt,sep=':') ## end here refers to the position of the actual variant
                ][
                    ,c('seqnames','start','end','name')
]

##----------------------------------------------------------------
## Manipulate deletions
## In case of deletions your start and end coordinates must span the alleles being deleted. You cannot have start = end (unless its a single nt deletion, see below), as is currently reported in the file
## Instead you need to specify the entire region affected by the deletion (see rs541631272 example in https://bioconductor.org/packages/release/bioc/vignettes/motifbreakR/inst/doc/motifbreakR-vignette.html)
## So for example, consider this multi-nucleotide deletion as an example: chr11:65736582-65736582  ref= TAA  alt = T
## The deleted nucleotides are 65736583 = A and 65736584 = A and 65736582 = T is actual nucleotide that remains, 
## The entire region affected by the deletion hence is chr11:65736582-65736584
## However, because the above requirement of a start-1 position is also valid for deletions, the final region to be specified in the bed file will be chr11:65736581-65736584

deletions.multinucleotide <- copy(indels)[nchar(ref) > nchar(alt)][ ! alt %like% '-'] ## this is to retain deletions from the indels subset (ie, entries with the number of nt in the alt column is < than those in the ref column)

deletions.multinucleotide.motifbreakR <- copy(deletions.multinucleotide)[
    ,width_deletion := nchar(ref)-1 ## calculate the width of the region affected by the deletion
        ][
            ,start := chrom_start-1 ## offset start coordinated by 1 as above
            ][
                ,end := chrom_start + width_deletion ## specify the coordinates of the region affected by the deletion, chrom_start in the original file is the beginning of the deleted region
                ][
                    ,seqnames := paste('chr',chr,sep='')
                    ][
                        ,name := paste(seqnames,end,ref,alt,sep=':') 
                        ][
                    ,c('seqnames','start','end','name')
]

# getSeq(BSgenome.Hsapiens.UCSC.hg38, 'chr6', start=135422619, end=135422627) ## check whether the specified region (with start +1 if you run this command) contains the same nucleotides as those in the ref column

## Single nucleotide deletions (ie, all the entries where there is no alt allele in the data)
## The solution to manipulate these entries is to add to the reference allele the nucleotide located just downstream of it (ie, chrom_start + 1)
## so consider this single nucleotide deletion as an example: chr1:116738074-116738074, ref = C alt = - ; instead of having alt = - , you need to code the deletion as chr1:116738073-116738074, ref = CT alt = T
deletions.singlenucleotide <-  copy(indels)[ alt %like% '-']
deletions.singlenucleotide$seqnames <- paste('chr',deletions.singlenucleotide$chr,sep='')
deletions.singlenucleotide <- deletions.singlenucleotide[,ref_new := as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, deletions.singlenucleotide$seqnames, start= deletions.singlenucleotide$chrom_start, end= deletions.singlenucleotide$chrom_end+1))]
deletions.singlenucleotide <- deletions.singlenucleotide[,alt := stringr::str_remove(ref_new,ref)]## now make the alt allele by deleting the ref nucleotide from the ref sequence
deletions.singlenucleotide <- deletions.singlenucleotide[,ref:= ref_new][,c('ref_new'):=NULL]

## Given now I've recoded the alleles for single-nt deletions, the lines below are the same as those used above for making the deletions.multinucleotide.motifbreakR table 
deletions.singlenucleotide.motifbreakR <- copy(deletions.singlenucleotide)[
    ,width_deletion := nchar(ref)-1 ## calculate the width of the region affected by the deletion
        ][
            ,start := chrom_start-1 ## offset start coordinated by 1 as above
            ][
                ,end := chrom_start + width_deletion ## specify the coordinates of the region affected by the deletion, chrom_start in the original file is the beginning of the deleted region
                ][
                    ,name := paste(seqnames,end,ref,alt,sep=':') 
                        ][
                    ,c('seqnames','start','end','name')
]

deletions.motifbreakR <- rbind(deletions.singlenucleotide.motifbreakR,deletions.multinucleotide.motifbreakR)

##----------------------------------------------------------------
## Multinucleotide insertions
## Conversely to what done for the deletions, the start and end coordinates of the table below have to reflect just the width of the region being affected by the insertion in the ref column
## So if ref is just a nucleotide then then width of the region must be 1. So this equation must be true chrom_end-chrom_start + 1 = nchar(ref) ---I've checked it and it is indeed true
insertions <- copy(indels)[ ! SNP %in% c(deletions.singlenucleotide$SNP,deletions.multinucleotide$SNP)]
insertions.multinucleotide <- copy(insertions)[ref !='-'][,seqnames := paste('chr',chr,sep='')]

insertions.multinucleotide.motifbreakR <- copy(insertions.multinucleotide)[
    ,c('SNP','seqnames','chrom_start','chrom_end','ref','alt')
    ][
        ,start := chrom_start-1
            ][
                ,end := chrom_end
                ][
                    ,name := paste(seqnames,end,ref,alt,sep=':') 
                        ][
                    ,c('seqnames','start','end','name')
]

## Single nucleotide insertions (ie, those lacking the reference allele)
## The solution is specular to that applied to the single-nucleotide deletions, ie add to the reference allele the nucleotide located just upstread of it (ie, chrom_start - 1) and then add to the current alt allele the ref one as well
## Example: chr2:102449062-102449062, ref = - alt = G ; you need to code the insertions as chr2:102449061-102449062, ref = A alt = AG 

# getSeq(BSgenome.Hsapiens.UCSC.hg38, 'chr2', start=102449062, end=102449063)

insertions.singlenucleotide <- copy(insertions)[ref =='-'][,seqnames := paste('chr',chr,sep='')]
insertions.singlenucleotide <- insertions.singlenucleotide[,ref := as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, insertions.singlenucleotide$seqnames, start= insertions.singlenucleotide$chrom_end, end=insertions.singlenucleotide$chrom_end))] ## note that chrom_end is specified in both arguments as this is the position of the insertion and chrom_start = chrom_end +1 
insertions.singlenucleotide <- insertions.singlenucleotide[,alt := paste(ref,alt,sep='')]

insertions.singlenucleotide.motifbreakR <- copy(insertions.singlenucleotide)[
    ,c('SNP','seqnames','chrom_start','chrom_end','ref','alt')
    ][
        ,start := chrom_end-1
            ][
                ,end := chrom_end
                ][
                    ,name := paste(seqnames,end,ref,alt,sep=':') ##chromosome:start:REF:ALT
                        ][
                    ,c('seqnames','start','end','name')
]


insertions.motifbreakR <- rbind(insertions.singlenucleotide.motifbreakR,insertions.multinucleotide.motifbreakR)

##----------------------------------------------------------------
## Double check that all entries are correct by inspecting that:
## 1) the REF allele (in name) corresponds to the actual ref allele in the BSgenome 
## 2) the ref and the alt alleles are the correct ones given the list of SNPs received
## 3) the position 

##----------------------------------------------------------------
## export file as bed with no header to then be read back in using the motifbreakR snps.from.file() function

fwrite(
    rbind(deletions.motifbreakR,insertions.motifbreakR,snps.motifbreakR),
    file = paste(data.dir,'motifbreakR-input-AI-GWAS-Variants.bed',sep='/'),col.names = F, row.names = F, sep = "\t", quote = F
)

##----------------------------------------------------------------
## read back GWAS hits into motifbreakR input format

motifbreak.input = snps.from.file(
    paste(data.dir,'motifbreakR-input-AI-GWAS-Variants.bed',sep='/'),
    search.genome = BSgenome.Hsapiens.UCSC.hg38, 
    indels = TRUE,
    format = "bed"
)

##----------------------------------------------------------------
## select the PWMs 
hocomoco.cores = c(
  'HOCOMOCOv11-core-A','HOCOMOCOv11-core-B','HOCOMOCOv11-core-C',
  'HOCOMOCOv11-secondary-A','HOCOMOCOv11-secondary-B','HOCOMOCOv11-secondary-C'
)
hocomocov11 = subset(MotifDb, organism=='Hsapiens' & dataSource %in% hocomoco.cores)

jaspar2022 = subset(MotifDb, organism=='Hsapiens' & dataSource=='jaspar2022')

##----------------------------------------------------------------
## Find broken motifs using Jaspar2022 list of PWMs
motifsBroken.jaspar2022 <- motifbreakR(
    snpList = motifbreak.input, 
    filterp = TRUE,
    pwmList = jaspar2022,
    threshold = 1e-4,
    method = "default",
    verbose = TRUE,
    bkg = c(A=0.3, C=0.2, G=0.2, T=0.3), ## use nt frequencies for non-coding regions human genome
    BPPARAM = MulticoreParam(workers=10)
)%>%as.data.table()

# fwrite(motifsBroken.jaspar2022,file = paste(out.dir,'motifbreakR-out-jaspar2022-AI-GWAS-Variants.txt',sep='/'),col.names = T, row.names = F, sep = "\t", quote = F)

## Find broken motifs using HOCOMOCO v11 list of PWMs
motifsBroken.hocomocov11 <- motifbreakR(
    snpList = motifbreak.input, 
    filterp = TRUE,
    pwmList = hocomocov11,
    threshold = 1e-4,
    method = "default",
    verbose = TRUE,
    bkg = c(A=0.3, C=0.2, G=0.2, T=0.3), ## use nt frequencies for non-coding regions human genome
    BPPARAM = MulticoreParam(workers=10)
)%>%as.data.table()

fwrite(motifsBroken.hocomocov11,file = paste(out.dir,'motifbreakR-out-hocomocov11-AI-GWAS-Variants.txt',sep='/'),col.names = T, row.names = F, sep = "\t", quote = F)

##----------------------------------------------------------------
## Add rsID info 
hits.snpID <- copy(hits.motifbreakR)[,seqnames := paste('chr',chr,sep='')][,c('allele','ref','alt','chr'):=NULL]

## HOCOMOCO v11
hocomoco.out <- fread(paste(out.dir,'motifbreakR-out-hocomocov11-AI-GWAS-Variants.txt',sep='/'))
hocomoco.snpID <- copy(hocomoco.out)[,c('SNP_id','seqnames','start','end','REF','ALT','varType')]%>%unique()

## Insertions
insertions.hocomoco <- copy(hocomoco.snpID)[varType=='Insertion']
insertions.rsid.chromStart <- copy(insertions.hocomoco)[,chrom_start := start +1][hits.snpID,on=c('seqnames','chrom_start'),nomatch=0,allow.cartesian=T][,c('SNP_id','SNP')]%>%setnames('SNP','rsid')%>%unique()
insertions.rsid.chromEnd <- copy(insertions.hocomoco)[,chrom_end := end][hits.snpID,on=c('seqnames','chrom_end'),nomatch=0,allow.cartesian=T][,c('SNP_id','SNP')]%>%setnames('SNP','rsid')%>%unique()
insertions.rsid <- rbind(insertions.rsid.chromStart,insertions.rsid.chromEnd)
insertions.hocomoco <-insertions.hocomoco[insertions.rsid,on='SNP_id',nomatch=0,allow.cartesian=T]%>%setorderv('rsid')

## Deletions
deletions.hocomoco <- copy(hocomoco.snpID)[varType=='Deletion']
deletions.rsid <- copy(deletions.hocomoco)[,chrom_start := start][hits.snpID,on=c('seqnames','chrom_start'),nomatch=0,allow.cartesian=T][,c('SNP_id','SNP')]%>%setnames('SNP','rsid')%>%unique()
deletions.hocomoco <-deletions.hocomoco[deletions.rsid,on='SNP_id',nomatch=0,allow.cartesian=T]%>%setorderv('rsid')

## SNPs
snps.hocomoco <- copy(hocomoco.snpID)[varType=='SNV']
snps.rsid <- copy(snps.hocomoco)[,chrom_start := start][hits.snpID,on=c('seqnames','chrom_start'),nomatch=0,allow.cartesian=T][,c('SNP_id','SNP')]%>%setnames('SNP','rsid')%>%unique()
snps.hocomoco <- snps.hocomoco[snps.rsid,on='SNP_id',nomatch=0,allow.cartesian=T]%>%setorderv('rsid')

## MNV 
mnv.hocomoco <- copy(hocomoco.snpID)[varType=='Other']
mnv.rsid <- copy(mnv.hocomoco)[,chrom_start := start][hits.snpID,on=c('seqnames','chrom_start'),nomatch=0,allow.cartesian=T][,c('SNP_id','SNP')]%>%setnames('SNP','rsid')%>%unique()
mnv.hocomoco <- mnv.hocomoco[mnv.rsid,on='SNP_id',nomatch=0,allow.cartesian=T]%>%setorderv('rsid')

hocomoco.rsid <- rbind(snps.hocomoco,deletions.hocomoco,insertions.hocomoco,mnv.hocomoco)[,c('SNP_id','rsid','REF','ALT','varType')]

## combine all information into final table
hocomoco.final <- hocomoco.out[hocomoco.rsid,on=c('SNP_id','REF','ALT','varType'),nomatch=0,allow.cartesian=T]

## check that you've reained all variants originally identified by hocomoco
hocomoco.original = fread(paste(out.dir,'motifbreakR-out-hocomocov11-AI-GWAS-Variants.txt',sep='/'))
setdiff(hocomoco.original$SNP_id,hocomoco.final$SNP_id)
setdiff(hocomoco.final$SNP_id,hocomoco.original$SNP_id)


##----------------------------------------------------------------
## JASPAR2022
jaspar.out <- fread(paste(out.dir,'motifbreakR-out-jaspar2022-AI-GWAS-Variants.txt',sep='/'))
jaspar.snpID <- copy(jaspar.out)[,c('SNP_id','seqnames','start','end','REF','ALT','varType')]%>%unique()

## Insertions
insertions.jaspar <- copy(jaspar.snpID)[varType=='Insertion']
insertions.rsid.chromStart <- copy(insertions.jaspar)[,chrom_start := start +1][hits.snpID,on=c('seqnames','chrom_start'),nomatch=0,allow.cartesian=T][,c('SNP_id','SNP')]%>%setnames('SNP','rsid')%>%unique()
insertions.rsid.chromEnd <- copy(insertions.jaspar)[,chrom_end := end][hits.snpID,on=c('seqnames','chrom_end'),nomatch=0,allow.cartesian=T][,c('SNP_id','SNP')]%>%setnames('SNP','rsid')%>%unique()
insertions.rsid <- rbind(insertions.rsid.chromStart,insertions.rsid.chromEnd)
insertions.jaspar <-insertions.jaspar[insertions.rsid,on='SNP_id',nomatch=0,allow.cartesian=T]%>%setorderv('rsid')

## Deletions
deletions.jaspar <- copy(jaspar.snpID)[varType=='Deletion']
deletions.rsid <- copy(deletions.jaspar)[,chrom_start := start][hits.snpID,on=c('seqnames','chrom_start'),nomatch=0,allow.cartesian=T][,c('SNP_id','SNP')]%>%setnames('SNP','rsid')%>%unique()
deletions.jaspar <-deletions.jaspar[deletions.rsid,on='SNP_id',nomatch=0,allow.cartesian=T]%>%setorderv('rsid')

## SNPs
snps.jaspar <- copy(jaspar.snpID)[varType=='SNV']
snps.rsid <- copy(snps.jaspar)[,chrom_start := start][hits.snpID,on=c('seqnames','chrom_start'),nomatch=0,allow.cartesian=T][,c('SNP_id','SNP')]%>%setnames('SNP','rsid')%>%unique()
snps.jaspar <- snps.jaspar[snps.rsid,on='SNP_id',nomatch=0,allow.cartesian=T]%>%setorderv('rsid')

## MNV
mnv.jaspar <- copy(jaspar.snpID)[varType=='Other']
mnv.rsid <- copy(mnv.jaspar)[,chrom_start := start][hits.snpID,on=c('seqnames','chrom_start'),nomatch=0,allow.cartesian=T][,c('SNP_id','SNP')]%>%setnames('SNP','rsid')%>%unique()
mnv.jaspar <- mnv.jaspar[mnv.rsid,on='SNP_id',nomatch=0,allow.cartesian=T]%>%setorderv('rsid')

jaspar.rsid <- rbind(snps.jaspar,deletions.jaspar,insertions.jaspar,mnv.jaspar)[,c('SNP_id','rsid','REF','ALT','varType')]

## combine all information into final table
jaspar.final <- jaspar.out[jaspar.rsid,on=c('SNP_id','REF','ALT','varType'),nomatch=0,allow.cartesian=T]

## check that you've reained all variants originally identified by Jaspar
jaspar.original = fread(paste(out.dir,'motifbreakR-out-jaspar2022-AI-GWAS-Variants.txt',sep='/'))
setdiff(jaspar.original$SNP_id,jaspar.final$SNP_id)
setdiff(jaspar.final$SNP_id,jaspar.original$SNP_id)


##----------------------------------------------------------------
## Export final results
fwrite(jaspar.final,file = paste(out.dir,'AI-GWAS-Variants-motifbreakR-jaspar2022.txt',sep='/'),col.names = T, row.names = F, sep = "\t", quote = F)
fwrite(hocomoco.final,file = paste(out.dir,'AI-GWAS-Variants-motifbreakR-hocomocov11.txt',sep='/'),col.names = T, row.names = F, sep = "\t", quote = F)

length(unique(jaspar.final$rsid))
length(unique(hocomoco.final$rsid))
length(unique(hits.motifbreakR$rsid))
