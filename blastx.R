pacman::p_load(ggplot2, ggnewscale, GGally, ggtree, tidyr, plyr, dplyr, stringr,ggtree, magrittr, sjmisc, ggseqlogo,gridExtra, seqinr, ape, rentrez, reutils, colorspace, xlsx, openxlsx)
list_sample <- ls()


for(i in 1:length(list_sample)){
  enz<-get(list_sample[i])
  names(enz)<-c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore",
                  "annotation","species","superkingdom","kingdom","phylum","class","order","family","genus","species",
                    "duplicate","unique_read","duplicated_read","total_read")
  assign(list_sample[i], enz)
   enz$duplicate<-gsub("중복", "duplicate", enz$duplicate)
      enz$duplicate<-gsub("단일", "unique", enz$duplicate)
        assign(paste0(list_sample[i],"_filtered"), enz[grep("unique",enz$duplicate),])
          rm(enz)
}

#1.    qseqid    query (e.g., unknown gene) sequence id
#2.    sseqid    subject (e.g., reference genome) sequence id
#3.    pident    percentage of identical matches
#4.    length    alignment length (sequence overlap)
#5.    mismatch    number of mismatches
#6.    gapopen    number of gap openings
#7.    qstart    start of alignment in query
#8.    qend    end of alignment in query
#9.    sstart    start of alignment in subject
#10.    send    end of alignment in subject
#11.    evalue    expect value
#12.    bitscore    bit score