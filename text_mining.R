pacman::p_load(dplyr, plyr, tidyr, stringr, stringi, reshape, reshape2, magrittr, rlist, xlsx, openxlsx, data.table, reutils, rentrez) 


{
  `1.11.1.6`<- entrez_search(db="protein", term="1.11.1.6[EC/RN Number]",retmode = "xml", retmax = 99999, use_history = TRUE)
  for( seq_start in seq(1,length(`1.11.1.6`$ids),500)){
    recs <- entrez_fetch(db="protein", web_history=`1.11.1.6`$web_history, rettype="fasta", retmax=500, retstart=seq_start)
    cat(recs, file="1.11.1.6.fasta", append=TRUE)
    cat(seq_start+499, "sequences downloaded\r")
  }}




{
  `gabT`<- entrez_search(db="protein", term="4-aminobutyrate aminotransferase",retmode = "xml", retmax = 99999, use_history = TRUE)
  for( seq_start in seq(1,length(`gabT`$ids),500)){
    recs <- entrez_fetch(db="protein", web_history=`gabT`$web_history, rettype="fasta", retmax=500, retstart=seq_start)
    cat(recs, file="gabT_183215.fasta", append=TRUE)
    cat(seq_start+499, "sequences downloaded\r")
  }}


{
  `phytase`<- entrez_search(db="protein", term="phytase",retmode = "xml", retmax = 99999, use_history = TRUE)
  for( seq_start in seq(1,length(`phytase`$ids),500)){
    recs <- entrez_fetch(db="protein", web_history=`phytase`$web_history, rettype="fasta", retmax=500, retstart=seq_start)
    cat(recs, file="phytase_86952.fasta", append=TRUE)
    cat(seq_start+499, "sequences downloaded\r")
  }}



{
  `wecB`<- entrez_search(db="protein", term="UDP-N-acetylglucosamine 2-epimerase",retmode = "xml", retmax = 999999, use_history = TRUE)
  for( seq_start in seq(1,length(`wecB`$ids),500)){
    recs <- entrez_fetch(db="protein", web_history=`wecB`$web_history, rettype="fasta", retmax=500, retstart=seq_start)
    cat(recs, file="wecB_189716_new.fasta", append=TRUE)
    cat(seq_start+499, "sequences downloaded\r")
  }}
{
  `tpiA`<- entrez_search(db="protein", term="Triosephosphate isomerase",retmode = "xml", retmax = 999999, use_history = TRUE)
  for(seq_start in seq(127001,length(`tpiA`$ids),500)){
    recs <- entrez_fetch(db="protein", web_history=`tpiA`$web_history, rettype="fasta", retmax=500, retstart=seq_start)
    cat(recs, file="tpiA_127001-214471.fasta", append=TRUE)
    cat(seq_start+499, "sequences downloaded\r")
  }}

{
  `eglS`<- entrez_search(db="protein", term="glutamyl aminopeptidase",retmode = "xml", retmax = 999999, use_history = TRUE)
  for(seq_start in seq(1,length(`eglS`$ids),500)){
    recs <- entrez_fetch(db="protein", web_history=`eglS`$web_history, rettype="fasta", retmax=500, retstart=seq_start)
    cat(recs, file="ytoP_77985.fasta", append=TRUE)
    cat(seq_start+499, "sequences downloaded\r")
  }}

