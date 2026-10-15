hu <- load("/Users/henno/Documents/Skripte/R_Scripts/DataAnalysis/20210129_AK_ProteinIsoforms/20250821_Rebuttal/GithubUpdate/Proteomics/Fasta_Database/compiled_gff_20260224_UniProt_version_2022_v4.rda")
hu

class(gff)
table(gff$DatabaseType)
gff$UniProtType

gffNew <- gff[Msequence!="Peptide",{
  
  temp <<- .SD
  gr <<- .BY
  if(dim(temp)[1]>1){
    if(any(UniProtType=="sp")){
      temp <- temp[UniProtType=="sp"]
    }else{
      if(any(UniProtType=="tr")){
        temp <- temp[UniProtType=="tr"]
      }else{
        stop()
      }
    }
  }
  if(dim(temp)[1]>1){
    stop()
  }
temp
  
  
  #
},Msequence]

unlink(iunlink)

gffNew[,{
  cat("\r",.GRP)
  write(paste(">",.BY$ProteinID,sep=""),"compiled_gff_20260224_UniProt_version_2022_v4.fasta",append=T)
  write(.BY$Msequence,"compiled_gff_20260224_UniProt_version_2022_v4.fasta",append = T)
  NULL
},.(Msequence,ProteinID)]
