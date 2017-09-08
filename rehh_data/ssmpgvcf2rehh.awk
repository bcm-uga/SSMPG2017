#gawk -f ssmpgvcf2rehh.awk outputIndAll_pop_sim1a.txt sim1a.vcf
BEGIN{
OFS=ORS=""
}
{
if(FNR==NR){
 if(FNR>1){
  nind[$NF]++
  listind[$NF,nind[$NF]]=$2
  if($NF>npop){npop=$NF}
 }
}else{
 if(FNR==1){
  split(FILENAME,tmp,"\.")
  snpfile=tmp[1]".detsnp"
  prefix=tmp[1]
  }
 if($7=="PASS"){
  hapfileALL=prefix".chr"$1".thap"
  for(i=1;i<=npop;i++){
   hapfile=prefix".chr"$1".pop"i".thap"
   for(j=1;j<=nind[i];j++){
     split($(listind[i,j]+9),tmp,"|")
     zz=(tmp[1]+1)" "(tmp[2]+1)
     if(j==1){line=zz}else{line=line" "zz}   
     if(j==1 && i==1){lineALL=zz}else{lineALL=lineALL" "zz}   
   }
   {print line"\n">hapfile}
 }
   {print line"\n">hapfileALL} 
 nsnp++
 {print "rs"nsnp" "$1" "$2" 1 2\n">prefix".rehhmap.inp"}
 }
 }
}

