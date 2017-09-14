setwd("cichlid_data")
library(rehh)

code.chr=c(237,262,263)

cnt=0
for(i in 1:length(code.chr)){
  cnt=cnt+1
  tmp.hap.file=paste0("../cichlid.Contig",code.chr[i],".pop1.thap")
  tmp.data=data2haplohh(hap_file=tmp.hap.file,map_file="../cichlid.rehhmap.inp",min_maf = 0.01,chr.name=code.chr[i],haplotype.in.columns=T)
  tmp.pop1.scan=scan_hh(tmp.data,maxgap = 5000)
  tmp.hap.file=paste0("../cichlid.Contig",code.chr[i],".pop2.thap")
  tmp.data=data2haplohh(hap_file=tmp.hap.file,map_file="../cichlid.rehhmap.inp",min_maf = 0.01,chr.name=code.chr[i],haplotype.in.columns=T)
  tmp.pop2.scan=scan_hh(tmp.data,maxgap = 5000) 
  tmp.hap.file=paste0("../cichlid.Contig",code.chr[i],".pop4.thap")
  tmp.data=data2haplohh(hap_file=tmp.hap.file,map_file="../cichlid.rehhmap.inp",min_maf = 0.01,chr.name=code.chr[i],haplotype.in.columns=T)
  tmp.pop4.scan=scan_hh(tmp.data,maxgap = 5000)  
  if(cnt==1){
    scan.pop1=tmp.pop1.scan ;   scan.pop2=tmp.pop2.scan ;   scan.pop4=tmp.pop4.scan
  }else{
    scan.pop1=rbind(scan.pop1,tmp.pop1.scan)
    scan.pop2=rbind(scan.pop2,tmp.pop2.scan)
    scan.pop4=rbind(scan.pop4,tmp.pop4.scan)    
  }
}
############
###########
ihs.pop1=ihh2ihs(scan.pop1)
ihs.pop2=ihh2ihs(scan.pop2)
ihs.pop4=ihh2ihs(scan.pop4)

##############
snp.sel=rownames(scan.pop1)[rownames(scan.pop1) %in% rownames(scan.pop2)]
xpehh.pop1_vs_pop2=ies2xpehh(scan.pop1[snp.sel,],scan.pop2[snp.sel,])
snp.sel=rownames(scan.pop1)[rownames(scan.pop1) %in% rownames(scan.pop4)]
xpehh.pop1_vs_pop4=ies2xpehh(scan.pop1[snp.sel,],scan.pop4[snp.sel,])
snp.sel=rownames(scan.pop2)[rownames(scan.pop2) %in% rownames(scan.pop4)]
xpehh.pop2_vs_pop4=ies2xpehh(scan.pop2[snp.sel,],scan.pop4[snp.sel,])






