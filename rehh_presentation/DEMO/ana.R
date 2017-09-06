setwd("/home/mathieu/Bureau/PartageWINDOWS/MES_DOCS/EXPOSES/MG/COURS/2017_AUSSOIS/REHH/DEMO")
library(rehh)

code.chr=c(paste0("0",1:9),10:29)
###################
##For each chromosome in turn: 
## i)  Read haplotype data
## ii) Compute iHH_A, iHH_D, iES_Tang and iES_Sabeti
###################

#####
##For the CGU pop.
#####
cnt=0
for(i in 10:15){
 cnt=cnt+1
 tmp.hapfile=bzfile(paste0("EHH_",code.chr[i],"_CGU.hap.bz2"))
 tmp.hap=data2haplohh(hap_file=tmp.hapfile, map_file="map.all.inp", chr.name=i)
 tmp.scan=scan_hh(tmp.hap,threads=4)
 if(cnt==1){wgscan.cgu=tmp.scan}else{wgscan.cgu=rbind(wgscan.cgu,tmp.scan)}
}

head(wgscan.cgu)
#why so many NA's?
i=10
tmp.hap=data2haplohh(bzfile(paste0("EHH_",code.chr[i],"_CGU.hap.bz2")),"map.all.inp",chr.name=i)
calc_ehh(tmp.hap,mrk=8) ;abline(h=0.05)
calc_ehh(tmp.hap,mrk=8,limehh=0.2)
graphics.off()
head(wgscan.cgu,25)

###################
#Computing and plotting iHS within CGU
###################
ihs.cgu=ihh2ihs(wgscan.cgu,minmaf=0.05,freqbin=0.025)

head(ihs.cgu$iHS,25)
ihs.cgu$frequency.class
distribplot(ihs.cgu$iHS$iHS)

ihsplot(ihs.cgu)

###################
#Comparing CGU vs. EUT with rSB and XP-EHH
###################

#####
##Computing iES for for European Taurines
#####
cnt=0
for(i in 10:15){ 
 cnt=cnt+1
 tmp.hap=data2haplohh(bzfile(paste0("EHH_",code.chr[i],"_EUT.hap.bz2")),"map.all.inp",chr.name=i)
 tmp.scan=scan_hh(tmp.hap,threads=4)
 if(cnt==1){wgscan.eut=tmp.scan}else{wgscan.eut=rbind(wgscan.eut,tmp.scan)}
}

####
#Computing Rsb and XpEHH
####

rsb=ies2rsb(wgscan.cgu,wgscan.eut)
xpehh=ies2xpehh(wgscan.cgu,wgscan.eut)

rsbplot(rsb)
xpehhplot(xpehh)

##################
#Vizualising haplotype structure around a core SNP
##################

i=12
tmp.hap=data2haplohh(bzfile(paste0("EHH_",code.chr[i],"_CGU.hap.bz2")),
                     "map.all.inp",chr.name=i)
layout(matrix(1:2,2,1))
bifurcation.diagram(tmp.hap,mrk_foc=456,all_foc=1,nmrk_l=20,nmrk_r=20,
main="Bifurcation diagram (RXFP2 SNP on BTA12): Ancestral Allele")
bifurcation.diagram(tmp.hap,mrk_foc=456,all_foc=2,nmrk_l=20,nmrk_r=20,
main="Bifurcation diagram (RXFP2 SNP on BTA12): Derived Allele")



