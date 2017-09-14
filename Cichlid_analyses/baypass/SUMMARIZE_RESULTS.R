source("standardize.xtx.R")
source("baypass_utils.R")

detsnp=read.table("cichlid.detsnp.filtmaf1p")
tmp=1:3 ; names(tmp)=unique(detsnp[,1])
chr.col=tmp[detsnp[,1]]

om=as.matrix(read.table("anacovis1_mat_omega.out"))

for(i in 1:3){
  tmp.xtx=read.table(paste0("anacovis",i,"_summary_pi_xtx.out"),h=T)$M_XtX
  tmp.bfis=read.table(paste0("anacovis",i,"_summary_betai_reg.out"),h=T)$BF.dB.
  if(i==1){
    xtx=tmp.xtx ; bfis=tmp.bfis
  }else{
    xtx=cbind(xtx,tmp.xtx)
    bfis=cbind(bfis,tmp.bfis)
  }
}

cor(xtx)
cor(bfis)

###Generate PODs
pi.beta.coef=read.table("anacovis1_summary_beta_params.out",h=T)$Mean
geno.data<-geno2YN("cichlid.genobaypass.filtmaf1p")
simu<-simulate.baypass(omega.mat=om,nsnp=100000,sample.size=geno.data$NN,beta.pi=pi.beta.coef,pi.maf=0,suffix="cichlid.pods")

xtx.pods=read.table(paste0("anapod_summary_pi_xtx.out"),h=T)$M_XtX
bfis.pods=read.table(paste0("anapod_summary_betai_reg.out"),h=T)$BF.dB.
om.pods=as.matrix(read.table("anapod_mat_omega.out"))
plot(om.pods,om) ; abline(a=0,b=1)
fmd.dist(om.pods,om)
####################
bfmc=read.table(paste0("anaaux_summary_betai.out"),h=T)$BF.dB.
pip.is=read.table(paste0("anaauxisb1_summary_betai.out"),h=T)$M_Delta

layout(matrix(1:4,4,1))
plot(xtx[,1],col=chr.col,pch=16,main="XtX")
abline(h=quantile(xtx.pods,probs=0.999),lty=2)
plot(bfis[,1],col=chr.col,pch=16,main="BFis (asso. with scaled fish. press.)")
abline(h=20,lty=2)
plot(bfmc,col=chr.col,pch=16,main="BFmc (asso. with scaled fish. press.)")
abline(h=20,lty=2)
plot(pip.is,col=chr.col,pch=16,main="PIP (Ising model to account for spatial dep. of SNPs")
abline(h=0.5,lty=2)

##################
contrasts=read.table("anacontrasts_summary_contrast.out",h=T)
cPOP1_vs_ALL=contrasts[contrasts$CONTRAST==1,"C_std"]
cPOP4_vs_ALL=contrasts[contrasts$CONTRAST==2,"C_std"]
cPOP2_vs_ALL=contrasts[contrasts$CONTRAST==3,"C_std"]

layout(matrix(1:3,3,1))
Y_LIM=c(min(c(cPOP1_vs_ALL,cPOP2_vs_ALL,cPOP4_vs_ALL)),max(c(cPOP1_vs_ALL,cPOP2_vs_ALL,cPOP4_vs_ALL)))
plot((cPOP1_vs_ALL),main="POP1 vs. (POP2,POP4)",ylim=Y_LIM,col=chr.col,pch=16)
abline(h=-3,lty=2) ;abline(h=3,lty=2) 
plot((cPOP2_vs_ALL),main="POP2 vs. (POP1,POP4)",ylim=Y_LIM,col=chr.col,pch=16)
abline(h=-3,lty=2) ;abline(h=3,lty=2) 
plot((cPOP4_vs_ALL),main="POP4 vs. (POP1,POP2)",ylim=Y_LIM,col=chr.col,pch=16)
abline(h=-3,lty=2) ;abline(h=3,lty=2) 
