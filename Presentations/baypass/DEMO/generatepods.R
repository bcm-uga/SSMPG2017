setwd("~/Bureau/2015_PAPIER_BAYPASS/ANALYSES/CATTLE_DATA/ANA_BAYPASS_COV")
require(MCMCpack) ; require(mvtnorm) ;require(lattice)
source("~/Bureau/NICH_MNV_FORTRAN/utils/nich_mnv_utils.R")
om.bta       = as.matrix(read.table("../ANA_BAYPASS/i_rho1pi_mat_omega.out"))
pibetaparams = read.table("../ANA_BAYPASS/i_rho1pi_summary_beta_params.out",h=T)[,2]
yn.data=geno2YN("geno.dat")

simulate.nichmnv(omega.mat=om.bta,nsnp=100000,sample.size=yn.data$NN,beta.pi=pibetaparams,pi.maf=0.01,remove.fixed.loci=TRUE,suffix="btapods")

