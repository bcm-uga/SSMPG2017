standardize.xtx<-function(xtx,npop){
  ##Cdt1: E(lda*xtx+aa)=npop=lda*E(X) + aa
  ##Cdt2: V(lda*xtx+aa)=2npop=lda^2 V(X)
  lda=sqrt((2*npop)/var(xtx))
  aa=npop - lda*mean(xtx)
  xtx.std=xtx*lda + aa
  ######
  #Wilson–Hilferty transformation
  #####
  #If X~Chi2(ddf=k); \sqrt[3]{X/k} is approximately normally distributed with µ=1−2/9k and var=2/9k
  xtx.std=scale((xtx.std/npop)**(1/3))
  xtx.pval.pos=(1-pnorm(xtx.std))
  xtx.pval.bal=(pnorm(xtx.std))
  list(xtx.std=xtx.std,xtx.pval.pos=xtx.pval.pos,xtx.pval.bal=xtx.pval.bal)
}
