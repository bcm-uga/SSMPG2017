################################################################
##### Modèle discret de diffusion de la frequence allelique ####
################################################################

WFsuivant<-function(X,N,s=0,h=0,u=0,v=0)
{
  
  ## X est un vecteur de réalisations de la diffusion pour lesquelles
  ## il faut calculer les etats suivants (X est la frequence)
  ## N est tel que la population est de taille 2N
  ## s est le parametre de fitness
  ## h est le parametre de heterozygote fitness
  ## u est le taux de mutation de l'allele 1 vers 2
  ## v est le taux de mutation de l'allele 2 vers 1
  
  ## Voir Ewens - Mathematical Population Genetics p26 (1.67)
  
  Xpop<-2*N*X 
  ## population correspondante à la frequence donnée

  varAux<-((1+s)*Xpop^2+(1+s*h)*Xpop*(2*N-Xpop))
  ## variable auxiliaire pour le calcul
  
  probAux<-varAux/(varAux+(1+s*h)*Xpop*(2*N-Xpop)+(2*N-Xpop)^2) 
  ## proba sans mutation
  
  proba<-(1-u)*probAux+(1-probAux)*v
  ## proba finale
  
  return(rbinom(n=length(X),size=2*N,prob=proba)/(2*N))

}

WFtrajectoires<-function(X0,N,n,s=0,h=0,u=0,v=0)
{
  
  ## Genere des trajectoires partant des valeurs de X0
  ## N est tel que la population est de taille 2N
  ## n est le nombre de générations à réaliser
  ## s est le parametre de fitness
  ## h est le parametre de heterozygote fitness
  ## u est le taux de mutation de l'allele 1 vers 2
  ## v est le taux de mutation de l'allele 2 vers 1
  
  ## Voir Ewens - Mathematical Population Genetics p26 (1.67)
  
  ## Renvoie une matrice où chaque colonne est une trajectoire
  
  X<-matrix(nrow=1,ncol=length(X0),X0)
  for (i in 1:n)
  {
    ## print(i)
    X<-rbind(X,WFsuivant(X=X[i,],N=N,s=s,h=h,u=u,v=v))
  }
  return(X)
}
