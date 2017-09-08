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

WFexempleTrajectoires<-function(X0=c(0.1,0.5,0.8),n=10000,N=10000)
{

  ## n est le nombre de générations à simuler
  ## X0 est le vecteur contenant les points de départ des simulations
  
  plot(0,xlim=c(0,n),ylim=c(0,1),type="n")
  p<-WFtrajectoires(X0,N,n)
  for (j in 1:length(X0))
  {
    points(0:n,p[,j],type="l",col=rainbow(length(X0))[j])
  }
}

WFsimuleTempsAbso<-function(p,n,N,x=c(0,1),s=0,h=0,u=0,v=0,norm=FALSE)
{
  
  ## p est la valeur initiale de la chaine
  ## n est le nombre de temps à simuler
  ## x est l'ensemble des etats absorbants
  ## N le facteur de taille de population
  ## s est le parametre de fitness
  ## h est le parametre de heterozygote fitness
  ## u est le taux de mutation de l'allele 1 vers 2
  ## v est le taux de mutation de l'allele 2 vers 1
  
  ## Voir Ewens - Mathematical Population Genetics p26 (1.67)
  
  
  t<-NULL
  
  for (i in 1:n)
  {
    y<-p
    compteur<-0
    while (!is.element(y,x))
    {
      compteur<-compteur+1
      y<-WFsuivant(y,N=N,s=s,h=h,u=u,v=v)
    }
    #print(c(i,compteur))
    t<-c(t,compteur)
  }
  if (norm)
  {
    return(t/(2*N))
  }
  else
  {
    return(t)
  }
}

#n<-1000
#N<-10000
#Donnees<-WFestimeTempsAbso(0.1,n=n,N=N)

quelquesDensites<-function(X0=c(0.1,0.3,0.5),n=100,N=10000)
{
  plot(0,0,type="n",xlab="temps",ylab="proba d'absorption",main="Quelques densités du temps d'absorption",xlim=c(-1,8),ylim=c(0,1.5))
  i<-0
  for (p in X0)
  {
    i<-i+1
    Donnees<-WFsimuleTempsAbso(p,n=n,N=N)
    lines(density(Donnees/(2*N),bw=(1/n)^0.2,adjust=1,kernel="epanechnikov",from=0,to=8),col=rainbow(length(X0))[i]);
  }
}

quelquesFRE<-function(X0=c(0.1,0.3,0.5),n=100,N=10000)
{
  plot(0,0,type="n",xlab="temps",ylab="F(T<t)",main="Quelques FRE du temps d'absorption",xlim=c(-1,8),ylim=c(0,1))
  i<-0
  for (p in X0)
  {
    i<-i+1
    Donnees<-WFsimuleTempsAbso(p,n=n,N=N)
    points(c(0,sort(Donnees/(2*N)),10^99),seq(from=0,to=1,by=1/(length(Donnees)+1)),type="l",col=rainbow(length(X0))[i]);
  }
}

WFcalculeTempsAbso<-function(p,N,s=0,h=0,u=0,v=0)
{
  # p est un vecteur. Les temps d'absorbtion sont calculés 
  # pour des trajectoires partantes des éléments de p
  
  ## Voir Ewens - Mathematical Population Genetics p137 (4.3-4.5) et p141 (4.21)
  
  a<-function(x)
  {
    # a est la fonction de shift
    return(2*N*s*x*(1-x)*(x+h*(1-2*x))-2*N*u*x+2*N*v*(1-x))
  }
  
  b<-function(x)
  {
    # b est la fonction de variance moyenne
    return(x*(1-x))
  }
  
  rapport<-function(x)
  {
    retour<-NULL
    for (y in x)
    {
      retour<-c(retour,a(y)/b(y))
    }
    return(retour)
  }
  
  phi<-function(y)
  {
    retour<-NULL
    for (x in y)
    {
      retour<-c(retour,exp(-2*integrate(rapport,0,y)$value))
    }
    return(retour)

  }
  
  temps<-NULL
  
  P0<-function(p)
  {
    return(integrate(phi,p,1)$value/integrate(phi,0,1)$value)
  }
  
  P1<-function(p)
  {
    return(integrate(phi,0,p)$value/integrate(phi,0,1)$value)
  }
  
  
  for (k in p)
  {

    t<-function(x)
    {
      retour<-NULL
      
      for (y in x)
      {
        if (y<k)
        {
          retour<-c(retour,(2*P0(k)/(b(y)*phi(y))*integrate(phi,0,y)$value))
        }
        else
        {
          retour<-c(retour,2*P1(k)/(b(y)*phi(y))*integrate(phi,y,1)$value)
        }
      }
    return(retour)
    }
    
     temps<-c(temps,integrate(t,0,1)$value)
  }
  return(temps)
  

}

testTempsAbso<-function(N=1000,x=c(0,1),s=0,h=0,u=0,v=0,n=100,pas=0.1)
{
  plot(0,0,type="n",xlab="p",ylab="temps moyen d'absorption",main="Test de l'approximation de diffusion",xlim=c(0,1),ylim=c(0,2))

  tempsTH<-WFcalculeTempsAbso(p=pas,N=N,s=s,h=h,u=u,v=v)
  simulation<-WFsimuleTempsAbso(p=pas,n=n,N=N,x=x,s=s,h=h,u=u,v=v,norm=TRUE)
  tempsEX<-mean(simulation)
  aux<-var(simulation)/(sqrt(n))
  ICup95<-tempsEX+2*aux
  ICdown95<-tempsEX-2*aux
  ICup99<-tempsEX+3*aux
  ICdown99<-tempsEX-3*aux
  
  for (p in seq(from=2*pas,to=1-pas,by=pas))
  {
    tempsTH1<-WFcalculeTempsAbso(p=p,N=N,s=s,h=h,u=u,v=v)
    simulation<-WFsimuleTempsAbso(p=p,n=n,N=N,x=x,s=s,h=h,u=u,v=v,norm=TRUE)
    tempsEX1<-mean(simulation)
    aux<-var(simulation)/(sqrt(n))
    ICup951<-tempsEX1+2*aux
    ICdown951<-tempsEX1-2*aux
    ICup991<-tempsEX1+3*aux
    ICdown991<-tempsEX1-3*aux
    
    segments(p-pas,tempsTH,p,tempsTH1,col="red")
    segments(p-pas,tempsEX,p,tempsEX1,col="blue")
    segments(p-pas,ICup95,p,ICup951,col="green")
    segments(p-pas,ICdown95,p,ICdown951,col="green")
    segments(p-pas,ICup99,p,ICup991,col="cyan")
    segments(p-pas,ICdown99,p,ICdown991,col="cyan")
    
    tempsTH<-tempsTH1
    tempsEX<-tempsEX1
    ICup95<-ICup951
    ICdown95<-ICdown951
    ICup99<-ICup991
    ICdown99<-ICdown991
  }
}

WFNeutreP<-function(p,t)
{
  ## Renvoie P1(p,t)
  ## Voir Ewens - Mathematical Population Genetics p162 (5.28)
  ## Pour obtenir P0, faire P1(1-p,t)
  
    C<-acos(2*p-1)^2/4
    print(C)
    
    return((p*(1-p)/C)^0.25*exp(-2*C/t))
}


WFechantillonne<-function(tk,nk,X0,N,n,s=0,h=0,u=0,v=0)
{
  ## Genere des trajectoires de Wright Fisher puis les échantillonnes 
  ## à des temps fixés à l'avance
  ## Le retour est une matrice de 4 colonnes
  ## La première colonne est la colonne des temps tk (en génération)
  ## La deuxieme colonne est la colonne des tailles des échantillonnages nk
  ## La troisièmecolonne est la colonne des fréquences empiriques
  ## La quatrième colonne est la colonne des fréquences réelles
  
  ## X0 est la fréquence allélique initiale
  ## N est tel que la population est de taille 2N
  ## n est le nombre de générations à réaliser
  ## s est le parametre de fitness
  ## h est le parametre de heterozygote fitness
  ## u est le taux de mutation de l'allele 1 vers 2
  ## v est le taux de mutation de l'allele 2 vers 1
  
  ## Voir Ewens - Mathematical Population Genetics p26 (1.67)
  
  ## tk est le vecteur contenant les temps (en générations) à échantillonner
  ## nk est la taille de l'échantillon à la date tk
  
  frequenceReelle<-WFtrajectoires(X0,N,n,s=0,h=0,u=0,v=0)
  resultat<-matrix(0,nrow=length(tk),ncol=4)
  
  echantillonne<-function(pas){return(rbinom(1,size=nk[pas],prob=frequenceReelle[tk[pas]]))}
  
  
  
  
  resultat[,1]<-tk
  resultat[,2]<-nk
  resultat[,3]<-unlist(lapply(1:length(nk),FUN=echantillonne))
  resultat[,4]<-frequenceReelle[tk]
  
  return(resultat)
}

#write.table(WFechantillonne(tk=seq(from=1,to=10000,by=100),nk=rep(100,length(seq(from=1,to=10000,by=100))),X0=0.5,N=10000,n=10000),file="fichierTest.txt")

WFplotTrajectoireEchIncertitude<-function(x0,N,n,s,h,u,v,tk,nk)
{
  ## affiche une trajectoire de Wright Fisher partant de x0 avec les parametres N,n,s,h,u,v
  ## affiche les valeurs des échantillonages aux temps tk avec des tailles d'échantillon nk
  ## affiche les intervalle de confiance de ces echantillonages
  frequenceReelle<-WFtrajectoires(x0,N,n,s,h,u,v)
  
  ech<- matrix(0,nrow=length(tk),ncol=3)
  
  echantillonne<-function(pas){return(rbinom(1,size=nk[pas],prob=frequenceReelle[tk[pas]]))}
    
  obs<-unlist(lapply(1:length(nk),FUN=echantillonne))
  
  plot(frequenceReelle,
       main = 'Allele frequency evolution',
       xlab = 'time in generation',
       ylab = 'allele frequency',
       col = 'red',
       type='l',
       ylim = c(0,1))
  
  mean = frequenceReelle[tk]
  
  
  variance = mean*(1-mean)
  
  
  upperlimit = mean + 1.96*sqrt(variance)/sqrt(nk)
  lowerlimit = mean - 1.96*sqrt(variance)/sqrt(nk)
  
  
  points(tk,
         obs/nk,
         col='blue',
         pch=4)
  
  
  

  #install.packages("plotrix")
  require(plotrix)
  plotCI(tk,
         y=mean,
         uiw=upperlimit-mean,
         liw=mean-lowerlimit,
         err="y",
         col='darkgreen',
         add=TRUE,
         pch=NA)
  
  legend("bottomright",
         legend=c("Real allele frequency","sample","variability of the sampling"),
         col=c('red','blue','darkgreen'),
         lty=c(1,0,1),
         pch=c(NA,4,NA)
         )
  
  

  
  
  
  
}

#tk = seq(1,1000,200)
## #nk = c(10,12,20,30,100)

## nombrePoints = 5

## tk =(1000/nombrePoints)*1:nombrePoints
## nk = rep(100,nombrePoints)


## WFplotTrajectoireEchIncertitude(0.1,1000,1000,0.01,0.5,0,0,tk,nk)
#dev.print(device = pdf , file = "IMAGEAFE")

