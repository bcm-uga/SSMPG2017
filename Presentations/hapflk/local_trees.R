library(ape)
library(phangorn)

######################## USER INPUT #########################
### Pre filled with files from the tutorial
## tree file : the whole genome tree as obtained
## from the hapflk analysis
tree_file='practical/hapflk/NEU_tree.txt'
## wg_dist_file : the Global Reynolds distances between population,
## as obtained from the hapflk analysis
wg_dist_file='practical/hapflk/NEU_reynolds.txt'
## Name of the outgroup in the previous file (if any)
outgroup='Soay'
## local reynolds distances for SNPs and Haplotype clusters
## To be computed using the local_reynolds.py script
loc_dist_file_snp='practical/hapflk/hapflk_snp_reynolds.txt'
loc_dist_file_hap='practical/hapflk/hapflk_hap_reynolds.txt'
## prefix for outputfiles
prefix='practical/hapflk/mstn'
##################### END USER INPUT #########################

## Now run : R CMD BATCH local_trees.R 

####################### UTILS ##################################
parse.path=function(v,noms,ntips) {
  ### recontruit un chemin a partir de segments
  my.segs=noms[v!=0]
  my.nodes=simplify2array(sapply(my.segs,strsplit,split='<->'))
  my.nodes=t(apply(my.nodes,1,as.integer))
  ## get leaves in path 
  my.tips=as.integer(my.nodes[!(my.nodes>ntips)])
  if (length(my.tips)!=2) {
    print("error")
    print(list(segments=my.tips,ntips=ntips))
    return(NA)
  }
  my.tips
}

get.tips=function(X) {
  K=dim(X)[1]
  ntips=0.5*(1+sqrt(1+8*K))
  if (ntips != floor(ntips)) {
    print("stop the BS")
    return(NA)
  }
  branches=colnames(X)
  segments=t(apply(X,1,parse.path,noms=branches,ntips=ntips))
  segments
}

build.y=function(X,dm) {
  Nd=dim(X)[1]
  y.idx=get.tips(X)
  y=vector(length=Nd)
  for (i in 1:Nd) {
    y[i]=dm[y.idx[i,1],y.idx[i,2]]
  }
  y
}

####################### UTILS ##################################


#### Tree
mytree=unroot(read.tree(tree_file))
mytree=phangorn:::reorderPruning(mytree)

### Global tree plot
png(file=paste(prefix,'_WholeGenomeTree.png',sep=''),width=3.4,height=3.4
    ,unit='in',res=300)
par(mar=c(1,1,2,1),lwd=2)
plot(mytree,type='clado',label.offset=0.02,edge.width=2,cex=0.7)
tiplabels(bg='white',frame='circle',cex=0.7)
nodelabels(bg='white',frame='circle',cex=0.7)
title(main='Whole Genome tree')
dev.off()

#### Whole genome distances
wgd=read.table(wg_dist_file,row.names=1)
colnames(wgd)=rownames(wgd)
## drop outgroup
wgd=wgd[rownames(wgd)!=outgroup,colnames(wgd)!=outgroup]
## align wgd with tip labels
order.pop=match(mytree$tip.label,rownames(wgd))
wgd=wgd[order.pop,order.pop]

#### linear model of branch length
Xtree=designTree(mytree)
wgy=build.y(Xtree,wgd)
## you can check fit of the tree to the dist matrix
## summary(lm(wgy~Xtree-1))

### local distances
get_local_tree=function(loc_dist_file,prefix) {
    locd=read.table(loc_dist_file,row.names=1)
    ## drop outgroup
    locd=locd[rownames(locd)!=outgroup,colnames(locd)!=outgroup]

    order.pop=match(mytree$tip.label,rownames(locd))
    locd=locd[order.pop,order.pop]
    
    locy=build.y(Xtree,locd)
    locmod=lm(locy~Xtree-1,offset=wgy)
    mycoef=summary(locmod)$coefficients
    mycoef[,1]=round(mycoef[,1],digits=3)
    mycoef[,2]=round(mycoef[,2],digits=5)
    mycoef[,3]=round(mycoef[,3],digits=2)
    mycoef[,4]=format(mycoef[,4],scientific=TRUE,digits=2)
    write.table(mycoef,quote=FALSE,sep=',',file=paste(prefix,'_lmtab.csv',sep=''))

    loctree=mytree
    loctree$edge.length=mytree$edge.length+locmod$coefficients
    loctree$edge.length[loctree$edge.length<0]=0
    loctree$pvalue=summary(locmod)$coefficients[,4]
    loctree$coul=floor(-log10(loctree$pvalue))+1
    loctree$coul[loctree$coul>6]=6
    loctree$coul[locmod$coefficients<0]=1
    return(loctree)
}
snp_local_tree=get_local_tree(loc_dist_file_snp,prefix=paste(prefix,'_snp',sep=''))
hap_local_tree=get_local_tree(loc_dist_file_hap,prefix=paste(prefix,'_hap',sep=''))

##summary(locmod)

mycoul.snp=c("#EFF3FF","#C6DBEF","#9ECAE1","#6BAED6","#3182BD","#08519C")
mycoul.hap=c("#FEE5D9","#FCBBA1","#FC9272","#FB6A4A","#DE2D26","#A50F15")
##brewer.pal(6,'Blues')
##brewer.pal(6,'Reds')

### Local tree plot
png(file=paste(prefix,'_LocalSnpTree.png',sep=''),width=3.4,height=3.4
    ,unit='in',res=300)
par(mar=c(1,1,2,1),lwd=2)
plot(snp_local_tree,type='clado',label.offset=0.02,cex=0.7,
     edge.width=2,edge.color=mycoul.snp[snp_local_tree$coul])
title(main='SNP Local tree')
tiplabels(bg='white',frame='circle',cex=0.7)
nodelabels(bg='white',frame='circle',cex=0.7)
dev.off()

png(file=paste(prefix,'_LocalHapTree.png',sep=''),width=3.4,height=3.4
    ,unit='in',res=300)
par(mar=c(1,1,2,1),lwd=2)
plot(hap_local_tree,type='clado',label.offset=0.02,cex=0.7,
     edge.width=2,edge.color=mycoul.hap[hap_local_tree$coul])
title(main='Haplotypes Local tree')
tiplabels(bg='white',frame='circle',cex=0.7)
nodelabels(bg='white',frame='circle',cex=0.7)
dev.off()

## lm results

