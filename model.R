#!usr/bin/env Rscript

library(gtools)
library(stringi)
library(plyr)
library(fdrtool)
library(MASS)
library(BiasedUrn)

runmodel<-function(nohostspecies=1, pinf=0.3, pclear=0.05, reinfect=0, distribution="uniform", corr=0, meet=0.5) 
{
  agentID <- 1
  size <- 2000
  n <- 1000
  ##defined globally - why?
  N <<- 5000
  initialh1 <- n/2
  initialinfections <- 10
  initialinfectionsh1 <- initialinfections/2
  its <- 200
  repeats <- 1
  seqlength <- 10
  infectionsethost1 <- matrix(NA,nrow=(its+1),ncol=repeats) ##initialise data storage
  infectionsethost2 <- matrix(NA,nrow=(its+1),ncol=repeats) ##initialise data storage
  fitnesses1 <- array(NA,dim=c(n,(its+1),repeats))
  fitnesses2 <- array(NA,dim=c(n,(its+1),repeats))
  variantshost1 <- matrix(NA,nrow=(its+1),ncol=repeats) ##initialise data storage
  variantshost2 <- matrix(NA,nrow=(its+1),ncol=repeats) ##initialise data storage
  variants <- matrix(NA,nrow=(its+1),ncol=repeats) ##initialise data storage
  for (r in 1:repeats) {
    tlist <- worldgen(n,size,agentID,seqlength,nohostspecies,initialh1,initialinfections,initialinfectionsh1,corr,distribution) ##generate initial model state
    model <- tlist[[1]]
    agentID <- tlist[[2]]
    landscape <- tlist[[3]]
    rm(tlist)
    s <- 0
    tlist <- recording(s,r,its,model,landscape,nohostspecies,infectionsethost1,fitnesses1,variantshost1,infectionsethost2,fitnesses2,variantshost2,variants)
    infectionsethost1 <- tlist[[1]]
    infectionsethost2 <- tlist[[2]]
    fitnesses1 <- tlist[[3]]
    fitnesses2 <- tlist[[4]]
    variantshost1 <- tlist[[5]]
    variantshost2 <- tlist[[6]]
    variants <- tlist[[7]]
    rm(tlist)
    rm(s)
    for (s in (1:its)) { ##loop through model interations
      model<-lapply(seq_along(model), function (x) { ##unassign all agents
        model[[x]]$positioned <- 0
        return(model[[x]])
      })
      model<-lapply(seq_along(model), function (x) { ##set all agents as old
        model[[x]]$new <- 0
        return(model[[x]])
      })
      unassigned<-length(model) ##set number of unassigned agents
      space<-size
      for (q in seq_along(model)) { ##loop though agents
        if (unassigned==0) { ##if all agents assigned...
          break() ##move to next step of model
        }
        if (model[[q]]$positioned==1) {
          next()
        }
        model[[q]]$positioned <- 1 ##mark focal agent positioned
        unassigned <- unassigned-1 ##mark focal agent assigned
        tlist <- hostconditioning(model[[q]]$ID,model,agentID,landscape,pclear,reinfect) ##test for host survival and viral clearance
        model <- tlist[[1]] ##reassign model
        agentID <- tlist[[3]] ##update agentID
        if (!tlist[[2]]) { ##if the host survives (now always returns false as no host death (functional metapopulation))
          rm(tlist)
          numsharers <- rbinom(1,unassigned,1/space)
          if (numsharers!=0) {
          	unassigned <- unassigned-numsharers ##assign selected agents
          	if (nohostspecies!=1 & meet!=0.5) { ##if biased contact ratios
          		sharersh1 <- lapply(seq_along(model), function (x) { ##subset out assignable agents
              		if (model[[x]]$positioned==0 & model[[x]]$hostspecies==1) {
                	return(model[[x]])
              		}
            	})
            	sharersh1<-Filter(Negate(function (x) is.null(unlist(x))), sharersh1) ##filter out NULLs
            	sharersh2<-lapply(seq_along(model),function (x) { ##subset out assignable agents
              		if (model[[x]]$positioned==0&model[[x]]$hostspecies==2) {
                	return(model[[x]])
              		}
            	})
            	sharersh2<-Filter(Negate(function (x) is.null(unlist(x))), sharersh2) ##filter out NULLs
            	if (model[[q]]$hostspecies==1) { ##use wullenius' noncentral hypergeometric to generate host distribution over sharers given differential contact probabilities
            		numh1<-rWNCHypergeo(1,length(sharersh1),length(sharersh2), numsharers,(meet/(1-meet)))
            	} else {
            		numh1<-rWNCHypergeo(1,length(sharersh1),length(sharersh2), numsharers,((1-meet)/meet))
            		}
            	if (!(numh1==0|numh1==numsharers)) { ## draw sharers following host distribution
            		sharersh1<-sample(sharersh1,numh1,replace=F)
            		sharersh2<-sample(sharersh2,(numsharers-numh1),replace=F)
            		sharers<-c(sharersh1,sharersh2)
            	} else if (numh1==0) {
            		sharers<-sample(sharersh2,(numsharers-numh1),replace=F)
            	} else {
            		sharers<-sample(sharersh1,numh1,replace=F)
            	}
            	sharers<-sample(sharers,length(sharers),replace=F) ##randomise
          	} else { ##if not biased contact ratios
            	sharers<-lapply(seq_along(model),function (x) { ##subset out assignable agents
              		if (model[[x]]$positioned==0) {
                	return(model[[x]])
              		}
            	})
            	sharers<-Filter(Negate(function (x) is.null(unlist(x))), sharers)
            	sharers<-sample(sharers,numsharers,replace=F)
            }
            for (y in 1:length(sharers)) {
              tlist<-hostconditioning(sharers[[y]]$ID,model,agentID,landscape,pclear,reinfect)
              model<-tlist[[1]] ##reassign model
              agentID<-tlist[[3]] ##update agentID
              sharers[[y]]$dead<-tlist[[2]]
              rm(tlist)
            }
            sharers<-lapply(seq_along(sharers), function (x) {
              if (!sharers[[x]]$dead) {return(sharers[[x]])}
            })
            sharers<-Filter(Negate(function (x) is.null(unlist(x))), sharers)
            sharers<-c(sharers,list(model[[q]]))
            test<-do.call("rbind",lapply(seq_along(sharers), function (x) {
              if (!is.na(sharers[[x]]$virusgenotype)) {
                return(1)
              } else {
                return(0)
              }
            }))
            if (sum(test)>0) {
              model<-infection(sharers,model,landscape,pinf)
            }
            rm(test)
          }
        } else {
          rm(tlist)
        }
      }
      space<-space-1 
      model<-sample(model,replace=F)
      model<-evolve(model,landscape)
      tlist<-recording(s,r,its,model,landscape,nohostspecies,infectionsethost1,fitnesses1,variantshost1,infectionsethost2,fitnesses2,variantshost2,variants)
      infectionsethost1<-tlist[[1]]
      infectionsethost2<-tlist[[2]]
      fitnesses1<-tlist[[3]]
      fitnesses2<-tlist[[4]]
      variantshost1<-tlist[[5]]
      variantshost2<-tlist[[6]]
      variants<-tlist[[7]]
      rm(tlist)
      print(paste("Repeat ",paste(r, paste(", Iteration ",s ,sep=""), sep=""), sep=""))
    }
  }
  return(list(infectionsethost1,infectionsethost2,fitnesses1,fitnesses2,variantshost1,variantshost2,variants,model,landscape,list(nohostspecies,pinf,pclear,reinfect,distribution,corr,meet)))
}

infection<-function(sharers,model,landscape,pinf)
{
  genotypes<-do.call("rbind",lapply(seq_along(sharers), function (x) {
    return(sharers[[x]]$virusgenotype)
  }))
  genotypes<-unique(genotypes[!is.na(genotypes)])
  if (length(genotypes)==0) {
    return(model)
  }
  infectionmatrix<-matrix(NA,nrow=length(sharers),ncol=length(genotypes))
  for (i in seq_along(genotypes)) {
    infectionmatrix[,i]<-do.call("rbind",lapply(seq_along(sharers), function (x) {
      if (runif(1,0,1)<=infectionprob(sharers[[x]],landscape,genotypes[i],pinf) && !genotypes[i]%in%sharers[[x]]$immunegenotypes) {
        return(1)
      } else {
        return(0)
      }
    }))
  }
  for (i in 1:nrow(infectionmatrix)) {
    if (max(infectionmatrix[i,]==0)) {
      next()
    } else {
      withinhostfitnesses<-do.call("rbind",lapply(genotypes[which(infectionmatrix[i,]==max(infectionmatrix[i,]))], function (x) {
        return(landscape[which(landscape[,1]==x),(sharers[[i]]$hostspecies+1)])
      }))
      tempgenos<-genotypes[which(infectionmatrix[i,]==max(infectionmatrix[i,]))]
      if (!is.na(sharers[[i]]$virusgenotype)) {
        if (max(withinhostfitnesses)>landscape[which(landscape[,1]==sharers[[i]]$virusgenotype),(sharers[[i]]$hostspecies+1)]) {
          sharers[[i]]$virusgenotype<-tempgenos[which.max(withinhostfitnesses)]
        }
      } else {
        sharers[[i]]$virusgenotype<-tempgenos[which.max(withinhostfitnesses)]
      }
      rm(tempgenos)
    }
  }
  for (i in seq_along(sharers)) {
    target<-which.max(do.call("rbind",lapply(seq_along(model), function (x) {
      if (model[[x]]$ID==sharers[[i]]$ID) {
        return(1)
      } else {
        return(0)
      }
    })))
    model[[target]]<-sharers[[i]]
  }
  return(model)
}       

hostconditioning<-function(hostID,model,agentID,landscape,pclear,reinfect) 
{
  target<-which.max(do.call("rbind",lapply(seq_along(model), function (x) {
    if (model[[x]]$ID==hostID) {
      return(1)
    } else {
      return(0)
    }
  })))
  model[[target]]$positioned<-1
  if (!is.na(model[[target]]$virusgenotype)) {
      if (runif(1,0,1)<=clearanceprob(model[[target]],landscape,pclear)) { ##if random number determines agent clearance...
      	if (reinfect==0) {
        	if (anyNA(model[[target]]$immunegenotypes)) { ##add cleared genotype to immune list
          	model[[target]]$immunegenotypes<-model[[target]]$virusgenotype
        	} else {
          	model[[target]]$immunegenotypes<-append(model[[target]]$virusgenotype,model[[target]]$immunegenotypes)
        	}
        }
        model[[target]]$virusgenotype<-NA ##clear virus
      }
    }
  dead<-FALSE
  return(list(model,dead,agentID))
}

infectionprob<-function(host,landscape,virusgenotype,pinf)
{
  if (landscape[which(landscape[,1]==virusgenotype),(host$hostspecies+1)]==0) {
    return(0)
  } else {
    return(pinf)
  }
}

clearanceprob<-function(host,landscape,pclear)
{
  if (landscape[which(landscape[,1]==host$virusgenotype),(host$hostspecies+1)]==0) {
    return(1)
  } else {
    return(pclear)
  }
}

worldgen<-function(n,size,agentID,sequlength,nohostspecies,initialh1,initialinfections,initialinfectionsh1,corr,distribution)
{
  model<-vector("list",length=n) #generate agents
  model<-lapply(seq_along(model),function(x) {
    newhost<-vector("list") ##generate each host
    if (nohostspecies==2) {
    	if (x>initialh1) { ##currently only defined for 2 host species
      		newhost[["hostspecies"]]<-2
    	} else {
      	newhost[["hostspecies"]]<-1
    	}
    } else {
    	newhost[["hostspecies"]]<-1
    }
    newhost[["condition"]]<-runif(1,0,1) ##condition currently not used, can be removed
    newhost[["virusgenotype"]]<-NA ##start uninfected
    newhost[["ID"]]<-x ##agent ID of new host
    newhost[["positioned"]]<-0
    newhost[["status"]]<-0
    newhost[["new"]]<-0
    newhost[["immunegenotypes"]]<-NA
    newhost[["dead"]]<-FALSE
    model[[x]]<-newhost
  })
  sequences<-seqgen(sequlength) ##generate all sequences of seqlength
  landscape<-matrix(NA,ncol=nohostspecies+1,nrow=length(sequences)) ##define fitness landscape
  landscape<-as.data.frame(landscape)
  landscape[,1]<-sequences
  if (nohostspecies==1) {
  	if (distribution=="uniform") {
  		landscape[,2]<-runif(nrow(landscape),0,1) ##uniform uncorrelated
  	} else if (distribution=="halfnormal") {
  		landscape[,2]<-rhalfnorm(nrow(landscape)) ##halfnormal uncorrelated
  	} else {
  		landscape[,2]<-rgamma(nrow(landscape),1) ##gamma uncorrelated
  	}
  	for (i in 1:initialinfectionsh1) {
    	model[[i]]$virusgenotype<-sample(landscape[,1],1) ##provide initial infections
  	}
  } else {
  	sigma<-matrix(c(1,corr,
  					corr,1),nrow=2,ncol=2) ##correlation matrix for copula
  	if (distribution=="uniform") {
  		landscape[,c(2,3)]<-mvrnorm(nrow(landscape),c(0,0),sigma) ##generate MV normal variates
  		landscape[,2]<-pnorm(landscape[,2]) ##convert to correlated uniform
  		landscape[,3]<-pnorm(landscape[,3]) ##convert to correlated uniform
  	} else if (distribution=="halfnormal") {
  		landscape[,c(2,3)]<-mvrnorm(nrow(landscape),c(0,0),sigma) ##generate MV normal variates
  		landscape[,2]<-pnorm(landscape[,2]) ##convert to correlated uniform
  		landscape[,3]<-pnorm(landscape[,3]) ##convert to correlated uniform
  		landscape[,2]<-qhalfnorm(landscape[,2]) ##convert to from correlated uniform to correlated halfnormal
  		landscape[,3]<-qhalfnorm(landscape[,3]) ##convert to from correlated uniform to correlated halfnormal
  	} else {
  		landscape[,c(2,3)]<-mvrnorm(nrow(landscape),c(0,0),sigma) ##generate MV normal variates
  		landscape[,2]<-pnorm(landscape[,2]) ##convert to correlated uniform
  		landscape[,3]<-pnorm(landscape[,3]) ##convert to correlated uniform
  		landscape[,2]<-qgamma(landscape[,2],1) ##convert to from correlated uniform to correlated gamma
  		landscape[,3]<-qgamma(landscape[,3],1) ##convert to from correlated uniform to correlated gamma
  	}
  	for (i in 1:initialinfectionsh1) {
    	model[[i]]$virusgenotype<-sample(landscape[,1],1) ##provide initial infections
  	}
  	for (i in 1:initialinfectionsh1) {
    	model[[i]]$virusgenotype<-sample(landscape[,1],1)
  	}
  	for (i in (initialh1+1):((initialh1+1)+(initialinfections-initialinfectionsh1))) {
    	model[[i]]$virusgenotype<-sample(landscape[,1],1)
  	}
  }
  model<-sample(model,replace=F)
  agentID<-n+1
  return(list(model,agentID,landscape))
}

seqgen<-function(seqlength)
{
  bases<-c('A','T','C','G')
  seqs<-permutations(4,seqlength,bases,repeats.allowed = T)
  final<-NA
  for (i in 1:(ncol(seqs)-1)) {
    if (i==1) {
      final<-paste(seqs[,i],seqs[,i+1],sep="")
    } else {
      final<-paste(final,seqs[,i+1],sep="")
    }
    next()
  }
  return(final)
}

onestep<-function(inputseq)
{
  bases<-c('A','T','C','G')
  k<-sample(c(1:nchar(inputseq)),1)
  bases<-bases[!bases%in%substr(inputseq,k,k)]
  stri_sub(inputseq,k,k)<-bases[sample(1:3,1)]
  return(inputseq)
}

newagent<-function(agentID,hostspecies) 
{
  newhost<-vector("list") ##generate each host
  newhost[["hostspecies"]]<-hostspecies ##give host species of previous agent
  newhost[["condition"]]<-runif(1,0,1) ##condition currently not used, can be removed
  newhost[["virusgenotype"]]<-NA ##start uninfected
  newhost[["ID"]]<-agentID ##agent ID of new host
  newhost[["positioned"]]<-1
  newhost[["status"]]<-0
  newhost[["new"]]<-1
  newhost[["immunegenotypes"]]<-NA
  newhost[["dead"]]<-FALSE
  agentID<-agentID+1
  return(list(newhost,agentID))
}

evolve<-function(model,landscape)
{
  evolved<-lapply(seq_along(model), function (x) {
    if (is.na(model[[x]]$virusgenotype)) {
      return(model[[x]])
    } else {
      newseq<-sample(onestep(model[[x]]$virusgenotype),1)
      if (landscape[which(landscape[,1]==newseq),(model[[x]]$hostspecies+1)]==landscape[which(landscape[,1]==model[[x]]$virusgenotype),(model[[x]]$hostspecies+1)]) {
      	if ((1/N)>=runif(1,0,1)) {
      		model[[x]]$virusgenotype<-newseq
      	}
      } else {
      	selectco<-(landscape[which(landscape[,1]==newseq),(model[[x]]$hostspecies+1)]/landscape[which(landscape[,1]==model[[x]]$virusgenotype),(model[[x]]$hostspecies+1)])-1
      	if ((1-exp(-selectco))/(1-exp(-selectco*N))>=runif(1,0,1)) {
      		model[[x]]$virusgenotype<-newseq
      	}
      }
      return(model[[x]])
    }
  })
  return(evolved)
}

recording<-function(s,r,its,model,landscape,nohostspecies,infectionsethost1,fitnesses1,variantshost1,infectionsethost2=NULL,fitnesses2=NULL,variantshost2=NULL,variants=NULL)
{
  if (nohostspecies==1) {
    infectionsethost1[(s+1),r]<-sum(do.call("rbind",lapply(seq_along(model), function (x) {
      if (!is.na(model[[x]]$virusgenotype) && model[[x]]$hostspecies==1) {
        return(1)
      } else {
        return(0)
      }
    })))
    fitnesseshost1<-do.call("rbind",lapply(seq_along(model), function (x) {
      if (!is.na(model[[x]]$virusgenotype) && model[[x]]$hostspecies==1) {
        return(landscape[which(landscape[,1]==model[[x]]$virusgenotype),model[[x]]$hostspecies+1])
      } else {
        return(NA)
      }
    }))
    fitnesses1[c(1:length(fitnesseshost1)),(s+1),r]<-fitnesseshost1
    variants1<-do.call("rbind",lapply(seq_along(model), function (x) {
      if (!is.na(model[[x]]$virusgenotype) && model[[x]]$hostspecies==1) {
        return(model[[x]]$virusgenotype)
      } else {
        return(NA)
      }
    }))
    variantshost1[(s+1),r]<-length(unique(variants1[!is.na(variants1)]))
    return(list(infectionsethost1,infectionsethost2,fitnesses1,fitnesses2,variantshost1,variantshost2,variants))
  } else {
    infectionsethost1[(s+1),r]<-sum(do.call("rbind",lapply(seq_along(model), function (x) {
      if (!is.na(model[[x]]$virusgenotype) && model[[x]]$hostspecies==1) {
        return(1)
      } else {
        return(0)
      }
    })))
    infectionsethost2[(s+1),r]<-sum(do.call("rbind",lapply(seq_along(model), function (x) {
      if (!is.na(model[[x]]$virusgenotype) && model[[x]]$hostspecies==2) {
        return(1)
      } else {
        return(0)
      }
    })))
    fitnesseshost1<-do.call("rbind",lapply(seq_along(model), function (x) {
      if (!is.na(model[[x]]$virusgenotype) && model[[x]]$hostspecies==1) {
        return(landscape[which(landscape[,1]==model[[x]]$virusgenotype),model[[x]]$hostspecies+1])
      } else {
        return(NA)
      }
    }))
    fitnesses1[c(1:length(fitnesseshost1)),(s+1),r]<-fitnesseshost1
    variants1<-do.call("rbind",lapply(seq_along(model), function (x) {
      if (!is.na(model[[x]]$virusgenotype) && model[[x]]$hostspecies==1) {
        return(model[[x]]$virusgenotype)
      } else {
        return(NA)
      }
    }))
    variantshost1[(s+1),r]<-length(unique(variants1[!is.na(variants1)]))
    fitnesseshost2<-do.call("rbind",lapply(seq_along(model), function (x) {
      if (!is.na(model[[x]]$virusgenotype) && model[[x]]$hostspecies==2) {
        return(landscape[which(landscape[,1]==model[[x]]$virusgenotype),model[[x]]$hostspecies+1])
      } else {
        return(NA)
      }
    }))
    fitnesses2[c(1:length(fitnesseshost2)),(s+1),r]<-fitnesseshost2
    variants2<-do.call("rbind",lapply(seq_along(model), function (x) {
      if (!is.na(model[[x]]$virusgenotype) && model[[x]]$hostspecies==2) {
        return(model[[x]]$virusgenotype)
      } else {
        return(NA)
      }
    }))
    variantshost2[(s+1),r]<-length(unique(variants2[!is.na(variants2)]))
    variantsb<-do.call("rbind",lapply(seq_along(model), function (x) {
      if (!is.na(model[[x]]$virusgenotype)) {
        return(model[[x]]$virusgenotype)
      } else {
        return(NA)
      }
    }))
    variants[(s+1),r]<-length(unique(variantsb[!is.na(variantsb)]))
    return(list(infectionsethost1,infectionsethost2,fitnesses1,fitnesses2,variantshost1,variantshost2,variants))
  }
}



data<-runmodel()
save(data,file=paste("onehostbase",paste(Sys.time(),"Rdata",sep="."),sep="."))