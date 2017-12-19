#!usr/bin/env Rscript
library(fdrtool)
files<-list.files(pattern="*.Rdata")
results<-matrix(NA,ncol=22,nrow=length(files))
results<-as.data.frame(results)
colnames(results)<-c("Hosts","Infection Probability","Clearance Probability","Reinfection","Landscape Distribution","Correlation","Contact Bias","Time to Max Fitness Host 1","Time to Max Fitness Host 2","Number of Variants Host 1","Number of Variants Host 2","Number of Variants Total","Variance in Fitness Host 1","Variance in Fitness Host 2","Variance in Fitness Total","Number Infected Host 1","Number Infected Host 2","Correlation","Proportion exceeding the median in both host types","Max fitness host 1","Max fitness host 2","Time to 0.999 H1")
for (i in 1:length(files)) {
	load(files[i])
	results[i,1]<-data[[10]][[1]]
	results[i,2]<-data[[10]][[2]]
	results[i,3]<-data[[10]][[3]]
	results[i,4]<-data[[10]][[4]]
	results[i,5]<-data[[10]][[5]]
	results[i,6]<-data[[10]][[6]]
	results[i,7]<-data[[10]][[7]]
	if (sum(!is.na(data[[3]][,,1]))>0) {
		target<-max(data[[3]][,,1],na.rm=T)
		results[i,20]<-target
		if (data[[10]][[5]]=="uniform") {
			target2<-qunif(0.999)
			if (target2>target) {
				results[i,22]<-200
			}
			hit<-FALSE
			for (q in 1:dim(data[[3]])[2]) {
				if ((sum(data[[3]][,q,]>=target2,na.rm=T)>0)&(hit==FALSE)) {
					hit<-TRUE
					results[i,22]<-q
				}
				if (sum(data[[3]][,q,]%in%target)>0) {
					results[i,8]<-q
					break()
				}
			}
		} else if (data[[10]][[5]]=="gamma") {
			target2<-qgamma(0.999,1)
			if (target2>target) {
				results[i,22]<-200
			}
			hit<-FALSE
			for (q in 1:dim(data[[3]])[2]) {
				if ((sum(data[[3]][,q,]>=target2,na.rm=T)>0)&(hit==FALSE)) {
					hit<-TRUE
					results[i,22]<-q
				}
				if (sum(data[[3]][,q,]%in%target)>0) {
					results[i,8]<-q
					break()
				}
			}
		} else {
			target2<-qhalfnorm(0.999)
			if (target2>target) {
				results[i,22]<-200
			}
			hit<-FALSE
			for (q in 1:dim(data[[3]])[2]) {
				if ((sum(data[[3]][,q,]>=target2,na.rm=T)>0)&(hit==FALSE)) {
					hit<-TRUE
					results[i,22]<-q
				}
				if (sum(data[[3]][,q,]%in%target)>0) {
					results[i,8]<-q
					break()
				}
			}
		}
	}
	if (sum(!is.na(data[[4]][,,1]))>0) {
		target<-max(data[[4]][,,1],na.rm=T)
		results[i,21]<-target
		for (z in 1:dim(data[[4]])[2]) {
			if (sum(data[[4]][,z,1]%in%target)>0) {
				results[i,9]<-z
				break()
			}
		}
	}
	results[i,10]<-data[[5]][201,1]
	results[i,11]<-data[[6]][201,1]
	results[i,12]<-data[[7]][201,1]
	results[i,13]<-var(data[[3]][,201,1],na.rm=T)
	results[i,14]<-var(data[[4]][,201,1],na.rm=T)
	results[i,15]<-var(c(data[[3]][,201,1],data[[4]][,201,1]),na.rm=T)
	results[i,16]<-data[[1]][201,1]
	results[i,17]<-data[[2]][201,1]
	if (data[[10]][[1]]==2) {
	if ((data[[1]]+data[[2]])[201]!=0) {
		fh1<-unique(data[[3]][!(data[[3]][,201,1]%in%NA),201,1])
		fh2<-unique(data[[4]][!(data[[4]][,201,1]%in%NA),201,1])
		genotypes<-unique(c(data[[9]][(data[[9]][,2]%in%fh1),1],data[[9]][(data[[9]][,3]%in%fh2),1]))
		if (length(genotypes>1)) {
			results[i,18]<-cor(data[[9]][(data[[9]][,1]%in%genotypes),2],data[[9]][(data[[9]][,1]%in%genotypes),3],method="kendall")
		}
		if (data[[10]][[5]]=="uniform") {
			results[i,19]<-sum((data[[9]][(data[[9]][,1]%in%genotypes),2]>qunif(0.5))&(data[[9]][(data[[9]][,1]%in%genotypes),3]>qunif(0.5)))
		} else if (data[[10]][[5]]=="gamma") {
			results[i,19]<-sum((data[[9]][(data[[9]][,1]%in%genotypes),2]>qgamma(0.5,1))&(data[[9]][(data[[9]][,1]%in%genotypes),3]>qgamma(0.5,1)))
		} else {
			results[i,19]<-sum((data[[9]][(data[[9]][,1]%in%genotypes),2]>qhalfnorm(0.5))&(data[[9]][(data[[9]][,1]%in%genotypes),3]>qhalfnorm(0.5)))
		}
		}
	}
	print(results[i,22])
	print(i)
}

write.csv(results,"results.csv")