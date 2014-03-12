#Calcul of CVNND for one trait with our without division by the range of the trait
CVNND<-function(trait, div_range=FALSE){
	
	r=sort(trait)
	if(length(r)<1){
	nnd=NA}
	
	else{nnd=vector(length=length(r)-1)
		for(j in 2:length(r)){
			nnd[j-1]=r[j]-r[j-1]
		}
	}
	
	CVNND<-sd(nnd,na.rm=T)/mean(nnd, na.rm=T)
	
	if (div_range==T) {CVNND<-CVNND/range(trait)} 
	else {}
	
	return(CVNND)
}

#Calcul of four statistics (mean, range, CVNND and kurtosis) to test community assembly using null models
#For each statistic this function calcul the mean, the standard deviation and two quantiles

#Proba can be use to take into account for abundance (c.f. com.i et al. 2010)
#if pop.aggregate=FALSE, ?????

com.index.c<-function(traits=NULL, ind.plot=NULL, sp=NULL, proba=NULL,nperm=99, val.quant=c(0.025,0.975), printprogress=TRUE, pop.aggregate=TRUE){
	
	name_sp_sites=paste(sp, ind.plot,sep="_")
	comm=NULL
	comm<-t(table(ind.plot,1:length(ind.plot)))
		
	if(is.null(proba) & pop.aggregate==FALSE)
		{proba=rep(1, times=dim(traits)[1])}  
	else {}
	
	if(is.null(proba) & pop.aggregate==TRUE) 
		{proba=rep(1, times=length(unique(name_sp_sites)))}  
	else {}
	
	S = colSums(comm>0)

	
	######################################### 
	#### 	  Calcul on null models  	 ####
	######################################### 
	if(pop.aggregate==T) {
			traits<-apply(traits, 2 , function (x) tapply(x, name_sp_sites, mean , na.rm=T))
			S<- table(unlist(strsplit(levels(as.factor(name_sp_sites)),split="_"))[3*(1:nlevels(as.factor(name_sp_sites)))])
	}
		
	if(is.numeric(nperm)){
		null_com=list()
		name=list(1:nperm,S)
		Null_range=list()
		Null_mean=list()
		Null_CV=list()
		Null_kurtosis=list()
		#__________
		Range_Mean<-matrix(nrow=dim(traits)[2], ncol=length(S))
		Range_sd<-matrix(nrow=dim(traits)[2], ncol=length(S))
		Range_qsup<-matrix(nrow=dim(traits)[2], ncol=length(S))
		Range_qinf<-matrix(nrow=dim(traits)[2], ncol=length(S))
		#__________
		Mean_Mean<-matrix(nrow=dim(traits)[2], ncol=length(S))
		Mean_sd<-matrix(nrow=dim(traits)[2], ncol=length(S))
		Mean_qsup<-matrix(nrow=dim(traits)[2], ncol=length(S))
		Mean_qinf<-matrix(nrow=dim(traits)[2], ncol=length(S))
		#__________
		CV.NND_Mean<-matrix(nrow=dim(traits)[2], ncol=length(S))
		CV.NND_sd<-matrix(nrow=dim(traits)[2], ncol=length(S))
		CV.NND_qsup<-matrix(nrow=dim(traits)[2], ncol=length(S))
		CV.NND_qinf<-matrix(nrow=dim(traits)[2], ncol=length(S))
		#__________
		kurtosis_Mean<-matrix(nrow=dim(traits)[2], ncol=length(S))
		kurtosis_sd<-matrix(nrow=dim(traits)[2], ncol=length(S))
		kurtosis_qsup<-matrix(nrow=dim(traits)[2], ncol=length(S))
		kurtosis_qinf<-matrix(nrow=dim(traits)[2], ncol=length(S))
		
		tik<-list()
			
		for (t in 1: dim(traits)[2]){
			null_com[[t]]<-list(length=length(S))
			name[[t]]<-list(1:nperm,S)
			Null_range[[t]]<-as.data.frame(matrix(nrow=length(name[[t]][[1]]),ncol=length(name[[t]][[2]]),dimnames=name[[t]]))
			Null_mean[[t]]<-as.data.frame(matrix(nrow=length(name[[t]][[1]]),ncol=length(name[[t]][[2]]),dimnames=name[[t]]))
			Null_CV[[t]]<-as.data.frame(matrix(nrow=length(name[[t]][[1]]),ncol=length(name[[t]][[2]]),dimnames=name[[t]]))
			Null_kurtosis[[t]]<-as.data.frame(matrix(nrow=length(name[[t]][[1]]),ncol=length(name[[t]][[2]]),dimnames=name[[t]]))
		   
			for(i in 1:nperm){
				for (s in 1:length(S)){
					tik[[t]]<-sample(traits[,t],S[s],prob=proba)
					null_com[[t]][[S[s]]]=tik[[t]]
					# Range
					Null_range[[t]][i,s]<-max(null_com[[t]][[S[s]]],na.rm=T)-min(null_com[[t]][[S[s]]],na.rm=T)
					# Mean
					Null_mean[[t]][i,s]<-mean(null_com[[t]][[S[s]]],na.rm=T)
					# kurtosis
					Null_kurtosis[[t]][i,s]<-kurtosis(null_com[[t]][[S[s]]],na.rm=T)
					# coefficient of variation of Nearest-Neighbor Distances (CV NND)
					Null_CV[[t]][i,s]<-CVNND(null_com[[t]][[S[s]]],)
				}
			}
			
			Range_Mean[t,]<-apply(Null_range[[t]],2,mean)
			Range_sd[t,]<-apply(Null_range[[t]],2,sd)
			Range_qsup[t,]<-apply(Null_range[[t]],2,function(x) quantile(x, prob=val.quant[2], na.rm=T))
			Range_qinf[t,]<-apply(Null_range[[t]],2,function(x) quantile(x, prob=val.quant[1], na.rm=T))
			#__________
			Mean_Mean[t,]<-apply(Null_mean[[t]],2,mean)
			Mean_sd[t,]<-apply(Null_mean[[t]],2,sd)
			Mean_qsup[t,]<-apply(Null_mean[[t]],2,function(x) quantile(x, prob=val.quant[2], na.rm=T))
			Mean_qinf[t,]<-apply(Null_mean[[t]],2,function(x) quantile(x, prob=val.quant[1], na.rm=T))
			#__________
			CV.NND_Mean[t,]<-apply(Null_CV[[t]],2,mean)
			CV.NND_sd[t,]<-apply(Null_CV[[t]],2,sd)
			CV.NND_qsup[t,]<-apply(Null_CV[[t]],2,function(x) quantile(x, prob=val.quant[2], na.rm=T))
			CV.NND_qinf[t,]<-apply(Null_CV[[t]],2,function(x) quantile(x, prob=val.quant[1], na.rm=T))
			#__________
			kurtosis_Mean[t,]<-apply(Null_kurtosis[[t]],2,mean)
			kurtosis_sd[t,]<-apply(Null_kurtosis[[t]],2,sd)
			kurtosis_qsup[t,]<-apply(Null_kurtosis[[t]],2,function(x) quantile(x, prob=val.quant[2], na.rm=T))
			kurtosis_qinf[t,]<-apply(Null_kurtosis[[t]],2,function(x) quantile(x, prob=val.quant[1], na.rm=T))
			
			if(printprogress==T){print(paste(round(t/dim(traits)[2]*100, 2), "%", sep=" "))} else {} 
		} 
	}
	
	else{}
	
	  
	######################################## 
	####	Calcul of observed values	####
	######################################## 
	
		
	range_obs<-matrix(nrow=dim(traits)[2] ,ncol=length(S))
	mean_obs<-matrix(nrow=dim(traits)[2] ,ncol=length(S))
	CV_obs<-matrix(nrow=dim(traits)[2] ,ncol=length(S))
	kurtosis_obs<-matrix(nrow=dim(traits)[2] ,ncol=length(S))
	
	if(pop.aggregate==T) 
		{interm<-lapply(strsplit(paste(rownames(traits), sep="_"), split="_"), function(x) x[3])}
	else if(pop.aggregate==F)
		{interm<-lapply(strsplit(paste(name_sp_sites, sep="_"), split="_"), function(x) x[3])}
	else {print("Error: pop.aggregate must be TRUE or FALSE")}
			
	# Range
	range_obs<-apply(traits, 2, function (x) tapply(x, unlist(interm), max, na.rm=T))-apply(traits, 2, function (x) tapply(x, unlist(interm), min, na.rm=T))

	# Mean
	mean_obs<-apply(traits, 2, function (x) tapply(x, unlist(interm), mean, na.rm=T))

	# kurtosis
	kurtosis_obs<-apply(traits, 2, function (x) tapply(x, unlist(interm), kurtosis, na.rm=T))

	# coefficient of variation of Nearest-Neighbor Distances (CV NND)
	CV_obs<-apply(traits, 2, function (x) tapply(x, unlist(interm), CVNND))
	
	
	######################################## 
	####		Create results list		####
	######################################## 
	
	com.i<-list()  
	
	com.i$Range_obs<-t(range_obs)
	com.i$Mean_obs<-t(mean_obs)
	com.i$CV.NND_obs<-t(CV_obs)
	com.i$kurtosis_obs<-t(kurtosis_obs)
	
	com.i$Range_obs[!is.finite(com.i$Range_obs)]<-NA
	com.i$Mean_obs[!is.finite(com.i$Mean_obs)]<-NA
	com.i$CV.NND_obs[!is.finite(com.i$CV.NND_obs)]<-NA
	com.i$kurtosis_obs[!is.finite(com.i$kurtosis_obs)]<-NA
	
	
	if(is.numeric(nperm)){
		com.i$Range_Mean<-Range_Mean
		com.i$Mean_Mean<-Mean_Mean
		com.i$CV.NND_Mean<-CV.NND_Mean
		com.i$kurtosis_Mean<-kurtosis_Mean
			
		com.i$Range_sd<-Range_sd
		com.i$Mean_sd<-Mean_sd
		com.i$CV.NND_sd<-CV.NND_sd
		com.i$kurtosis_sd<-kurtosis_sd
		
		com.i$Range_qsup<-Range_qsup
		com.i$Mean_qsup<-Mean_qsup
		com.i$CV.NND_qsup<-CV.NND_qsup
		com.i$kurtosis_qsup<-kurtosis_qsup
		
		com.i$Range_qinf<-Range_qinf
		com.i$Mean_qinf<-Mean_qinf
		com.i$CV.NND_qinf<-CV.NND_qinf
		com.i$kurtosis_qinf<-kurtosis_qinf
		
		com.i$Null_range<-Null_range
		com.i$Null_mean<-Null_mean
		com.i$Null_CV<-Null_CV
		com.i$Null_kurtosis<-Null_kurtosis
		
	}
	
	com.i$sites_richness<-S
	com.i$namestraits<-colnames(traits)
	
	return(com.i)
}


plot.com<-function(com.object=NULL, namestraits=NULL, sites_names=NULL, plot.ask=TRUE, ord=NULL){
	
	com.i<-com.object
	namestraits=com.i$namestraits		
	
	if(is.null(sites_names)) 
		{sites_names=colnames(com.i$Null_mean[[1]])} 
	else{}
	
	oldpar<-par(no.readonly = TRUE)
	
	par(ask=plot.ask)
	par(mfrow=c(2,2))
	par(mar=c(2,4,1,2))
	
	if (is.null(ord))
		{ord<-order(sites_names)}
	else{}	
	
	for(t in 1:length(namestraits)) {

		yli=c(min(((com.i$Range_obs[t,]-com.i$Range_Mean[t])/com.i$Range_sd[t])[ord], 0, na.rm=T)-3, max(((com.i$Range_obs[t,]-com.i$Range_Mean[t])/com.i$Range_sd[t])[ord], 0, na.rm=T)+3)
		plot(((com.i$Range_obs[t,]-com.i$Range_Mean[t])/com.i$Range_sd[t])[ord], bty="n", xaxt="n", pch=16, ylim=yli, main=namestraits[t] ,ylab="SES Range", xlab=NA)
		text(0.5+((com.i$Range_obs[t,]-com.i$Range_Mean[t])/com.i$Range_sd[t])[ord], labels=sites_names[ord])
		abline(h=0, lty=1)
		lines(((com.i$Range_qsup[t,] -com.i$Range_Mean[t])/com.i$Range_sd[t])[ord] , lty=2)
		lines(((com.i$Range_qinf[t,] -com.i$Range_Mean[t])/com.i$Range_sd[t])[ord] , lty=2)

		yli=c(min(((com.i$Mean_obs[t,]-com.i$Mean_Mean[t])/com.i$Mean_sd[t])[ord], 0, na.rm=T)-3, max(((com.i$Mean_obs[t,]-com.i$Mean_Mean[t])/com.i$Mean_sd[t])[ord], 0, na.rm=T)+3)
		plot(((com.i$Mean_obs[t,]-com.i$Mean_Mean[t])/com.i$Mean_sd[t])[ord], bty="n", xaxt="n", pch=16, main=namestraits[t], ylim=yli, ylab="SES Mean", xlab=NA)
		text(0.5+((com.i$Mean_obs[t,]-com.i$Mean_Mean[t])/com.i$Mean_sd[t])[ord], labels=sites_names[ord])
		abline(h=0, lty=1)
		lines(((com.i$Mean_qsup[t,] -com.i$Mean_Mean[t])/com.i$Mean_sd[t])[ord] , lty=2)
		lines(((com.i$Mean_qinf[t,] -com.i$Mean_Mean[t])/com.i$Mean_sd[t])[ord] , lty=2)

		yli=c(min(((com.i$CV.NND_obs[t,]-com.i$CV.NND_Mean[t])/com.i$CV.NND_sd[t])[ord], 0, na.rm=T)-3, max(((com.i$CV.NND_obs[t,]-com.i$CV.NND_Mean[t])/com.i$CV.NND_sd[t])[ord], 0, na.rm=T)+3)
		plot(((com.i$CV.NND_obs[t,]-com.i$CV.NND_Mean[t])/com.i$CV.NND_sd[t])[ord], bty="n", xaxt="n", pch=16, main=namestraits[t], ylim=yli, ylab="SES CV_NND", xlab=NA)
		text(0.5+((com.i$CV.NND_obs[t,]-com.i$CV.NND_Mean[t])/com.i$CV.NND_sd[t])[ord], labels=sites_names[ord])
		abline(h=0, lty=1)
		lines(((com.i$CV.NND_qsup[t,] -com.i$CV.NND_Mean[t])/com.i$CV.NND_sd[t])[ord] , lty=2)
		lines(((com.i$CV.NND_qinf[t,] -com.i$CV.NND_Mean[t])/com.i$CV.NND_sd[t])[ord] , lty=2)

		yli=c(min(((com.i$kurtosis_obs[t,]-com.i$kurtosis_Mean[t])/com.i$kurtosis_sd[t])[ord], 0, na.rm=T)-3, max(((com.i$kurtosis_obs[t,]-com.i$kurtosis_Mean[t])/com.i$kurtosis_sd[t])[ord], 0, na.rm=T)+3)
		plot(((com.i$kurtosis_obs[t,]-com.i$kurtosis_Mean[t])/com.i$kurtosis_sd[t])[ord], bty="n", xaxt="n", main=namestraits[t], pch=16, ylim=yli, ylab="SES kurtosis", xlab=NA)
		text(0.5+((com.i$kurtosis_obs[t,]-com.i$kurtosis_Mean[t])/com.i$kurtosis_sd[t])[ord],labels=sites_names[ord])
		abline(h=0, lty=1)
		lines(((com.i$kurtosis_qsup[t,] -com.i$kurtosis_Mean[t])/com.i$kurtosis_sd[t])[ord] , lty=2)
		lines(((com.i$kurtosis_qinf[t,] -com.i$kurtosis_Mean[t])/com.i$kurtosis_sd[t])[ord] , lty=2)
	}

par(oldpar)

}


mat.comm<-function(com.object=NULL, tstats=NULL,  val.quant=c(0.025,0.975), sites_names=NULL, method="p.val") {
	
	val.quanti<-val.quant
	
	if(!is.null(com.object) & !is.null(tstats))
		{nindex<-7 ; com.i<-com.object} 
	else if(is.null(com.object) & !is.null(tstats))
		{nindex<-3}
	else if(!is.null(com.object) & is.null(tstats))
		{nindex<-4 ; com.i<-com.object}
	else {stop("Both com.i.object and tstats are NULL")}
	
	if( nindex==7 | nindex==4 ) 
		{nplots<- length(com.i$sites_richness)}
	
	else{nplots<- dim(tstats$T_IP.IC)[1]}
	
	if(is.null(sites_names) & nindex==4)
		{sites_names<-as.factor(1:length(com.i$sites_richness))}
	 
	else if (is.null(sites_names))
		{sites_names<-as.factor(1:dim(tstats$T_IP.IC)[1])}	
		
	
	if(nindex==7){
		#________________________________________
		#Calcul of standardised effect size
		
		ses.T_IP.IC<-ses(tstats$T_IP.IC,tstats$T_IP.IC_nm, val.quant=val.quanti)
		ses.T_IC.IR<-ses(tstats$T_IC.IR,tstats$T_IC.IR_nm, val.quant=val.quanti)
		ses.T_PC.PR<-ses(tstats$T_PC.PR,tstats$T_PC.PR_nm, val.quant=val.quanti)
		
		ses.range<-ses(com.i$Range_obs, com.i$Null_range, val.quant=val.quanti)
		ses.mean<-ses(com.i$Mean_obs, com.i$Null_mean, val.quant=val.quanti)
		ses.CV<-ses(com.i$CV.NND_obs, com.i$Null_CV, val.quant=val.quanti)
		ses.kurtosis<-ses(com.i$kurtosis_obs, com.i$Null_kurtosis, val.quant=val.quanti)
	 	
		
		#________________________________________
		#Condition to be significantly different from null models with respect to values of quantile choose
		cond.T_IP.IC.inf<-ses.T_IP.IC$ses<ses.T_IP.IC$ses.inf
		cond.T_IC.IR.inf<-ses.T_IC.IR$ses<ses.T_IC.IR$ses.inf
		cond.T_PC.PR.inf<-ses.T_PC.PR$ses<ses.T_PC.PR$ses.inf
			
		cond.T_IP.IC.sup<-ses.T_IP.IC$ses>ses.T_IP.IC$ses.sup
		cond.T_IC.IR.sup<-ses.T_IC.IR$ses>ses.T_IC.IR$ses.sup
		cond.T_PC.PR.sup<-ses.T_PC.PR$ses>ses.T_PC.PR$ses.sup
		
		cond.range.inf<-t(ses.range$ses<ses.range$ses.inf)
		cond.mean.inf<-t(ses.mean$ses<ses.mean$ses.inf)
		cond.CV.NND.inf<-t(ses.CV$ses<ses.CV$ses.inf)
		cond.kurtosis.inf<-t(ses.kurtosis$ses<ses.kurtosis$ses.inf)
		
		cond.range.sup<-t(ses.range$ses>ses.range$ses.sup)
		cond.mean.sup<-t(ses.mean$ses>ses.mean$ses.sup)
		cond.CV.NND.sup<-t(ses.CV$ses>ses.CV$ses.sup)
		cond.kurtosis.sup<-t(ses.kurtosis$ses>ses.kurtosis$ses.sup)
				
		
	}
	else if(nindex==3){
		#________________________________________
		#Calcul of standardised effect size
		
		ses.T_IP.IC<-ses(tstats$T_IP.IC,tstats$T_IP.IC_nm, val.quant=val.quanti)
		ses.T_IC.IR<-ses(tstats$T_IC.IR,tstats$T_IC.IR_nm, val.quant=val.quanti)
		ses.T_PC.PR<-ses(tstats$T_PC.PR,tstats$T_PC.PR_nm, val.quant=val.quanti)
		
		#________________________________________
		#Condition to be significantly different from null models with respect to values of quantile choose
		cond.T_IP.IC.inf<-ses.T_IP.IC$ses<ses.T_IP.IC$ses.inf
		cond.T_IC.IR.inf<-ses.T_IC.IR$ses<ses.T_IC.IR$ses.inf
		cond.T_PC.PR.inf<-ses.T_PC.PR$ses<ses.T_PC.PR$ses.inf
			
		cond.T_IP.IC.sup<-ses.T_IP.IC$ses>ses.T_IP.IC$ses.sup
		cond.T_IC.IR.sup<-ses.T_IC.IR$ses>ses.T_IC.IR$ses.sup
		cond.T_PC.PR.sup<-ses.T_PC.PR$ses>ses.T_PC.PR$ses.sup
		
	}
	
	else if(nindex==4){
		#________________________________________
		#Calcul of standardised effect size
		
		ses.range<-ses(com.i$Range_obs, com.i$Null_range, val.quant=val.quanti)
		ses.mean<-ses(com.i$Mean_obs, com.i$Null_mean, val.quant=val.quanti)
		ses.CV<-ses(com.i$CV.NND_obs, com.i$Null_CV, val.quant=val.quanti)
		ses.kurtosis<-ses(com.i$kurtosis_obs, com.i$Null_kurtosis, val.quant=val.quanti)
			
		#________________________________________
		#Condition to be significantly different from null models with respect to values of quantile choose
			
		cond.range.inf<-t(ses.range$ses<ses.range$ses.inf)
		cond.mean.inf<-t(ses.mean$ses<ses.mean$ses.inf)
		cond.CV.NND.inf<-t(ses.CV$ses<ses.CV$ses.inf)
		cond.kurtosis.inf<-t(ses.kurtosis$ses<ses.kurtosis$ses.inf)
		
		cond.range.sup<-t(ses.range$ses>ses.range$ses.sup)
		cond.mean.sup<-t(ses.mean$ses>ses.mean$ses.sup)
		cond.CV.NND.sup<-t(ses.CV$ses>ses.CV$ses.sup)
		cond.kurtosis.sup<-t(ses.kurtosis$ses>ses.kurtosis$ses.sup)
		
	}
	
	else{}
	
	
	if(nindex==7 | nindex==4){
		ntr<-dim(com.i$Range_obs)[1]
	}
	
	else { ntr<-dim(tstats$T_PC.PR)[2] }
	
	if(method=="p.val"){
		mat_ass_com<-matrix(0, nrow=nplots*nindex, ncol=ntr)
		
		for (t in 1: ntr){
			for (s in 1:nplots){
						
				if(nindex==7 | nindex==4){
					#Range
					if(is.na(cond.range.inf[s,t])) 
						{mat_ass_com[s,t]<- 0}
					else if(cond.range.inf[s,t])
						{mat_ass_com[s,t]<- 1}
					else if(cond.range.sup[s,t])
						{mat_ass_com[s,t]<- -1}
					else{mat_ass_com[s,t]<- 0}
				
					#Mean
					if(is.na(cond.mean.inf[s,t]))
						{mat_ass_com[s+nplots,t]<- 0}
					else if(cond.mean.inf[s,t])
						{mat_ass_com[s+nplots,t]<- 1}
					else if(cond.mean.sup[s,t])
						{mat_ass_com[s+nplots,t]<- -1}
					else{mat_ass_com[s+nplots,t]<- 0}
					
					#CVNND
					if(is.na(cond.CV.NND.inf[s,t])) 
						{mat_ass_com[s+2*nplots,t]<- 0}
					else if(cond.CV.NND.inf[s,t])
						{mat_ass_com[s+2*nplots,t]<- 1}    
					else if(cond.CV.NND.sup[s,t])
						{mat_ass_com[s+2*nplots,t]<- -1}
					else{mat_ass_com[s+2*nplots,t]<- 0}
					
					#Kurtosis
					if(is.na(cond.kurtosis.inf[s,t])) 
						{mat_ass_com[s+3*nplots,t]<- 0}
					else if(cond.kurtosis.inf[s,t])
						{mat_ass_com[s+3*nplots,t]<- 1}    
					else if(cond.kurtosis.sup[s,t])
						  {mat_ass_com[s+3*nplots,t]<- -1}
					else{mat_ass_com[s+3*nplots,t]<- 0}
				}
				
				if(nindex==7){
					#T.IP_IC
					if(is.na(cond.T_IP.IC.inf[s,t])) 
						{mat_ass_com[s+4*nplots,t]<-0}
					else if(cond.T_IP.IC.inf[s,t])
						{mat_ass_com[s+4*nplots,t]<- 1}
					else if(cond.T_IP.IC.sup[s,t])
						{mat_ass_com[s+4*nplots,t]<- -1}
					else{mat_ass_com[s+4*nplots,t]<- 0}
					
					#T.IC_IR
					if(is.na(cond.T_IC.IR.inf[s,t])) 
						{mat_ass_com[s+5*nplots,t]<-  0}
					else if(cond.T_IC.IR.inf[s,t])
						{mat_ass_com[s+5*nplots,t]<-  1}
					else if(cond.T_IC.IR.sup[s,t])
						{mat_ass_com[s+5*nplots,t]<- -1}
					else{mat_ass_com[s+5*nplots,t]<-  0}
						
					#T.PC_PR
					if(is.na(cond.T_PC.PR.inf[s,t])) 
						{mat_ass_com[s+6*nplots,t]<-0}
					else if(cond.T_PC.PR.inf[s,t])
						{mat_ass_com[s+6*nplots,t]<- 1}
					else if(cond.T_PC.PR.sup[s,t])
						{mat_ass_com[s+6*nplots,t]<- -1}
					else{mat_ass_com[s+6*nplots,t]<- 0}
				}
				
				else if (nindex==3){
					#T.IP_IC
					if(is.na(cond.T_IP.IC.inf[s,t])) 
						{mat_ass_com[s,t]<-0}
					else if(cond.T_IP.IC.inf[s,t])
						{mat_ass_com[s,t]<-1}
					else if(cond.T_IP.IC.sup[s,t])
						{mat_ass_com[s,t]<- -1}
					else{mat_ass_com[s,t]<- 0}
					
					#T.IC_IR
					if(is.na(cond.T_IC.IR.inf[s,t])) 
						{mat_ass_com[s+nplots,t]<-0}
					else if(cond.T_IC.IR.inf[s,t])
						{mat_ass_com[s+nplots,t]<- 1}
					else if(cond.T_IC.IR.sup[s,t])
						{mat_ass_com[s+nplots,t]<- -1}
					else{mat_ass_com[s+nplots,t]<- 0}
						
					#T.PC_PR
					if(is.na(cond.T_PC.PR.inf[s,t])) 
						{mat_ass_com[s+2*nplots,t]<-0}
					else if(cond.T_PC.PR.inf[s,t])
						{mat_ass_com[s+2*nplots,t]<- 1}
					else if(cond.T_PC.PR.sup[s,t])
						{mat_ass_com[s+2*nplots,t]<- -1}
					else{mat_ass_com[s+2*nplots,t]<- 0}
				}
				else{}
			}
		}
	}
	
	else if(method=="ses"){
		mat_ass_com<-matrix(0, nrow=nplots*nindex, ncol=ntr)
		
		for (t in 1: ntr){
			for (s in 1:nplots){
				if(nindex==7){
					mat_ass_com[s,t]<-ses.range$ses[t,s]
					mat_ass_com[s+nplots,t]<-ses.mean$ses[t,s]
					mat_ass_com[s+2*nplots,t]<-ses.CV$ses[t,s]
					mat_ass_com[s+3*nplots,t]<-ses.kurtosis$ses[t,s]
					mat_ass_com[s+4*nplots,t]<-ses.T_IP.IC$ses[s,t]
					mat_ass_com[s+5*nplots,t]<-ses.T_IC.IR$ses[s,t]
					mat_ass_com[s+6*nplots,t]<-ses.T_PC.PR$ses[s,t]		
				}
				
				else if (nindex==3){
					mat_ass_com[s,t]<-ses.T_IP.IC$ses[s,t]
					mat_ass_com[s+nplots,t]<-ses.T_IC.IR$ses[s,t]
					mat_ass_com[s+2*nplots,t]<-ses.T_PC.PR$ses[s,t]
				}
				
				else if (nindex==4){
					mat_ass_com[s,t]<-ses.range$ses[t,s]
					mat_ass_com[s+nplots,t]<-ses.mean$ses[t,s]
					mat_ass_com[s+2*nplots,t]<-ses.CV$ses[t,s]
					mat_ass_com[s+3*nplots,t]<-ses.kurtosis$ses[t,s]
				}
				else{}
			}
		}
	}
	
	else{print("Error: method need to be `ses` or `p.val` ")}
	
	if(!is.null(tstats)) 
		{colnames(mat_ass_com)<- colnames(tstats$T_IP.IC)}
	else {colnames(mat_ass_com)<- colnames(com.i$namestraits)}

	if(nindex==4)
		{rownames(mat_ass_com)<-c(paste(unique(sites_names),"range", sep=" "), paste(unique(sites_names),"mean", sep=" "), paste(unique(sites_names),"CV_NND", sep=" "), paste(unique(sites_names),"kurtosis", sep=" "))}
	
	else if(nindex==3)
		{rownames(mat_ass_com)<-c(paste(unique(sites_names),"T_IP.IC", sep=" "), paste(unique(sites_names),"T_IC.IR", sep=" "), paste(unique(sites_names),"T_PC.PR", sep=" "))}
	
	else if(nindex==7)
	{rownames(mat_ass_com)<-c(paste(unique(sites_names),"range", sep=" "), paste(unique(sites_names),"mean", sep=" "), paste(unique(sites_names),"CV_NND", sep=" "), paste(unique(sites_names),"kurtosis", sep=" "), paste(unique(sites_names),"T_IP.IC", sep=" "), paste(unique(sites_names),"T_IC.IR", sep=" "), paste(unique(sites_names),"T_PC.PR", sep=" "))}
	
	else{}
	



return(mat_ass_com)

}

### Function to represent standardised effect size of all indices using null models
# index.list is a list of index associate with a list of corresponding null models in this order: [1] index 1 obs - [2] index 1 null model - [3] index 2 obs - [4] index 2 null model ...
#e.g. index.list<-list(T_IP.IC=res.finch$T_IP.IC, T_IP.IC_nm=res.finch$T_IP.IC_nm, T_PC.PR=res.finch$T_PC.PR, T_PC.PR_nm=res.finch$T_PC.PR_nm)
#observed matrix of values need to be of the same dimension

#You can transpose the observed matrix to represent either the ses by traits or by plots

plot.ses<-function(index.list, col.index=c("red","purple","green"), type="normal", add.conf=TRUE, color.cond=TRUE, val.quant=c(0.025,0.975), grid.v=TRUE, grid.h=TRUE, xlim=NULL, ylim=NULL){
	#possible type = "simple",  "simple_range", "normal" and "barplot"	
	
	namesindex.all<-names(index.list)
	nindex<-length(names(index.list))/2
	namesindex<-names(index.list)[seq(1,nindex*2, by=2)]
	namestraits<-colnames(index.list[[1]])
	namescommunity<-rownames(index.list[[1]])
	
	ncom<-dim(index.list[[1]])[1]
	ntr<-dim(index.list[[1]])[2]
	
	
	if(length(col.index)<nindex){
		col.index<-palette(rainbow(nindex))
	}
	
	#________________________________________
	#Calcul of standardised effect size
	
	res<-list()
	for (i in seq(1,nindex*2, by=2)){
		res[[eval(namesindex.all[i])]] <- ses(obs=index.list[[i]], nullmodel=index.list[[i+1]], val.quant=val.quant)
	}
	
	if(is.null(ylim)) { ylim=c(0,5.5+(nindex+1)*ntr)}
	if(is.null(xlim)) { xlim=c(min(c(unlist(res),-2), na.rm=T),max(c(unlist(res),2), na.rm=T))}
	
	par(mar=c(5, 7, 4, 2))
	plot(0, ylab="Traits",yaxt= "n", xlab="Standardized Effect Size", ylim=ylim, xlim=xlim, col="black", type="l")
	axis(side=2, seq(from=5.5, to=4.5+(nindex+1)*ntr, by=nindex+1)+(nindex+1)/2, labels=namestraits, las=1, cex.axis=0.7 ) 
	abline(v=0)
	
	if(grid.v==T) {
		range.<-max(c(unlist(res),2), na.rm=T)-min(c(unlist(res),-2), na.rm=T)
		
		vect.grid<-seq(min(unlist(res), na.rm=T),max(unlist(res), na.rm=T), by=round(range.,2)/9)
		for(j in vect.grid){
			abline(v=j, lty=2, col="lightgray")	
		}
	}
	
	if(grid.h==T) {
		for(j in seq(5.5,5.5+(nindex+1)*ntr)  ){
			abline(h=j, lty=2, col="lightgray")
		}
	}
	
	#________________________________________
	#plot : possible type = "simple", "simple_range", "normal" and "barplot"	
	
	#__________
	if(type=="simple"){
			
		for(i in 1:nindex){
			for(t in 1:ntr){
				
				if(color.cond==F){
					points(res[[eval(namesindex[i])]]$ses [,t], rep(5.5+(nindex+1)*t-i, times=ncom), pch=20, col=col.index[i])
				}
				
				if(color.cond==T){
					if(length(col.index)!=2*nindex) {col.index[(nindex+1):(nindex*2)]<-"grey"} 
					points(res[[eval(namesindex[i])]]$ses [,t], rep(5.5+(nindex+1)*t-i, times=ncom), pch=20, col=col.index[nindex+i])
					condition<-res[[eval(namesindex[i])]]$ses [,t] > res[[eval(namesindex[i])]]$ses.sup [,t]  |   res[[eval(namesindex[i])]]$ses [,t] < res[[eval(namesindex[i])]]$ses.inf [,t]
					condition[is.na(condition)]<-FALSE					
					
					if(sum(condition)>0){
						points(res[[eval(namesindex[i])]]$ses [condition,t], rep(5.5+(nindex+1)*t-i, times=sum(condition)), pch=20, col=col.index[i])
					}
				}
								
				points(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch=17, col=col.index[i])
												
				if(add.conf==T){
					points(mean(res[[eval(namesindex[i])]]$ses.sup [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="|", col=col.index[i])
					points(mean(res[[eval(namesindex[i])]]$ses.inf [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="|", col=col.index[i])
				}
		
				abline(seq(from=5.5, to=4.5+(nindex+1)*ntr, by=nindex+1)[t],b=0, lty=4, lwd=0.2)
				
			}		
		}
	}
	
	#__________
	else if(type=="simple_range"){
		
		for(i in 1:nindex){
			for(t in 1:ntr){
				
				if(color.cond==F){
					points(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch=17, col=col.index[i])
				}
				
				if(color.cond==T){
					if(length(col.index)!=2*nindex) {col.index[(nindex+1):(nindex*2)]<-"grey"} 
					points(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch=17, col=col.index[nindex+i])
					
					condition<-mean(res[[eval(namesindex[i])]]$ses [,t], na.rm=T) > mean(res[[eval(namesindex[i])]]$ses.sup [,t], na.rm=T)  |  mean( res[[eval(namesindex[i])]]$ses [,t], na.rm=T) < mean(res[[eval(namesindex[i])]]$ses.inf [,t], na.rm=T)
					if(sum(condition)>0){
						points(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch=17, col=col.index[i])
					}
				}
				
				if(add.conf==T){
					points(mean(res[[eval(namesindex[i])]]$ses.sup [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="|", col=col.index[i])
					points(mean(res[[eval(namesindex[i])]]$ses.inf [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="|", col=col.index[i])
				}
				
				segments(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T) + sd(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i , mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T) - sd(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, col=col.index[i])
				
				points(min(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="*", col=col.index[i])
				points(max(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="*", col=col.index[i])
				
				abline(seq(from=5.5, to=4.5+(nindex+1)*ntr, by=nindex+1)[t],b=0, lty=4, lwd=0.2)
		
			}		
		}
	}
	
	#__________
	else if(type=="normal"){
		
		
		for(i in 1:nindex){
			for(t in 1:ntr){
				
				if(color.cond==F){
					points(res[[eval(namesindex[i])]]$ses [,t], rep(5.5+(nindex+1)*t-i, times=ncom), pch=20, col=col.index[i])
				}
				
				if(color.cond==T){
					if(length(col.index)!=2*nindex) {col.index[(nindex+1):(nindex*2)]<-"grey"} 
					points(res[[eval(namesindex[i])]]$ses [,t], rep(5.5+(nindex+1)*t-i, times=ncom), pch=20, col=col.index[nindex+i])
					condition<-res[[eval(namesindex[i])]]$ses [,t] > res[[eval(namesindex[i])]]$ses.sup [,t]  |   res[[eval(namesindex[i])]]$ses [,t] < res[[eval(namesindex[i])]]$ses.inf [,t]
					condition[is.na(condition)]<-FALSE					
					
					if(sum(condition)>0){
						points(res[[eval(namesindex[i])]]$ses [condition,t], rep(5.5+(nindex+1)*t-i, times=sum(condition)), pch=20, col=col.index[i])
					}
				}
		
				if(add.conf==T){
					points(mean(res[[eval(namesindex[i])]]$ses.sup [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="|", col=col.index[i])
					points(mean(res[[eval(namesindex[i])]]$ses.inf [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="|", col=col.index[i])
				}
				
				segments(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T) + sd(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i , mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T) - sd(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, col=col.index[i])
				points(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch=10, col=col.index[i])
				
				abline(seq(from=5.5, to=4.5+(nindex+1)*ntr, by=nindex+1)[t],b=0, lty=4, lwd=0.2)
		
			}		
		}
		
	}

	#__________
	else if(type=="barplot"){
		for(i in 1:nindex){
			for(t in 1:ntr){
				
				segments(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, 0, 5.5+(nindex+1)*t-i, pch=17, col=col.index[i], lwd=8)
				segments(mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T) + sd(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i , mean(res[[eval(namesindex[i])]]$ses [,t],na.rm=T) - sd(res[[eval(namesindex[i])]]$ses [,t],na.rm=T), 5.5+(nindex+1)*t-i, col=col.index[i])
								
				if(add.conf==T){
					points(mean(res[[eval(namesindex[i])]]$ses.sup [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="|", col=col.index[i], cex=2)
					points(mean(res[[eval(namesindex[i])]]$ses.inf [,t],na.rm=T), 5.5+(nindex+1)*t-i, pch="|", col=col.index[i], cex=2)
				}
				
				abline(seq(from=5.5, to=4.5+(nindex+1)*ntr, by=nindex+1)[t],b=0, lty=4, lwd=0.2)
		
			}		
		}
		
	}
	else{print(paste("Error:",type,"is not a valid type of plot"))}
	
	legend("bottom", inset=.005, namesindex, fill=col.index, ncol=round(nindex/3) ,cex=0.6, bty="0")

	par(mar=c(5, 4, 4, 2) + 0.1) #return to default parameter
}




#Calcul of statistics (e.g. mean, range, CVNND and kurtosis) to test community assembly using null models
#For each statistic this function return observed value and correspondant Null distribution
#This function implement three null models which keep unchanged the number of individual per community
#Models 1 correspond to randomization of individual values whithin community
#Models 2 correspond to randomization of individual values whithin region
#Models 3 correspond to randomization of population values whithin region

#In most case, model 1 and 2 correspond to index at the individual level and the model 3 to index at the species (or any other aggregate variable like genus or family) level


com.index<-function(traits=NULL, index=NULL, namesindex=NULL, nullmodels=NULL, ind.plot=NULL, sp=NULL, nperm=99, printprogress=TRUE){
	
	nindex<-length(index)
	
	if(is.null(namesindex)) {  namesindex<-index }
	ntr<-dim(traits)[2]
	namestraits<-colnames(traits)
	
	traits<-traits[order(ind.plot),]
	ind.plot<-ind.plot[order(ind.plot)]
	sp<-sp[order(ind.plot)]
	
	name_sp_sites=paste(sp, ind.plot, sep="_")
	comm=NULL
	comm<-t(table(ind.plot,1:length(ind.plot)))
	
	S = colSums(comm>0)
	ncom=length(S)
	
	if(is.numeric(nperm)){
		######################################### 
		#### 	  Calcul of null models  	 ####
		######################################### 
		#Creation of three null models 
		if(printprogress==T){ print("creating null models")}
		
		if(sum(nullmodels==1)>0){
			#________________________________________
			#modèle nul 1: permutation des valeurs de traits des individus dans la communauté   
			
			traits.nm1<-list()
			
			for (t in 1: ntr){	
				traits.nm1[[eval(namestraits[t])]]<-matrix(NA, nrow=dim(traits)[1], ncol=nperm)
				perm_ind.plot<-list()
				
				for(n in 1:nperm){
					for(s in 1:  ncom) {
						perm_ind.plot[[s]]<-sample(traits[ind.plot==levels(ind.plot)[s], t], table(ind.plot)[s])
					}
					
					traits.nm1[[eval(namestraits[t])]][,n]<-unlist(perm_ind.plot)
				} 
				if(printprogress==T){
					print(paste("nm.1",round(t/ntr*100,2),"%")) 
				} 
			}
		}
		
		
		if(sum(nullmodels==2)>0){	
			#________________________________________
			#modèle nul 2: permutation des valeurs de traits des individus de la région    
			traits.nm2<-list()
			
			for (t in 1: ntr){	
				traits.nm2[[eval(namestraits[t])]]<-matrix(NA, nrow=dim(traits)[1], ncol=nperm)
				perm_ind.plot<-list()
				
				for(n in 1:nperm){
					for(s in 1:  ncom) {
						perm_ind.plot[[s]]<-sample(traits[, t], table(ind.plot)[s])
					}
					
					traits.nm2[[eval(namestraits[t])]][,n]<-unlist(perm_ind.plot)
				} 
				if(printprogress==T){
					print(paste("nm.2",round(t/ntr*100,2),"%")) 
				} 
			}
		}
		
		
		if(sum(nullmodels==3)>0){
			#________________________________________  
			#modèle nul 3: permutation des espèces au niveau de la région   
			traits.nm3<-list()
			traits_by_sp<-apply(traits,2,function(x) tapply(x,name_sp_sites,mean, na.rm=T))  
			traits_by_pop<-traits_by_sp[match(name_sp_sites,rownames(traits_by_sp)),]
			
			for (t in 1: ntr){	
				traits.nm3[[eval(namestraits[t])]]<-matrix(NA, nrow=dim(traits)[1], ncol=nperm)
				perm_ind.plot<-list()
				
				for(n in 1:nperm){
					for(s in 1:  ncom) {
						perm_ind.plot[[s]]<-sample(traits_by_pop, table(ind.plot)[s])
					}
					
					traits.nm3[[eval(namestraits[t])]][,n]<-unlist(perm_ind.plot)
				} 		
				if(printprogress==T){
					print(paste("nm.3",round(t/ntr*100,2),"%")) 
				} 
			}
		}
		
		
		######################################## 
		####	 Calcul of random values   	####
		######################################## 
		Null<-list()
		
		if(printprogress==T){print("calcul of null values using null models")}
		
		for(i in 1:nindex){
			if(nullmodels[i]==1){nm.bis<-traits.nm1[[1]]}
			else if(nullmodels[i]==2){nm.bis<-traits.nm2[[1]]}
			else if(nullmodels[i]==3){nm.bis<-traits.nm3[[1]]}
			else{print("nullmodels need 1, 2 or 3")}
			
			functionindex= eval(index[i])
			
			dim2<-dim(apply(nm.bis, 2, function (x) eval(parse(text=functionindex))))[1]
			Null[[eval(namesindex[i])]] <- array(NA, dim=c(ntr, dim2, nperm) )
			
			if(is.null(dim2)) {
				Null[[eval(namesindex[i])]] <- array(NA, dim=c(ntr, 1, nperm) )
			}
			
			for (t in 1: ntr){
			
				if(nullmodels[i]==1){nm<-traits.nm1[[t]]}
				else if(nullmodels[i]==2){nm<-traits.nm2[[t]]}
				else if(nullmodels[i]==3){nm<-traits.nm3[[t]]}
				else{print("nullmodels need 1, 2 or 3")}
				
				Null[[eval(namesindex[i])]] [t,,] <- apply(nm, 2, function (x) eval(parse(text=functionindex)))				
		
				if(printprogress==T){
					print(paste(eval(namesindex[i]), round(t/ntr*100,2),"%")) 
				} 
			}
		}
	}
		  
	######################################## 
	####	Calcul of observed values	####
	######################################## 
	obs<-list()
	
	if(printprogress==T){print("calcul of observed values")}
	
	for(i in 1:nindex){
			functionindex= eval(index[i])
			obs[[eval(namesindex[i])]] <- array(dim=c(ntr, dim(apply(traits, 2, function (x) eval(parse(text=functionindex))))[1]))
	
		if(nullmodels[i]==3) {
			traits.pop<-apply(traits, 2 , function (x) tapply(x, name_sp_sites, mean , na.rm=T))
			obs[[eval(namesindex[i])]] <-  apply(traits.pop, 2, function (x) eval(parse(text=functionindex)))
		}
		
		else if(nullmodels[i]==1  |  nullmodels[i]==2) {
			obs[[eval(namesindex[i])]] <- apply(traits, 2, function (x) eval(parse(text=functionindex)))
			#obs[[eval(namesindex[i])]] [ !is.finite(obs[[eval(namesindex[i])]] )]<-NA
		}
		if(printprogress==T){
			print(paste(round(i/nindex*100,2),"%")) 
		} 
	}
			
		
	######################################## 
	####		Create results list		####
	######################################## 
	
	com.index<-list()  
	com.index$obs<-obs
		
	if(is.numeric(nperm)){
		com.index$Null<-Null
	}
	
	com.index$list.index<-list()
	com.index$list.index.t<-list()
	name.com.index_list.index<-vector()
	
	for(i in 1:nindex){
		com.index$list.index.t[[seq(1,nindex*2,by=2)[i]]]<-t(obs[[i]])
		com.index$list.index[[seq(1,nindex*2,by=2)[i]]]<-obs[[i]]
		name.com.index_list.index[seq(1,nindex*2,by=2)[i]]<-names(obs)[i]
		
		if(is.numeric(nperm)){
			com.index$list.index[[seq(1,nindex*2,by=2)[i]+1]]<-Null[[i]]
			com.index$list.index.t[[seq(1,nindex*2,by=2)[i]+1]]<-Null[[i]]
			name.com.index_list.index[seq(1,nindex*2,by=2)[i]+1]<-paste(names(Null)[i], "nm", sep="_")
		}
	}
	
	names(com.index$list.index.t)<-name.com.index_list.index
	names(com.index$list.index)<-name.com.index_list.index
		
	com.index$sites_richness<-S
	com.index$namestraits<-namestraits
	
	return(com.index)
}
