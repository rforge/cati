#Function for Tstats Package					

### Function to calcul Tstats 
Tstats<-function(Ttraits, ind_plot, sp, nperm=NULL, printprogress=TRUE, p.value=TRUE){
	#6 variances: I: individual, P: population, C: community, R: region
	#IP; IC; IR; PC; PR; CR
	
	#Ttraits is the matrix of individual traits, ind_plot is the name of the plot in which the individu is (factor type), and sp is the species name of each individual
	
	names_sp_ind_plot<-as.factor(paste(sp, ind_plot, sep="@")) 
	Tplosp=unlist(strsplit(levels(names_sp_ind_plot),split="@"))[2*(1:nlevels(names_sp_ind_plot))]; names(Tplosp)=levels(names_sp_ind_plot);
	#Tplosp is the plot in wich the population is
  
  
	######################################## 
	####	Calcul of observed values	####
	######################################## 

	#________________________________________
	#Objects creation
	mean_IP<-matrix(nrow=nlevels(names_sp_ind_plot), ncol=ncol(Ttraits)); rownames(mean_IP)=levels(names_sp_ind_plot);
	mean_PC<-matrix(nrow=nlevels(ind_plot), ncol=ncol(Ttraits))
	var_IP<-matrix(nrow=nlevels(names_sp_ind_plot), ncol=ncol(Ttraits))
	var_PC<-matrix(nrow=nlevels(ind_plot), ncol=ncol(Ttraits))
	var_CR<-vector()
	var_IC<-matrix(nrow=nlevels(ind_plot), ncol=ncol(Ttraits))
	var_PR<-vector()
	var_IR<-vector()
	T_IP.IC<-matrix(nrow=nlevels(ind_plot), ncol=ncol(Ttraits))
	T_IC.IR<-matrix(nrow=nlevels(ind_plot), ncol=ncol(Ttraits))
	T_PC.PR<-matrix(nrow=nlevels(ind_plot), ncol=ncol(Ttraits))
  
	for (t in 1: ncol(Ttraits)){
		mean_IP[,t]<-tapply(Ttraits[,t], names_sp_ind_plot  ,mean, na.rm=T)
		mean_PC[,t]<-tapply(mean_IP[,t], Tplosp , mean, na.rm=T)
		
		var_IP[,t]<-tapply(Ttraits[,t], names_sp_ind_plot, var, na.rm=T)
		var_PC[,t]<-tapply(mean_IP[,t], Tplosp  ,var, na.rm=T)
		var_CR[t]<-var(mean_PC[,t], na.rm=T)
		var_IC[,t]<-tapply(Ttraits[,t], ind_plot  ,var, na.rm=T)
		var_PR[t]<-var(as.vector(mean_IP[,t]), na.rm=T)
		var_IR[t]<-var(Ttraits[,t], na.rm=T)
		  
		for(s in 1 : nlevels(ind_plot)){
			T_IP.IC[s,t]<-mean(var_IP[grepl(levels(ind_plot)[s],Tplosp),t], na.rm=T)/var_IC[s,t]
			T_IC.IR[s,t]<-var_IC[s,t]/var_IR[t]
			T_PC.PR[s,t]<-var_PC[s,t]/var_PR[t]
		}
	}
	
	#________________________________________
	
	######################################### 
	#### 	  Creating null models  	 ####
	######################################### 
	
	if(is.numeric(nperm)){
		
		var_IP_nm1<-array(dim=c(nperm,ncol(Ttraits),nrow=length(Tplosp)))
		var_PC_nm3<-array(dim=c(nperm,ncol(Ttraits),nlevels(ind_plot)))
		var_IC_nm1<-array(dim=c(nperm,ncol(Ttraits),nlevels(ind_plot)))
		var_IC_nm2<-array(dim=c(nperm,ncol(Ttraits),nlevels(ind_plot)))
		var_PR_nm3<-array(dim=c(nperm,ncol(Ttraits)))
		var_IR_nm2<-array(dim=c(nperm,ncol(Ttraits)))
       
		mean_IP_nm3<-array(dim=c(nperm,ncol(Ttraits),length(Tplosp)))
		mean_PC_nm3<-array(dim=c(nperm,ncol(Ttraits),nlevels(ind_plot)))
       
		Ttraits.nm1<-list()
		Ttraits.nm2<-list()
		Ttraits.nm3<-list()
              
		T_IP.IC_nm1<-array(dim=c(nperm,ncol(Ttraits),nlevels(ind_plot)))
		T_IC.IR_nm2<-array(dim=c(nperm,ncol(Ttraits),nlevels(ind_plot)))
		T_PC.PR_nm3<-array(dim=c(nperm,ncol(Ttraits),nlevels(ind_plot)))
      
       
		#Creation of three null models 
		print("creating null models")
		
		#________________________________________
		#modèle nul 1: permutation des valeurs de traits des individus dans la communauté   
		for (t in 1: ncol(Ttraits)){
			Ttraits.nm1[[t]]<-list()
			for(s in 1:  nlevels(ind_plot)) {
				Ttraits.nm1[[t]][[s]]<-list()
				for(i in 1:nperm){
					if (length(Ttraits[ind_plot==levels(ind_plot)[s], t]) != 1) {
						perm_ind_plot1<-sample(Ttraits[ind_plot==levels(ind_plot)[s], t], table(ind_plot)[s])
						Ttraits.nm1[[t]][[s]][[i]]<-perm_ind_plot1
					}
					else {Ttraits.nm1[[t]][[s]][[i]]<-"NA"}
				}
			} 
			if(printprogress==T){print(paste(round(t/ncol(Ttraits)/3*100,2),"%")) } else {}
		}
		
		#________________________________________
		#modèle nul 2: permutation des valeurs de traits des individus de la région    
		for (t in 1: ncol(Ttraits)){
			Ttraits.nm2[[t]]<-list()
			for(s in 1:  nlevels(ind_plot)) {
				Ttraits.nm2[[t]][[s]]<-list()
				for(i in 1:nperm){
					perm_ind_plot2<-sample(Ttraits[, t], table(ind_plot)[s])
					Ttraits.nm2[[t]][[s]][[i]]<-perm_ind_plot2
				}
			}
			if(printprogress==T){print(paste(round(33.3+t/ncol(Ttraits)/3*100, 2),"%"))} else {}
		}
		
		#________________________________________  
		#modèle nul 3: permutation des espèces au niveau de la région   
		Ttraits_by_sp<-apply(Ttraits,2,function(x) tapply(x,names_sp_ind_plot,mean))  
		Ttraits_by_pop<-Ttraits_by_sp[match(names_sp_ind_plot,rownames(Ttraits_by_sp)),]
		#Ttraits_by_sp<-aggregate(Ttraits, by = list(names_sp_ind_plot), mean, na.rm = T)[,-1] 
				
		for (t in 1: ncol(Ttraits)){
			Ttraits.nm3[[t]]<-list()
			for(s in 1:  nlevels(ind_plot)){
				Ttraits.nm3[[t]][[s]]<-list()
				for(i in 1:nperm){
					perm_ind_plot3<-sample(Ttraits_by_pop, table(ind_plot)[s])
					Ttraits.nm3[[t]][[s]][[i]]<-perm_ind_plot3
				}
			} 
			if(printprogress==T){print(paste(round(66.6+t/ncol(Ttraits)/3*100, 2),"%"))} else {}
		}
		
		#________________________________________
	
		######################################### 
		#### calcul of Tstats on null models ####
		######################################### 

		print("calcul of Tstats using null models")
		
		yy<-length(names_sp_ind_plot)
		for (t in 1: ncol(Ttraits)){
			for(i in 1:nperm){ 
				mean_IP_nm3[i,t,]<-tapply(unlist(Ttraits.nm3[[t]])[(1+(i-1)*yy) : (i*yy)], names_sp_ind_plot  ,function(x) mean(x, na.rm=T))
				mean_PC_nm3[i,t,]<-tapply(mean_IP_nm3[i,t,], Tplosp, mean, na.rm=T)
			}
			if(printprogress==T){print(paste(round(t/ncol(Ttraits)/3*100, 2),"%"))} else {}
		} 
		   
		   
		for (t in 1: ncol(Ttraits)){
			for(i in 1:nperm){
				var_IP_nm1[i,t,]<-tapply(unlist(Ttraits.nm1[[t]])[(1+(i-1)*yy) : (i*yy)], names_sp_ind_plot  ,function(x) var(x, na.rm=T))
				var_PC_nm3[i,t,]<-tapply(mean_IP_nm3[i,t,], Tplosp  ,var, na.rm=T)
				var_IC_nm1[i,t,]<-tapply(unlist(Ttraits.nm1[[t]])[(1+(i-1)*yy) : (i*yy)], ind_plot  ,function(x) var(x, na.rm=T))
				var_IC_nm2[i,t,]<-tapply(unlist(Ttraits.nm2[[t]])[(1+(i-1)*yy) : (i*yy)], ind_plot  ,function(x) var(x, na.rm=T))
				var_PR_nm3[i,t]<-var(as.vector(mean_IP_nm3[i,t,]), na.rm=T)
				var_IR_nm2[i,t]<-var(unlist(Ttraits.nm2[[t]])[(1+(i-1)*yy) : (i*yy)], na.rm=T)
			}
			if(printprogress==T){print(paste(round(33.3+t/ncol(Ttraits)/3*100, 2),"%"))} else {}
		}
		 
		   
		for (t in 1: ncol(Ttraits)){
			for(i in 1:nperm){
				for(s in 1 : nlevels(ind_plot)){
					T_IP.IC_nm1[i,t,s]<-mean(var_IP_nm1[i,t,grepl(levels(ind_plot)[s],Tplosp)], na.rm=T)/var_IC_nm1[i,t,s] 
					T_IC.IR_nm2[i,t,s]<-var_IC_nm2[i,t,s]/var_IR_nm2[i,t]
					T_PC.PR_nm3[i,t,s]<-var_PC_nm3[i,t,s]/var_PR_nm3[i,t]
				}
			} 
			if(printprogress==T){print(paste(round(66.6+t/ncol(Ttraits)/3*100, 2),"%"))} else {}
		}       
		      
	}#end of calcul of Tstats using null models
         
	colnames(T_IP.IC)<-colnames(Ttraits)
    colnames(T_IC.IR)<-colnames(Ttraits)
    colnames(T_PC.PR)<-colnames(Ttraits)
  
	if(is.numeric(nperm)){
		colnames(T_IP.IC_nm1)<-colnames(Ttraits)
		colnames(T_IC.IR_nm2)<-colnames(Ttraits)
		colnames(T_PC.PR_nm3)<-colnames(Ttraits)
	}
	
	rownames(T_IP.IC)<-levels(as.factor(Tplosp))
  	rownames(T_IC.IR)<-levels(as.factor(Tplosp))
 	rownames(T_PC.PR)<-levels(as.factor(Tplosp))
 
	
	#________________________________________
    res<-list()
    res$T_IP.IC<-T_IP.IC
    res$T_IC.IR<-T_IC.IR
    res$T_PC.PR<-T_PC.PR
    res$var_IP<-var_IP
    res$var_PC<-var_PC
    res$var_CR<-var_CR
	res$var_IC<-var_IC
    res$var_PR<-var_PR
    res$var_IR<-var_IR
	
	if(is.numeric(nperm)){	 
		res$T_IP.IC_nm<-T_IP.IC_nm1
       	res$T_IC.IR_nm<-T_IC.IR_nm2
        res$T_PC.PR_nm<-T_PC.PR_nm3
    }   
    else{}
 	
 	#________________________________________
 	
 	######################################### 
	####		 calcul of p.value		 ####
	######################################### 
 
 	print("calcul of p.value")
 	
 	if(p.value==T){
		p.valueT_IP.IC.sup<-matrix(ncol=ncol(Ttraits), nrow= nlevels(ind_plot))
		p.valueT_IC.IR.sup<-matrix(ncol=ncol(Ttraits), nrow= nlevels(ind_plot))
		p.valueT_PC.PR.sup<-matrix(ncol=ncol(Ttraits), nrow= nlevels(ind_plot))
		
		p.valueT_IP.IC.inf<-matrix(ncol=ncol(Ttraits), nrow= nlevels(ind_plot))
		p.valueT_IC.IR.inf<-matrix(ncol=ncol(Ttraits), nrow= nlevels(ind_plot))
		p.valueT_PC.PR.inf<-matrix(ncol=ncol(Ttraits), nrow= nlevels(ind_plot))
		
		for (t in 1: ncol(Ttraits)){
			for(s in 1:  nlevels(ind_plot)){
 				p.valueT_IP.IC.sup[s,t]<-sum(res$T_IP.IC[s,t]<res$T_IP.IC_nm[,t,s], na.rm=T)/(1+length(res$T_IP.IC_nm[,t,s]))
 				p.valueT_IC.IR.sup[s,t]<-sum(res$T_IC.IR[s,t]<res$T_IC.IR_nm[,t,s], na.rm=T)/(1+length(res$T_IC.IR_nm[,t,s]))
 				p.valueT_PC.PR.sup[s,t]<-sum(res$T_PC.PR[s,t]<res$T_PC.PR_nm[,t,s], na.rm=T)/(1+length(res$T_PC.PR_nm[,t,s]))
		
				p.valueT_IP.IC.inf[s,t]<-sum(res$T_IP.IC[s,t]>res$T_IP.IC_nm[,t,s], na.rm=T)/(1+length(res$T_IP.IC_nm[,t,s]))
				p.valueT_IC.IR.inf[s,t]<-sum(res$T_IC.IR[s,t]>res$T_IC.IR_nm[,t,s], na.rm=T)/(1+length(res$T_IC.IR_nm[,t,s]))
				p.valueT_PC.PR.inf[s,t]<-sum(res$T_PC.PR[s,t]>res$T_PC.PR_nm[,t,s], na.rm=T)/(1+length(res$T_PC.PR_nm[,t,s]))
			}
		}	
	    
		colnames(p.valueT_IP.IC.sup)<-colnames(Ttraits)
		colnames(p.valueT_IC.IR.sup)<-colnames(Ttraits)
		colnames(p.valueT_PC.PR.sup)<-colnames(Ttraits)
		
		rownames(p.valueT_IP.IC.sup)<-levels(Tplosp)
		rownames(p.valueT_IC.IR.sup)<-levels(Tplosp)
		rownames(p.valueT_PC.PR.sup)<-levels(Tplosp)
  	
		colnames(p.valueT_IP.IC.inf)<-colnames(Ttraits)
		colnames(p.valueT_IC.IR.inf)<-colnames(Ttraits)
		colnames(p.valueT_PC.PR.inf)<-colnames(Ttraits)
		
		rownames(p.valueT_IP.IC.inf)<-levels(Tplosp)
		rownames(p.valueT_IC.IR.inf)<-levels(Tplosp)
		rownames(p.valueT_PC.PR.inf)<-levels(Tplosp)

 	}
 	else{}
    
    #________________________________________
    	
	if(p.value==T) {
		res$pval_T_IP.IC.inf<-p.valueT_IP.IC.inf
		res$pval_T_IC.IR.inf<-p.valueT_IC.IR.inf
		res$pval_T_PC.PR.inf<-p.valueT_PC.PR.inf
		
		res$pval_T_IP.IC.sup<-p.valueT_IP.IC.sup
		res$pval_T_IC.IR.sup<-p.valueT_IC.IR.sup
		res$pval_T_PC.PR.sup<-p.valueT_PC.PR.sup
	}
	else{}
	
   
    
    return(res)
}


### Function to represent standardised effect size of Tstats using null models
plot.ses.Tstats<-function(tstats=NULL, val.quant=c(0.025,0.975), col.Tstats=c("red","purple","green"), type="normal", add.conf=TRUE){
	#possible type = "color_cond", "simple", "simple_sd", "normal" and "barplot"	
	
	#________________________________________
	#Calcul of standardised effect size
	ses.T_IP.IC<-(tstats$T_IP.IC-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_IC.IR<-(tstats$T_IC.IR-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_PC.PR<-(tstats$T_PC.PR-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T))
	  
	ses.T_IP.IC.inf<-(apply(tstats$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_IC.IR.inf<-(apply(tstats$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_PC.PR.inf<-(apply(tstats$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T))
	  
	ses.T_IP.IC.sup<-(apply(tstats$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_IC.IR.sup<-(apply(tstats$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_PC.PR.sup<-(apply(tstats$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T))
	
	#________________________________________
	#Condition to be significantly different from null models with respect to values of quantile choose
	cond.T_IP.IC.inf<-ses.T_IP.IC<ses.T_IP.IC.inf
	cond.T_IC.IR.inf<-ses.T_IC.IR<ses.T_IC.IR.inf
	cond.T_PC.PR.inf<-ses.T_PC.PR<ses.T_PC.PR.inf
		
	cond.T_IP.IC.sup<-ses.T_IP.IC>ses.T_IP.IC.sup
	cond.T_IC.IR.sup<-ses.T_IC.IR>ses.T_IC.IR.sup
	cond.T_PC.PR.sup<-ses.T_PC.PR>ses.T_PC.PR.sup
	
	all=c(ses.T_IP.IC,ses.T_IC.IR,ses.T_PC.PR)
	par(mar=c(5, 7, 4, 2))
	plot(0, ylab="Traits",yaxt= "n", xlab="Tstats Standardized Effect Size", ylim=c(5*dim(tstats$T_IP.IC)[2],0), xlim=c(min(all, na.rm=T),max(all, na.rm=T)), col="black", type="l")
	axis(side=2, seq(from=5.5, to=4*dim(tstats$T_IP.IC)[2]+1.5, by=4), labels=colnames(tstats$T_IP.IC), las=1, cex.axis=0.7 ) 
	legend("bottom", inset=.005, title="Tstats", c("T_IP.IC","T_IC.IR","T_PC.PR"), fill=col.Tstats, horiz=TRUE, cex=0.7, bty="n")
	
	#________________________________________
	#plot : possible type = "color_cond", "simple", "simple_range", "normal" and "barplot"	
	
	#__________
	if(type=="color_cond"){
		
		if(length(col.Tstats)==3) {col.Tstats[4:6]<-"grey"} 
		if(length(col.Tstats)!=6) {print("Warnings: plot type color_cond need 3 or 6 colors in the argument col.Tstats")}
				
		for(t in 1:dim(tstats$T_IP.IC)[2]){
			points(ses.T_IP.IC[,t], rep(t*4, times=dim(tstats$T_IP.IC)[1]), pch=20, col=col.Tstats[4])
			points(ses.T_IC.IR[,t], rep(t*4+1, times=dim(tstats$T_IP.IC)[1]), pch=20, col=col.Tstats[5])
			points(ses.T_PC.PR[,t], rep(t*4+2, times=dim(tstats$T_IP.IC)[1]), pch=20, col=col.Tstats[6])
			
			points(mean(ses.T_IP.IC[,t], na.rm=T), t*4, pch=17, col=col.Tstats[1])
			points(mean(ses.T_IC.IR[,t], na.rm=T), t*4+1,pch=17, col=col.Tstats[2])
			points(mean(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch=17, col=col.Tstats[3])
			
			points(ses.T_IP.IC[,t][cond.T_IP.IC.inf[,t]], rep(t*4,times=length(ses.T_IP.IC.inf[,t][cond.T_IP.IC.inf[,t]])),pch=16, col=col.Tstats[1])
			points(ses.T_IC.IR[,t][cond.T_IC.IR.inf[,t]], rep(t*4+1,times=length(ses.T_IC.IR.inf[,t][cond.T_IC.IR.inf[,t]])), pch=16, col=col.Tstats[2])
			points(ses.T_PC.PR[,t][cond.T_PC.PR.inf[,t]], rep(t*4+2, times=length(ses.T_PC.PR.inf[,t][cond.T_PC.PR.inf[,t]])), pch=16, col=col.Tstats[3])
			
			points(ses.T_IP.IC[,t][cond.T_IP.IC.sup[,t]], rep(t*4,times=length(ses.T_IP.IC.sup[,t][cond.T_IP.IC.sup[,t]])), pch=16, col=col.Tstats[1])
			points(ses.T_IC.IR[,t][cond.T_IC.IR.sup[,t]], rep(t*4+1,times=length(ses.T_IC.IR.sup[,t][cond.T_IC.IR.sup[,t]])), pch=16, col=col.Tstats[2])
			points(ses.T_PC.PR[,t][cond.T_PC.PR.sup[,t]], rep(t*4+2, times=length(ses.T_PC.PR.sup[,t][cond.T_PC.PR.sup[,t]])), pch=16, col=col.Tstats[3])
			
			if (length(ses.T_IP.IC.inf[,t][cond.T_IP.IC.inf[,t]])>0) 	{text(ses.T_IP.IC[,t][cond.T_IP.IC.inf[,t]], rep(t*4-2,times=length(ses.T_IP.IC.inf[,t][cond.T_IP.IC.inf[,t]])), labels=rownames(tstats$T_IP.IC)[cond.T_IP.IC.inf[,t]], cex=0.6, srt=45, col=col.Tstats[1], pos=3)}
			if (length(ses.T_IC.IR.inf[,t][cond.T_IC.IR.inf[,t]])>0) 	{text(ses.T_IC.IR[,t][cond.T_IC.IR.inf[,t]], rep(t*4,times=length(ses.T_IC.IR.inf[,t][cond.T_IC.IR.inf[,t]])), labels=rownames(tstats$T_IP.IC)[cond.T_IC.IR.inf[,t]], cex=0.6, srt=45, col=col.Tstats[2], pos=2)}
			if (length(ses.T_PC.PR.inf[,t][cond.T_PC.PR.inf[,t]])>0) 	{text(ses.T_PC.PR[,t][cond.T_PC.PR.inf[,t]], rep(t*4+4, times=length(ses.T_PC.PR.inf[,t][cond.T_PC.PR.inf[,t]])), labels=rownames(tstats$T_IP.IC)[cond.T_PC.PR.inf[,t]], cex=0.6, srt=45, col=col.Tstats[3], pos=1)}
			
			if (length(ses.T_IP.IC.sup[,t][cond.T_IP.IC.sup[,t]])>0) 	{text(ses.T_IP.IC[,t][cond.T_IP.IC.sup[,t]], rep(t*4-2,times=length(ses.T_IP.IC.sup[,t][cond.T_IP.IC.sup[,t]])), labels=rownames(tstats$T_IP.IC)[cond.T_IP.IC.sup[,t]], cex=0.6, srt=45, col=col.Tstats[1], pos=3)}
			if (length(ses.T_IC.IR.sup[,t][cond.T_IC.IR.sup[,t]])>0) 	{text(ses.T_IC.IR[,t][cond.T_IC.IR.sup[,t]], rep(t*4,times=length(ses.T_IC.IR.sup[,t][cond.T_IC.IR.sup[,t]])), labels=rownames(tstats$T_IP.IC)[cond.T_IC.IR.sup[,t]], cex=0.6, srt=45, col=col.Tstats[2], pos=2)}
			if (length(ses.T_PC.PR.sup[,t][cond.T_PC.PR.sup[,t]])>0) 	{text(ses.T_PC.PR[,t][cond.T_PC.PR.sup[,t]], rep(t*4+4, times=length(ses.T_PC.PR.sup[,t][cond.T_PC.PR.sup[,t]])), labels=rownames(tstats$T_IP.IC)[cond.T_PC.PR.sup[,t]], cex=0.6, srt=45, col=col.Tstats[3], pos=1)}
			
			abline(a=t*4+3,b=0, lty=4, lwd=0.2)
		}
	}

	#__________
	else if(type=="simple"){
	
		for(t in 1:dim(tstats$T_IP.IC)[2]){
			
			points(ses.T_IP.IC[,t], rep(t*4, times=dim(tstats$T_IP.IC)[1]), pch=20, col=col.Tstats[1])
			points(ses.T_IC.IR[,t], rep(t*4+1, times=dim(tstats$T_IP.IC)[1]), pch=20, col=col.Tstats[2])
			points(ses.T_PC.PR[,t], rep(t*4+2, times=dim(tstats$T_IP.IC)[1]), pch=20, col=col.Tstats[3])
			
			points(mean(ses.T_IP.IC[,t], na.rm=T), t*4, pch=17, col=col.Tstats[1])
			points(mean(ses.T_IC.IR[,t], na.rm=T), t*4+1,pch=17, col=col.Tstats[2])
			points(mean(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch=17, col=col.Tstats[3])
			
			points(mean(ses.T_IP.IC.sup[,t], na.rm=T), t*4, pch="|", col=col.Tstats[1])
			points(mean(ses.T_IC.IR.sup[,t], na.rm=T), t*4+1,pch="|", col=col.Tstats[2])
			points(mean(ses.T_PC.PR.sup[,t], na.rm=T), t*4+2, pch="|", col=col.Tstats[3])
			
			points(mean(ses.T_IP.IC.inf[,t], na.rm=T), t*4, pch="|", col=col.Tstats[1])
			points(mean(ses.T_IC.IR.inf[,t], na.rm=T), t*4+1,pch="|", col=col.Tstats[2])
			points(mean(ses.T_PC.PR.inf[,t], na.rm=T), t*4+2, pch="|", col=col.Tstats[3])	
			
			abline(a=t*4+3,b=0, lty=4, lwd=0.2)
		}
	}
	
	#__________
	else if(type=="simple_sd"){
				
		for(t in 1:dim(tstats$T_IP.IC)[2]){
			
			points(mean(ses.T_IP.IC[,t], na.rm=T), t*4, pch=17, col=col.Tstats[1])
			points(mean(ses.T_IC.IR[,t], na.rm=T), t*4+1,pch=17, col=col.Tstats[2])
			points(mean(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch=17, col=col.Tstats[3])
			
			segments(mean(ses.T_IP.IC[,t], na.rm=T)+sd(ses.T_IP.IC[,t], na.rm=T), t*4, mean(ses.T_IP.IC[,t], na.rm=T)-sd(ses.T_IP.IC[,t], na.rm=T), t*4, pch=17, col=col.Tstats[1])
			segments(mean(ses.T_IC.IR[,t], na.rm=T)+sd(ses.T_IC.IR[,t], na.rm=T), t*4+1, mean(ses.T_IC.IR[,t], na.rm=T)-sd(ses.T_IC.IR[,t], na.rm=T), t*4+1, pch=17, col=col.Tstats[2])
			segments(mean(ses.T_PC.PR[,t], na.rm=T)+sd(ses.T_PC.PR[,t], na.rm=T), t*4+2, mean(ses.T_PC.PR[,t], na.rm=T)-sd(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch=17, col=col.Tstats[3])
					
			points(mean(ses.T_IP.IC.sup[,t], na.rm=T), t*4, pch="|", col=col.Tstats[1])
			points(mean(ses.T_IC.IR.sup[,t], na.rm=T), t*4+1,pch="|", col=col.Tstats[2])
			points(mean(ses.T_PC.PR.sup[,t], na.rm=T), t*4+2, pch="|", col=col.Tstats[3])
			
			points(mean(ses.T_IP.IC.inf[,t], na.rm=T), t*4, pch="|", col=col.Tstats[1])
			points(mean(ses.T_IC.IR.inf[,t], na.rm=T), t*4+1,pch="|", col=col.Tstats[2])
			points(mean(ses.T_PC.PR.inf[,t], na.rm=T), t*4+2, pch="|", col=col.Tstats[3])	
			
			abline(a=t*4+3,b=0, lty=4, lwd=0.2)
		}	
	}
	
	#__________
	else if(type=="normal"){
		
		for(t in 1:dim(tstats$T_IP.IC)[2]){
		
			points(mean(ses.T_IP.IC[,t], na.rm=T), t*4, pch=17, col=col.Tstats[1])
			points(mean(ses.T_IC.IR[,t], na.rm=T), t*4+1,pch=17, col=col.Tstats[2])
			points(mean(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch=17, col=col.Tstats[3])
			
			
			segments(mean(ses.T_IP.IC[,t], na.rm=T)+sd(ses.T_IP.IC[,t], na.rm=T), t*4, mean(ses.T_IP.IC[,t], na.rm=T)-sd(ses.T_IP.IC[,t], na.rm=T), t*4, pch=17, col=col.Tstats[1])
			segments(mean(ses.T_IC.IR[,t], na.rm=T)+sd(ses.T_IC.IR[,t], na.rm=T), t*4+1, mean(ses.T_IC.IR[,t], na.rm=T)-sd(ses.T_IC.IR[,t], na.rm=T), t*4+1, pch=17, col=col.Tstats[2])
			segments(mean(ses.T_PC.PR[,t], na.rm=T)+sd(ses.T_PC.PR[,t], na.rm=T), t*4+2, mean(ses.T_PC.PR[,t], na.rm=T)-sd(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch=17, col=col.Tstats[3])
	   
			points(min(ses.T_IP.IC[,t], na.rm=T), t*4, pch="*", col=col.Tstats[1])
			points(min(ses.T_IC.IR[,t], na.rm=T), t*4+1,pch="*", col=col.Tstats[2])
			points(min(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch="*", col=col.Tstats[3])
			
			points(max(ses.T_IP.IC[,t], na.rm=T), t*4, pch="*", col=col.Tstats[1])
			points(max(ses.T_IC.IR[,t], na.rm=T), t*4+1,pch="*", col=col.Tstats[2])
			points(max(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch="*", col=col.Tstats[3])
			
			abline(a=t*4+3,b=0, lty=4, lwd=0.2)
		}
		
		if(add.conf==T){
			points(colMeans(ses.T_IP.IC.sup, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4, type="l", col=col.Tstats[1])
			points(colMeans(ses.T_IC.IR.sup, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4+1, type="l", col=col.Tstats[2])
			points(colMeans(ses.T_PC.PR.sup, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4+2, type="l", col=col.Tstats[3])
			
			points(colMeans(ses.T_IP.IC.inf, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4, type="l", col=col.Tstats[1])
			points(colMeans(ses.T_IC.IR.inf, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4+1, type="l", col=col.Tstats[2])
			points(colMeans(ses.T_PC.PR.inf, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4+2, type="l", col=col.Tstats[3])  
		}	
		else {}
	}

	#__________
	else if(type=="barplot"){
	
		for(t in 1:dim(tstats$T_IP.IC)[2]){
		
			segments(mean(ses.T_IP.IC[,t], na.rm=T), t*4  , 0, t*4, pch=17, col=col.Tstats[1], lwd=8)
			segments(mean(ses.T_IC.IR[,t], na.rm=T), t*4+1, 0, t*4+1, pch=17, col=col.Tstats[2], lwd=8)
			segments(mean(ses.T_PC.PR[,t], na.rm=T), t*4+2, 0, t*4+2, pch=17, col=col.Tstats[3], lwd=8)
			
			segments(mean(ses.T_IP.IC[,t], na.rm=T)+sd(ses.T_IP.IC[,t], na.rm=T), t*4, mean(ses.T_IP.IC[,t], na.rm=T)-sd(ses.T_IP.IC[,t], na.rm=T), t*4, pch=17, col=col.Tstats[1])
			segments(mean(ses.T_IC.IR[,t], na.rm=T)+sd(ses.T_IC.IR[,t], na.rm=T), t*4+1, mean(ses.T_IC.IR[,t], na.rm=T)-sd(ses.T_IC.IR[,t], na.rm=T), t*4+1, pch=17, col=col.Tstats[2])
			segments(mean(ses.T_PC.PR[,t], na.rm=T)+sd(ses.T_PC.PR[,t], na.rm=T), t*4+2, mean(ses.T_PC.PR[,t], na.rm=T)-sd(ses.T_PC.PR[,t], na.rm=T), t*4+2, pch=17, col=col.Tstats[3])

			abline(a=t*4+3,b=0, lty=4, lwd=0.2)
		}
		
		if(add.conf==T){
			points(colMeans(ses.T_IP.IC.sup, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4, type="l", col=col.Tstats[1])
			points(colMeans(ses.T_IC.IR.sup, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4+1, type="l", col=col.Tstats[2])
			points(colMeans(ses.T_PC.PR.sup, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4+2, type="l", col=col.Tstats[3])
			
			points(colMeans(ses.T_IP.IC.inf, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4, type="l", col=col.Tstats[1])
			points(colMeans(ses.T_IC.IR.inf, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4+1, type="l", col=col.Tstats[2])
			points(colMeans(ses.T_PC.PR.inf, na.rm=T), (1:dim(tstats$T_IP.IC)[2])*4+2, type="l", col=col.Tstats[3])  
		}
		else {}
	}
	else{print(paste("Error:",type,"is not a valid type of plot"))}
	
	par(mar=c(5, 4, 4, 2) + 0.1) #return to default parameter
}


### Function to represent correlations between Tstats
plot.cor.Tstats<-function(tstats=NULL, val.quant=c(0.025,0.975), add.text=FALSE, bysites=FALSE, col.obj=NULL, plot.ask=TRUE) {
	
	oldpar<-par(no.readonly = TRUE)
	par(ask=plot.ask)
	
	#________________________________________
	ses.T_IP.IC.moy<-t(colMeans((tstats$T_IP.IC-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T)), na.rm=T))
	ses.T_IC.IR.moy<-t(colMeans((tstats$T_IC.IR-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T)), na.rm=T))
	ses.T_PC.PR.moy<-t(colMeans((tstats$T_PC.PR-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T)), na.rm=T))
	
	ses.T_IP.IC<-t((tstats$T_IP.IC-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T)))
	ses.T_IC.IR<-t((tstats$T_IC.IR-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T)))
	ses.T_PC.PR<-t((tstats$T_PC.PR-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T)))
	
	ses.T_IP.IC.inf<-t((apply(tstats$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T)))
	ses.T_IC.IR.inf<-t((apply(tstats$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T)))
	ses.T_PC.PR.inf<-t((apply(tstats$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T)))
	
	ses.T_IP.IC.sup<-t((apply(tstats$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T)))
	ses.T_IC.IR.sup<-t((apply(tstats$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T)))
	ses.T_PC.PR.sup<-t((apply(tstats$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T)))
	
	cond.T_IP.IC.inf<-ses.T_IP.IC<ses.T_IP.IC.inf
	cond.T_IC.IR.inf<-ses.T_IC.IR<ses.T_IC.IR.inf
	cond.T_PC.PR.inf<-ses.T_PC.PR<ses.T_PC.PR.inf
	
	cond.T_IP.IC.sup<-ses.T_IP.IC>ses.T_IP.IC.sup
	cond.T_IC.IR.sup<-ses.T_IC.IR>ses.T_IC.IR.sup
	cond.T_PC.PR.sup<-ses.T_PC.PR>ses.T_PC.PR.sup
	
	#________________________________________
	if(bysites==F){
		
		if(is.null(col.obj)) {col.obj<-rainbow(dim(tstats$T_IP.IC)[2])}
		else{}
		
		#__________
		#First panel of figures
		par(mfrow=c(sqrt(dim(tstats$T_IP.IC)[2])+1,sqrt(dim(tstats$T_IP.IC)[2])+1))
		par(mar=c(4,4,1,1))
		plot(0,0, xlim=c(-4,4), ylim=c(-4,4), cex.lab=1.2 ,ylab="ses.T_IP.IC", xlab="ses.T_IC.IR", type="n")
		abline(h=2)
		abline(v=2)
		abline(h=-2)
		abline(v=-2)
		text(0,0,"null \r\n model \r\n zone")		
		par(mar=c(1,2,1,1))
		
		for(t in 1:dim(tstats$T_IP.IC)[2]){
			plot(as.vector(ses.T_IP.IC)~as.vector(ses.T_IC.IR), col="grey", pch=20, main=rownames(ses.T_IC.IR)[t])
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_IC.IR.inf)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type="l")
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_IC.IR.sup)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type="l")
			points(sort(as.vector(ses.T_IC.IR)),as.vector(ses.T_IP.IC.inf)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type="l")
			points(sort(as.vector(ses.T_IC.IR)),as.vector(ses.T_IP.IC.sup)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type="l")
			segments(rep(rowMeans(ses.T_IC.IR, na.rm=T)[t], times=dim(tstats$T_IP.IC)[1]), rep(rowMeans(ses.T_IP.IC, na.rm=T)[t], times=dim(tstats$T_IP.IC)[1]) ,ses.T_IC.IR[t,],  ses.T_IP.IC[t,], col=col.obj[t])
		}
		
		#__________
		#Second panel of figures
		par(mfrow=c(sqrt(dim(tstats$T_IP.IC)[2])+1,sqrt(dim(tstats$T_IP.IC)[2])+1))
		par(mar=c(4,4,1,1))
		plot(0,0, xlim=c(-4,4), ylim=c(-4,4), cex.lab=1.2 ,ylab="ses.T_IP.IC", xlab="ses.T_PC.PR", type="n")
		abline(h=2)
		abline(v=2)
		abline(h=-2)
		abline(v=-2)
		text(0,0,"null \r\n model \r\n zone")		
		par(mar=c(1,2,1,1))
		
		for(t in 1:dim(tstats$T_IP.IC)[2]){
			plot(as.vector(ses.T_IP.IC)~as.vector(ses.T_PC.PR), col="grey", pch=20, main=rownames(ses.T_PC.PR)[t])
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_PC.PR.inf)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type="l")
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_PC.PR.sup)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type="l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IP.IC.inf)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type="l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IP.IC.sup)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type="l")
			segments(rep(rowMeans(ses.T_PC.PR, na.rm=T)[t], times=dim(tstats$T_IP.IC)[1]), rep(rowMeans(ses.T_IP.IC, na.rm=T)[t], times=dim(tstats$T_IP.IC)[1]) ,ses.T_PC.PR[t,],  ses.T_IP.IC[t,], col=col.obj[t])
		}
		
		#__________
		#Third panel of figures
		par(mfrow=c(sqrt(dim(tstats$T_IP.IC)[2])+1,sqrt(dim(tstats$T_IP.IC)[2])+1))
		par(mar=c(4,4,1,1))
		plot(0,0, xlim=c(-4,4), ylim=c(-4,4), cex.lab=1.2 ,ylab="ses.T_IC.IR", xlab="ses.T_PC.PR", type="n")
		abline(h=2)
		abline(v=2)
		abline(h=-2)
		abline(v=-2)
		text(0,0,"null \r\n model \r\n zone")	
		par(mar=c(1,2,1,1))
		
		for(t in 1:dim(tstats$T_IC.IR)[2]){
			plot(as.vector(ses.T_IC.IR)~as.vector(ses.T_PC.PR), col="grey", pch=20, main=rownames(ses.T_PC.PR)[t])
			points(sort(as.vector(ses.T_IC.IR))~as.vector(ses.T_PC.PR.inf)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type="l")
			points(sort(as.vector(ses.T_IC.IR))~as.vector(ses.T_PC.PR.sup)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type="l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IC.IR.inf)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type="l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IC.IR.sup)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type="l")
			segments(rep(rowMeans(ses.T_PC.PR, na.rm=T)[t], times=dim(tstats$T_IC.IR)[1]), rep(rowMeans(ses.T_IC.IR, na.rm=T)[t], times=dim(tstats$T_IC.IR)[1]) ,ses.T_PC.PR[t,],  ses.T_IC.IR[t,], col=col.obj[t])
		}
	}

	#________________________________________
	else if(bysites==T){
	
		if(is.null(col.obj)) {col.obj<-rainbow(dim(tstats$T_IP.IC)[1])}
		else{}
		
		#__________
		#First panel of figures
		par(mfrow=c(sqrt(dim(tstats$T_IP.IC)[1])+1,sqrt(dim(tstats$T_IP.IC)[1])+1))
		par(mar=c(4,4,1,1))
		plot(0,0, xlim=c(-4,4), ylim=c(-4,4), cex.lab=1.2 ,ylab="ses.T_IC.IR", xlab="ses.T_IC.IR", type="n")
		abline(h=2)
		abline(v=2)
		abline(h=-2)
		abline(v=-2)
		text(0,0,"null \r\n model \r\n zone")
		par(mar=c(1,2,1,1))
		
		for(s in 1:dim(tstats$T_IP.IC)[1]){
			plot(as.vector(ses.T_IP.IC)~as.vector(ses.T_IC.IR), col="grey", pch=20, main=colnames(ses.T_PC.PR)[s] )
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_IC.IR.inf)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type="l")
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_IC.IR.sup)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type="l")
			points(sort(as.vector(ses.T_IC.IR)),as.vector(ses.T_IP.IC.inf)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type="l")
			points(sort(as.vector(ses.T_IC.IR)),as.vector(ses.T_IP.IC.sup)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type="l")
			segments(rep(colMeans(ses.T_IC.IR, na.rm=T)[s], times=dim(tstats$T_IP.IC)[2]), rep(colMeans(ses.T_IP.IC, na.rm=T)[s], times=dim(tstats$T_IP.IC)[2]) ,ses.T_IC.IR[,s],  ses.T_IP.IC[,s], col=col.obj[s])
		}
		
		#__________
		#Second panel of figures
		par(mfrow=c(sqrt(dim(tstats$T_IP.IC)[1])+1,sqrt(dim(tstats$T_IP.IC)[1])+1))
		par(mar=c(4,4,1,1))
		
		plot(0,0, xlim=c(-4,4), ylim=c(-4,4), cex.lab=1.2 ,ylab="ses.T_IP.IC", xlab="ses.T_PC.PR", type="n")
		abline(h=2)
		abline(v=2)
		abline(h=-2)
		abline(v=-2)
		
		text(0,0,"null \r\n model \r\n zone")
				
		par(mar=c(1,2,1,1))
		
		for(s in 1:dim(tstats$T_IP.IC)[1]){
			plot(as.vector(ses.T_IP.IC)~as.vector(ses.T_PC.PR), col="grey", pch=20, main=colnames(ses.T_PC.PR)[s] )
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_PC.PR.inf)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type="l")
			points(sort(as.vector(ses.T_IP.IC))~as.vector(ses.T_PC.PR.sup)[order(as.vector(ses.T_IP.IC))][!is.na(as.vector(ses.T_IP.IC))], type="l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IP.IC.inf)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type="l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IP.IC.sup)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type="l")
			segments(rep(colMeans(ses.T_PC.PR, na.rm=T)[s], times=dim(tstats$T_IP.IC)[2]), rep(colMeans(ses.T_IP.IC, na.rm=T)[s], times=dim(tstats$T_IP.IC)[2]) ,ses.T_PC.PR[,s],  ses.T_IP.IC[,s], col=col.obj[s])
		}
		
		#__________
		#Third panel of figures
		par(mfrow=c(sqrt(dim(tstats$T_IP.IC)[1])+1,sqrt(dim(tstats$T_IP.IC)[1])+1))
		par(mar=c(4,4,1,1))
		
		plot(0,0, xlim=c(-4,4), ylim=c(-4,4), cex.lab=1.2 ,ylab="ses.T_IC.IR", xlab="ses.T_PC.PR", type="n")
		abline(h=2)
		abline(v=2)
		abline(h=-2)
		abline(v=-2)
		
		text(0,0,"null \r\n model \r\n zone")
				
		par(mar=c(1,2,1,1))

		for(s in 1:dim(tstats$T_IC.IR)[1]){
			plot(as.vector(ses.T_IC.IR)~as.vector(ses.T_PC.PR), col="grey", pch=20, main=colnames(ses.T_PC.PR)[s] )
			points(sort(as.vector(ses.T_IC.IR))~as.vector(ses.T_PC.PR.inf)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type="l")
			points(sort(as.vector(ses.T_IC.IR))~as.vector(ses.T_PC.PR.sup)[order(as.vector(ses.T_IC.IR))][!is.na(as.vector(ses.T_IC.IR))], type="l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IC.IR.inf)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type="l")
			points(sort(as.vector(ses.T_PC.PR)),as.vector(ses.T_IC.IR.sup)[order(as.vector(ses.T_PC.PR))][!is.na(as.vector(ses.T_PC.PR))], type="l")
			segments(rep(colMeans(ses.T_PC.PR, na.rm=T)[s], times=dim(tstats$T_IC.IR)[2]), rep(colMeans(ses.T_IC.IR, na.rm=T)[s], times=dim(tstats$T_IC.IR)[2]) ,ses.T_PC.PR[,s],  ses.T_IC.IR[,s], col=col.obj[s])
		}	
	}
	
	else{print("Error: obj need to be either traits or sites")}
	
	par(oldpar)
}


### Function to summarize traits and community which show a significative difference between observed and simulated value
summary.Tstats<-function(tstats=NULL, val.quant=c(0.025,0.975), type="all") {
	
	#________________________________________
	ses.T_IP.IC<-(tstats$T_IP.IC-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_IC.IR<-(tstats$T_IC.IR-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_PC.PR<-(tstats$T_PC.PR-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T))
	  
	ses.T_IP.IC.inf<-(apply(tstats$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_IC.IR.inf<-(apply(tstats$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_PC.PR.inf<-(apply(tstats$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T))
	  
	ses.T_IP.IC.sup<-(apply(tstats$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_IC.IR.sup<-(apply(tstats$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T))
	ses.T_PC.PR.sup<-(apply(tstats$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T))
	
	ses.T_IP.IC.mean<-t(colMeans((tstats$T_IP.IC-apply(tstats$T_IP.IC_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IP.IC_nm, c(3,2), function(x) sd(x, na.rm=T)), na.rm=T))
	ses.T_IC.IR.mean<-t(colMeans((tstats$T_IC.IR-apply(tstats$T_IC.IR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_IC.IR_nm, c(3,2), function(x) sd(x, na.rm=T)), na.rm=T))
	ses.T_PC.PR.mean<-t(colMeans((tstats$T_PC.PR-apply(tstats$T_PC.PR_nm, c(3,2), function(x) mean(x, na.rm=T)))/apply(tstats$T_PC.PR_nm, c(3,2), function(x) sd(x, na.rm=T)), na.rm=T))
	  
	ses.T_IP.IC.inf.mean<-apply(ses.T_IP.IC.inf,2, mean)
	ses.T_IC.IR.inf.mean<-apply(ses.T_IC.IR.inf,2, mean)
	ses.T_PC.PR.inf.mean<-apply(ses.T_PC.PR.inf,2, mean)	  
	
	ses.T_IP.IC.sup.mean<-apply(ses.T_IP.IC.sup,2, mean)
	ses.T_IC.IR.sup.mean<-apply(ses.T_IC.IR.sup,2, mean)
	ses.T_PC.PR.sup.mean<-apply(ses.T_PC.PR.sup,2, mean)
	
	#________________________________________
	#Condition to be significantly different from null models with respect to values of quantile choosen
	cond.T_IP.IC.inf<-ses.T_IP.IC<ses.T_IP.IC.inf
	cond.T_IC.IR.inf<-ses.T_IC.IR<ses.T_IC.IR.inf
	cond.T_PC.PR.inf<-ses.T_PC.PR<ses.T_PC.PR.inf
		
	cond.T_IP.IC.sup<-ses.T_IP.IC>ses.T_IP.IC.sup
	cond.T_IC.IR.sup<-ses.T_IC.IR>ses.T_IC.IR.sup
	cond.T_PC.PR.sup<-ses.T_PC.PR>ses.T_PC.PR.sup
	
	cond.T_IP.IC.inf.mean<-ses.T_IP.IC.mean<ses.T_IP.IC.inf.mean
	cond.T_IC.IR.inf.mean<-ses.T_IC.IR.mean<ses.T_IC.IR.inf.mean
	cond.T_PC.PR.inf.mean<-ses.T_PC.PR.mean<ses.T_PC.PR.inf.mean
		
	cond.T_IP.IC.sup.mean<-ses.T_IP.IC.mean>ses.T_IP.IC.sup.mean
	cond.T_IC.IR.sup.mean<-ses.T_IC.IR.mean>ses.T_IC.IR.sup.mean
	cond.T_PC.PR.sup.mean<-ses.T_PC.PR.mean>ses.T_PC.PR.sup.mean

	
	#________________________________________
	if(type=="binary"){
		summ.Tstats <-matrix("H0 not rejected",nrow=6, ncol=dim(cond.T_IP.IC.inf)[2])
		summ.Tstats <- rbind(cond.T_IP.IC.inf.mean, cond.T_IP.IC.sup.mean ,cond.T_IC.IR.inf.mean, cond.T_IC.IR.sup.mean ,cond.T_PC.PR.inf.mean, cond.T_IC.IR.sup.mean)
		rownames(summ.Tstats) <- c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		colnames(summ.Tstats) <- colnames(tstats$T_IP.IC)
	}
	
	#________________________________________
	else if(type=="percent"){
		
		summ.Tstats <-matrix("H0 not rejected",nrow=6, ncol=dim(cond.T_IP.IC.inf)[2])
		summ.Tstats[1,]<-paste(round(colSums(cond.T_IP.IC.inf, na.rm=T)/colSums(!is.na(cond.T_IP.IC.inf)),2)*100, "%", sep="")
		summ.Tstats[2,]<-paste(round(colSums(cond.T_IP.IC.sup, na.rm=T)/colSums(!is.na(cond.T_IP.IC.sup)),2)*100, "%", sep="")
		summ.Tstats[3,]<-paste(round(colSums(cond.T_IC.IR.inf, na.rm=T)/colSums(!is.na(cond.T_IC.IR.inf)),2)*100, "%", sep="")
		summ.Tstats[4,]<-paste(round(colSums(cond.T_IC.IR.sup, na.rm=T)/colSums(!is.na(cond.T_IC.IR.sup)),2)*100, "%", sep="")
		summ.Tstats[5,]<-paste(round(colSums(cond.T_PC.PR.inf, na.rm=T)/colSums(!is.na(cond.T_PC.PR.inf)),2)*100, "%", sep="")
		summ.Tstats[6,]<-paste(round(colSums(cond.T_PC.PR.sup, na.rm=T)/colSums(!is.na(cond.T_PC.PR.sup)),2)*100, "%", sep="")
				
		summ.Tstats[1,][cond.T_IP.IC.inf.mean]<-paste(round(colSums(cond.T_IP.IC.inf, na.rm=T)/colSums(!is.na(cond.T_IP.IC.inf)),2)[cond.T_IP.IC.inf.mean]*100, "%" ,"*", sep="")
		summ.Tstats[2,][cond.T_IP.IC.sup.mean]<-paste(round(colSums(cond.T_IP.IC.sup, na.rm=T)/colSums(!is.na(cond.T_IP.IC.sup)),2)[cond.T_IP.IC.sup.mean]*100, "%" ,"*", sep="")
		summ.Tstats[3,][cond.T_IC.IR.inf.mean]<-paste(round(colSums(cond.T_IC.IR.inf, na.rm=T)/colSums(!is.na(cond.T_IC.IR.inf)),2)[cond.T_IC.IR.inf.mean]*100, "%","*", sep="")
		summ.Tstats[4,][cond.T_IC.IR.sup.mean]<-paste(round(colSums(cond.T_IC.IR.sup, na.rm=T)/colSums(!is.na(cond.T_IC.IR.sup)),2)[cond.T_IC.IR.sup.mean]*100, "%","*", sep="")
		summ.Tstats[5,][cond.T_PC.PR.inf.mean]<-paste(round(colSums(cond.T_PC.PR.inf, na.rm=T)/colSums(!is.na(cond.T_PC.PR.inf)),2)[cond.T_PC.PR.inf.mean]*100, "%","*", sep="")
		summ.Tstats[6,][cond.T_PC.PR.sup.mean]<-paste(round(colSums(cond.T_PC.PR.sup, na.rm=T)/colSums(!is.na(cond.T_PC.PR.sup)),2)[cond.T_PC.PR.sup.mean]*100, "%","*", sep="")	
	
		rownames(summ.Tstats) <- c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		colnames(summ.Tstats) <- colnames(tstats$T_IP.IC)

	}
	
	#________________________________________
	else if(type=="site"){
	
		summ.Tstats <-matrix("H0 not rejected",nrow=6, ncol=dim(cond.T_IP.IC.inf)[2])
		for(t in 1: dim(cond.T_IP.IC.inf)[2]){
			
			if(sum(cond.T_IP.IC.inf[,t], na.rm=T)>0) 
				{summ.Tstats[1,t]<-paste( na.exclude(rownames(cond.T_IP.IC.inf)[cond.T_IP.IC.inf[,t]]), collapse=" ") }
			else{summ.Tstats[1,t]<-"H0 not rejected"}
				
			if(sum(cond.T_IP.IC.sup[,t], na.rm=T)>0) 
				{summ.Tstats[2,t]<-paste( na.exclude(rownames(cond.T_IP.IC.sup)[cond.T_IP.IC.sup[,t]]), collapse=" ")  }
			else{summ.Tstats[2,t]<-"H0 not rejected"}
								
			if(sum(cond.T_IC.IR.inf[,t], na.rm=T)>0)
				{summ.Tstats[3,t]<-paste( na.exclude(rownames(cond.T_IC.IR.inf)[cond.T_IP.IC.inf[,t]]), collapse=" ")  }
			else{summ.Tstats[3,t]<-"H0 not rejected"}
				
			if(sum(cond.T_IC.IR.sup[,t], na.rm=T)>0)
				{summ.Tstats[4,t]<-paste( na.exclude(rownames(cond.T_IC.IR.sup)[cond.T_IC.IR.sup[,t]]), collapse=" ")	}
			else{summ.Tstats[4,t]<-"H0 not rejected"}
				
			if(sum(cond.T_PC.PR.inf[,t], na.rm=T)>0)
				{summ.Tstats[5,t]<-paste( na.exclude(rownames(cond.T_PC.PR.inf)[cond.T_PC.PR.inf[,t]]), collapse=" ") 	}
			else{summ.Tstats[5,t]<-"H0 not rejected"}
				
			if(sum(cond.T_PC.PR.sup[,t], na.rm=T)>0)
				{summ.Tstats[6,t]<-paste( na.exclude(rownames(cond.T_PC.PR.sup)[cond.T_PC.PR.sup[,t]]), collapse=" ") 	}
			else{summ.Tstats[6,t]<-"H0 not rejected"}		
		}
		rownames(summ.Tstats) <- c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		colnames(summ.Tstats) <- colnames(tstats$T_IP.IC)

	}
	
	
	#________________________________________
	else if(type=="p.value"){
		summ.Tstats <-matrix("H0 not rejected",nrow=6, ncol=dim(cond.T_IP.IC.inf)[2])
		summ.Tstats <- rbind(tstats$pval_T_IP.IC.inf, tstats$pval_T_IP.IC.sup , tstats$pval_T_IC.IR.inf, tstats$pval_T_IC.IR.sup , tstats$pval_T_PC.PR.inf, tstats$pval_T_PC.PR.sup)
		rownames(summ.Tstats)<-c(paste(rep("T_IP.IC.inf",dim(tstats$T_IP.IC)[1]), rownames(tstats$T_IP.IC)), paste(rep("T_IP.IC.sup",dim(tstats$T_IP.IC)[1]), rownames(tstats$T_IP.IC)), paste(rep("T_IC.IR.inf",dim(tstats$T_IP.IC)[1]), rownames(tstats$T_IP.IC)), paste(rep("T_IC.IR.sup",dim(tstats$T_IP.IC)[1]), rownames(tstats$T_IP.IC)), paste(rep("T_PC.PR.inf",dim(tstats$T_IP.IC)[1]), rownames(tstats$T_IP.IC)), paste(rep("T_PC.PR.sup",dim(tstats$T_IP.IC)[1]), rownames(tstats$T_IP.IC)))
		colnames(summ.Tstats) <- colnames(tstats$T_IP.IC)
	}
	

	#________________________________________
	else if(type=="all"){
		summ.Tstats<-list()
		
		#__________
		##p.value
		summ.Tstats$p.value <-matrix("H0 not rejected", nrow=6, ncol=dim(cond.T_IP.IC.inf)[2])
		summ.Tstats$p.value <- rbind(tstats$pval_T_IP.IC.inf, tstats$pval_T_IP.IC.sup , tstats$pval_T_IC.IR.inf, tstats$pval_T_IC.IR.sup , tstats$pval_T_PC.PR.inf, tstats$pval_T_PC.PR.sup)
	
		#__________
		##percent
		summ.Tstats$percent <-matrix("H0 not rejected",nrow=6, ncol=dim(cond.T_IP.IC.inf)[2])
		summ.Tstats$percent[1,]<-paste(round(colSums(cond.T_IP.IC.inf, na.rm=T)/colSums(!is.na(cond.T_IP.IC.inf)),2)*100, "%", sep="")
		summ.Tstats$percent[2,]<-paste(round(colSums(cond.T_IP.IC.sup, na.rm=T)/colSums(!is.na(cond.T_IP.IC.sup)),2)*100, "%", sep="")
		summ.Tstats$percent[3,]<-paste(round(colSums(cond.T_IC.IR.inf, na.rm=T)/colSums(!is.na(cond.T_IC.IR.inf)),2)*100, "%", sep="")
		summ.Tstats$percent[4,]<-paste(round(colSums(cond.T_IC.IR.sup, na.rm=T)/colSums(!is.na(cond.T_IC.IR.sup)),2)*100, "%", sep="")
		summ.Tstats$percent[5,]<-paste(round(colSums(cond.T_PC.PR.inf, na.rm=T)/colSums(!is.na(cond.T_PC.PR.inf)),2)*100, "%", sep="")
		summ.Tstats$percent[6,]<-paste(round(colSums(cond.T_PC.PR.sup, na.rm=T)/colSums(!is.na(cond.T_PC.PR.sup)),2)*100, "%", sep="")
				
		summ.Tstats$percent[1,][cond.T_IP.IC.inf.mean]<-paste(round(colSums(cond.T_IP.IC.inf, na.rm=T)/colSums(!is.na(cond.T_IP.IC.inf)),2)[cond.T_IP.IC.inf.mean]*100, "%" ,"*", sep="")
		summ.Tstats$percent[2,][cond.T_IP.IC.sup.mean]<-paste(round(colSums(cond.T_IP.IC.sup, na.rm=T)/colSums(!is.na(cond.T_IP.IC.sup)),2)[cond.T_IP.IC.sup.mean]*100, "%" ,"*", sep="")
		summ.Tstats$percent[3,][cond.T_IC.IR.inf.mean]<-paste(round(colSums(cond.T_IC.IR.inf, na.rm=T)/colSums(!is.na(cond.T_IC.IR.inf)),2)[cond.T_IC.IR.inf.mean]*100, "%","*", sep="")
		summ.Tstats$percent[4,][cond.T_IC.IR.sup.mean]<-paste(round(colSums(cond.T_IC.IR.sup, na.rm=T)/colSums(!is.na(cond.T_IC.IR.sup)),2)[cond.T_IC.IR.sup.mean]*100, "%","*", sep="")
		summ.Tstats$percent[5,][cond.T_PC.PR.inf.mean]<-paste(round(colSums(cond.T_PC.PR.inf, na.rm=T)/colSums(!is.na(cond.T_PC.PR.inf)),2)[cond.T_PC.PR.inf.mean]*100, "%","*", sep="")
		summ.Tstats$percent[6,][cond.T_PC.PR.sup.mean]<-paste(round(colSums(cond.T_PC.PR.sup, na.rm=T)/colSums(!is.na(cond.T_PC.PR.sup)),2)[cond.T_PC.PR.sup.mean]*100, "%","*", sep="")	
		
		#__________
		##sites
		summ.Tstats$sites <-matrix("H0 not rejected",nrow=6, ncol=dim(cond.T_IP.IC.inf)[2])
		for(t in 1: dim(cond.T_IP.IC.inf)[2]){
			
			if(sum(cond.T_IP.IC.inf[,t], na.rm=T)>0) 
				{summ.Tstats$sites[1,t]<-paste( na.exclude(rownames(cond.T_IP.IC.inf)[cond.T_IP.IC.inf[,t]]), collapse=" ") }
			else{summ.Tstats$sites[1,t]<-"H0 not rejected"}
				
			if(sum(cond.T_IP.IC.sup[,t], na.rm=T)>0) 
				{summ.Tstats$sites[2,t]<-paste( na.exclude(rownames(cond.T_IP.IC.sup)[cond.T_IP.IC.sup[,t]]), collapse=" ")  }
			else{summ.Tstats$sites[2,t]<-"H0 not rejected"}
								
			if(sum(cond.T_IC.IR.inf[,t], na.rm=T)>0)
				{summ.Tstats$sites[3,t]<-paste( na.exclude(rownames(cond.T_IC.IR.inf)[cond.T_IP.IC.inf[,t]]), collapse=" ")  }
			else{summ.Tstats$sites[3,t]<-"H0 not rejected"}
				
			if(sum(cond.T_IC.IR.sup[,t], na.rm=T)>0)
				{summ.Tstats$sites[4,t]<-paste( na.exclude(rownames(cond.T_IC.IR.sup)[cond.T_IC.IR.sup[,t]]), collapse=" ")	}
			else{summ.Tstats$sites[4,t]<-"H0 not rejected"}
				
			if(sum(cond.T_PC.PR.inf[,t], na.rm=T)>0)
				{summ.Tstats$sites[5,t]<-paste( na.exclude(rownames(cond.T_PC.PR.inf)[cond.T_PC.PR.inf[,t]]), collapse=" ") 	}
			else{summ.Tstats$sites[5,t]<-"H0 not rejected"}
				
			if(sum(cond.T_PC.PR.sup[,t], na.rm=T)>0)
				{summ.Tstats$sites[6,t]<-paste( na.exclude(rownames(cond.T_PC.PR.sup)[cond.T_PC.PR.sup[,t]]), collapse=" ") 	}
			else{summ.Tstats$sites[6,t]<-"H0 not rejected"}		
		}
		
		#__________
		##binary
		summ.Tstats$binary <-matrix("H0 not rejected",nrow=6, ncol=dim(cond.T_IP.IC.inf)[2])
		summ.Tstats$binary <- rbind(cond.T_IP.IC.inf.mean, cond.T_IP.IC.sup.mean ,cond.T_IC.IR.inf.mean, cond.T_IC.IR.sup.mean ,cond.T_PC.PR.inf.mean, cond.T_IC.IR.sup.mean)
		
		#__________
		rownames(summ.Tstats$p.value)<-c(rep("T_IP.IC.inf",dim(tstats$T_IP.IC)[1]), rep("T_IP.IC.sup",dim(tstats$T_IP.IC)[1]), rep("T_IC.IR.inf",dim(tstats$T_IP.IC)[1]), rep("T_IC.IR.sup",dim(tstats$T_IP.IC)[1]), rep("T_PC.PR.inf",dim(tstats$T_IP.IC)[1]), rep("T_PC.PR.sup",dim(tstats$T_IP.IC)[1]))
		rownames(summ.Tstats$binary)<-c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		rownames(summ.Tstats$percent)<-c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		rownames(summ.Tstats$sites)<-c("T_IP.IC.inf", "T_IP.IC.sup", "T_IC.IR.inf", "T_IC.IR.sup", "T_PC.PR.inf", "T_PC.PR.sup")
		colnames(summ.Tstats$p.value) <- colnames(tstats$T_IP.IC)
		colnames(summ.Tstats$binary) <- colnames(tstats$T_IP.IC)
		colnames(summ.Tstats$sites) <- colnames(tstats$T_IP.IC)
		colnames(summ.Tstats$percent) <- colnames(tstats$T_IP.IC)
	}
	
	else{stop("Error: type must be 'binary', 'percent', 'p.value', 'site' or 'all'.")}
	
	return(summ.Tstats)
}



### Function to calcul SES
ses<-function(obs=NULL, nullmodel=NULL, val.quant=c(0.025,0.975) ){
	
	if(length(dim(obs))!=2 ) {
		obs<-as.matrix(obs)
	}
	
	if(dim(obs)[1]==dim(obs)[2]) {
		warnings("Observed matrix have the same number of rows and columns. The function is not able to detect automatically the correspondance between dimension of observed matrix and null model. You need to be sure that the null model is in the form of an array whithin the first and second dimension corresespond respectively to the first and second dimension of the observed matrix and the third dimension correspond to permutation")
		
		cond=c(1,2)
	}
	
	else{
	
		if (class(nullmodel)=="list"){
			if (class(nullmodel[[1]])=="list"){
				nullmodel<-array(unlist(nullmodel), dim=c(nrow(nullmodel[[1]][[1]]),ncol(nullmodel[[1]][[1]]), length(unlist(nullmodel))/nrow(nullmodel[[1]][[1]])/ncol(nullmodel[[1]][[1]])))
			}
			
			else {nullmodel<-array(unlist(nullmodel), dim=c(nrow(nullmodel[[1]]),ncol(nullmodel[[1]]), length(unlist(nullmodel))/nrow(nullmodel[[1]])/ncol(nullmodel[[1]])))}
		}
		
		if (class(obs)=="list"){
			obs<-matrix(obs[[1]], nrow=nrow(obs[[1]]), ncol=ncol(obs[[1]]))
		}
		
		if(!is.null(dim(obs))) {	
			cond<-c(NA,NA)
			
			if(dim(obs)[1]==dim(nullmodel)[1]){
			cond[1]<-1
			}
			
			if(dim(obs)[1]==dim(nullmodel)[2]){
			cond[1]<-2
			}
			
			if(length(dim(nullmodel))==3){
				if(dim(obs)[1]==dim(nullmodel)[3]){
				cond[1]<-3
				}
			}
			
			if(dim(obs)[2]==dim(nullmodel)[1]){
			cond[2]<-1
			}
			
			if(dim(obs)[2]==dim(nullmodel)[2]){
			cond[2]<-2
			}
			
			if(length(dim(nullmodel))==3){
				if(dim(obs)[2]==dim(nullmodel)[3]){
				cond[2]<-3
				}	
			}
		}
	}
	
	cond<-na.omit(cond)
	
	res<-list()
	res$ses<-(obs-apply(nullmodel, cond, function(x) mean(x, na.rm=T)))/apply(nullmodel, cond, function(x) sd(x, na.rm=T))
	res$ses.inf<-(apply(nullmodel, cond, function(x) quantile(x, na.rm=T, prob=val.quant[1]))-apply(nullmodel, cond, function(x) mean(x, na.rm=T)))/apply(nullmodel, cond, function(x) sd(x, na.rm=T))
	res$ses.sup<-(apply(nullmodel, cond, function(x) quantile(x, na.rm=T, prob=val.quant[2]))-apply(nullmodel, cond, function(x) mean(x, na.rm=T)))/apply(nullmodel, cond, function(x) sd(x, na.rm=T))

	return(res)
}



### Function to represent summarize Tstats
plot.bar.Tstats<-function(tstats=NULL, val.quant=c(0.025,0.975), col.Tstats=c("red","purple","green","white")){
    
  T_IP.IC.inf<-apply(tstats$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))
  T_IC.IR.inf<-apply(tstats$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))
  T_PC.PR.inf<-apply(tstats$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[1]))
  
  T_IP.IC.sup<-apply(tstats$T_IP.IC_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))
  T_IC.IR.sup<-apply(tstats$T_IC.IR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))
  T_PC.PR.sup<-apply(tstats$T_PC.PR_nm, c(3,2), function(x) quantile(x, na.rm=T, prob=val.quant[2]))
    
  df.bar<-barplot(rbind(colMeans(na.omit(tstats$T_IP.IC)), colMeans(na.omit(tstats$T_IC.IR)),colMeans(na.omit(tstats$T_PC.PR)),0), beside=T, plot=F)
  barplot(rbind(colMeans(na.omit(tstats$T_IP.IC)), colMeans(na.omit(tstats$T_IC.IR)),colMeans(na.omit(tstats$T_PC.PR)),0), col=col.Tstats, beside=T, ylim=c(min(c(T_IP.IC.inf,T_IC.IR.inf,T_PC.PR.inf,colMeans(na.omit(tstats$T_IP.IC))-apply(na.omit(tstats$T_IP.IC), 2,sd),colMeans(na.omit(tstats$T_IC.IR))-apply(na.omit(tstats$T_IC.IR), 2,sd),colMeans(na.omit(tstats$T_PC.PR))-apply(na.omit(tstats$T_PC.PR), 2,sd)), na.rm=T) , max(c(colMeans(na.omit(tstats$T_IP.IC))+apply(na.omit(tstats$T_IP.IC), 2,sd),colMeans(na.omit(tstats$T_IC.IR))+apply(na.omit(tstats$T_IC.IR), 2,sd),colMeans(na.omit(tstats$T_PC.PR))+apply(na.omit(tstats$T_PC.PR), 2,sd)), na.rm=T)))
  segments( df.bar[1,], colMeans(na.omit(tstats$T_IP.IC))+apply(na.omit(tstats$T_IP.IC), 2,sd),df.bar[1,],colMeans(na.omit(tstats$T_IP.IC))-apply(na.omit(tstats$T_IP.IC), 2,sd))
  segments( df.bar[2,], colMeans(na.omit(tstats$T_IC.IR))+apply(na.omit(tstats$T_IC.IR), 2,sd),df.bar[2,],colMeans(na.omit(tstats$T_IC.IR))-apply(na.omit(tstats$T_IC.IR), 2,sd))
  segments( df.bar[3,], colMeans(na.omit(tstats$T_PC.PR))+apply(na.omit(tstats$T_PC.PR), 2,sd),df.bar[3,],colMeans(na.omit(tstats$T_PC.PR))-apply(na.omit(tstats$T_PC.PR), 2,sd))
  
  points(type="l", df.bar[1,], colMeans(T_IP.IC.sup, na.rm=T), col=col.Tstats[1])
  points(type="l", df.bar[2,], colMeans(T_IC.IR.sup, na.rm=T), col=col.Tstats[2])
  points(type="l", df.bar[3,], colMeans(T_PC.PR.sup, na.rm=T), col=col.Tstats[3])
  
  points(type="l", df.bar[1,], colMeans(T_IP.IC.inf, na.rm=T), col=col.Tstats[1])
  points(type="l", df.bar[2,], colMeans(T_IC.IR.inf, na.rm=T), col=col.Tstats[2])
  points(type="l", df.bar[3,], colMeans(T_PC.PR.inf, na.rm=T), col=col.Tstats[3])
  
}


#~ if bysites=F, plot this function by traits
plot.filter<-function(index.list, color.cond=NULL, val.quant=c(0.025,0.975), cex.text =0.8, plot.ask=TRUE, srt.text=90, xlim=NULL, ylim=NULL, bysites=T){
		
	namesindex.all<-names(index.list)
	nindex<-length(names(index.list))/2
	namesindex<-names(index.list)[seq(1,nindex*2, by=2)]
	namestraits<-colnames(index.list[[1]])
	namesplots<-rownames(index.list[[1]]) 
	
	if(is.null(namesplots)) {print("rownames of index.list[[1]] is empty so names of plots cannot be plot")}
	if(is.null(namestraits)) {print("colnames of index.list[[1]] is empty so names of traits cannot be plot")}
	
	ncom<-dim(index.list[[1]])[1]
	ntr<-dim(index.list[[1]])[2]
	
	
	#________________________________________
	#Calcul of standardised effect size
		
	res<-list()
	for (i in seq(1,nindex*2, by=2)){
		res[[eval(namesindex.all[i])]] <- ses(obs=index.list[[i]], nullmodel=index.list[[i+1]], val.quant=val.quant)
	}

		
	if(is.null(ylim)){ ylim=c(0.5,nindex+0.5)}
	if(is.null(xlim)){ xlim=c(min(c(-2,unlist(res)), na.rm=T),max(c(2,unlist(res)), na.rm=T))}
	oldpar<-par(no.readonly = TRUE)
	par(ask=plot.ask)
	
	if(is.null(color.cond)) {color.cond=c("blue","orange")}
		
	if(bysites==T){
		for (s in 1: ncom){
			plot(mean(res[[eval(namesindex.all[i])]]$ses[s,], na.rm=T), (1:nindex)[i] ,bty="n", cex.lab=0.8, yaxt="n", xlab=paste("SES", namesplots[s]), ylim=ylim, xlim=xlim, pch=16, type="n")
			abline(v=0)	
					
			for(i in 1:nindex){
				abline(h=(1:nindex)[i], lty=2, col="lightgray")	
					
				segments(mean(res[[eval(namesindex[i])]]$ses.sup[s,], na.rm=T), (1:nindex)[i], mean(res[[eval(namesindex[i])]]$ses.inf[s,], na.rm=T), (1:nindex)[i])
				
				points(mean(res[[eval(namesindex[i])]]$ses.sup[s,], na.rm=T), (1:nindex)[i], pch="|")
				points(mean(res[[eval(namesindex[i])]]$ses.inf[s,], na.rm=T), (1:nindex)[i], pch="|")
				
				points(res[[eval(namesindex[i])]]$ses[s,], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[s,]) ), pch="*")
								
				cond.sup<-res[[eval(namesindex[i])]]$ses[s,]>res[[eval(namesindex[i])]]$ses.sup[s,]
				points(res[[eval(namesindex[i])]]$ses[s,][cond.sup], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[s,][cond.sup]) ), pch="*", cex=3, col=color.cond[2])
				
				cond.inf<-res[[eval(namesindex[i])]]$ses[s,]<res[[eval(namesindex[i])]]$ses.inf[s,]
				points(res[[eval(namesindex[i])]]$ses[s,][cond.inf], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[s,][cond.inf]) ), pch="*", cex=3, col=color.cond[1])
				
				
				points(mean(res[[eval(namesindex[i])]]$ses[s,], na.rm=T), (1:nindex)[i], col="red", pch=16)
				
				text(1, (1:nindex)[i]+0.3, namesindex[i], cex=0.8,  pos=4, font=2)
								
				chh <- par()$cxy[ 2 ]  ##  character height
				text(res[[eval(namesindex[i])]]$ses[s,], chh + rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[s,]) ), namestraits, cex=cex.text, srt=srt.text,)
				
			}
		}	
	}
	
	
	else if(bysites==F){
		for (t in 1: ntr){
			plot(mean(res[[eval(namesindex.all[i])]]$ses[,t], na.rm=T), (1:nindex)[i] ,bty="n", cex.lab=0.8, yaxt="n", xlab=paste("SES", namestraits[t]), ylim=ylim, xlim=xlim, pch=16, type="n")
			abline(v=0)	
					
			for(i in 1:nindex){
				
				abline(h=(1:nindex)[i], lty=2, col="lightgray")		
				
				segments(mean(res[[eval(namesindex[i])]]$ses.sup[,t], na.rm=T), (1:nindex)[i], mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm=T), (1:nindex)[i])
				
				points(mean(res[[eval(namesindex[i])]]$ses.sup[,t], na.rm=T), (1:nindex)[i], pch="|")
				points(mean(res[[eval(namesindex[i])]]$ses.inf[,t], na.rm=T), (1:nindex)[i], pch="|")
				
				points(res[[eval(namesindex[i])]]$ses[,t], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[,t]) ), pch="*")
								
				cond.sup<-res[[eval(namesindex[i])]]$ses[,t]>res[[eval(namesindex[i])]]$ses.sup[,t]
				points(res[[eval(namesindex[i])]]$ses[,t][cond.sup], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[,t][cond.sup]) ), pch="*", cex=3, col=color.cond[2])
				
				cond.inf<-res[[eval(namesindex[i])]]$ses[,t]<res[[eval(namesindex[i])]]$ses.inf[,t]
				points(res[[eval(namesindex[i])]]$ses[,t][cond.inf], rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[,t][cond.inf]) ), pch="*", cex=3, col=color.cond[1])
				
				
				points(mean(res[[eval(namesindex[i])]]$ses[,t], na.rm=T), (1:nindex)[i], col="red", pch=16)
				
				text(1, (1:nindex)[i]+0.3, namesindex[i], cex=0.8,  pos=4, font=2)
								
				chh <- par()$cxy[ 2 ]  ##  character height
				text(res[[eval(namesindex[i])]]$ses[,t], chh + rep( (1:nindex)[i], length(res[[eval(namesindex[i])]]$ses[,t]) ), namesplots, cex=cex.text, srt=srt.text,)
				
					
			}
		}	
	}
}




plot.filter.c<-function(com.object.ind=NULL, com.object.sp=NULL, tstats=NULL, color.cond=NULL, val.quant=c(0.025, 0.975), xlim=NULL, ylim=NULL, bysites=TRUE, cex.text =0.8, plot.ask=TRUE, srt.text=90){
  	
	
	index.list.internal.ind<-list(T_IP.IC=tstats$T_IP.IC, T_IP.IC_nm=tstats$T_IP.IC_nm, 
							  CV.ind=t(com.object.ind$CV.NND_obs), CV.ind.null=com.object.ind$Null_CV,
							  kurtosis.ind=t(com.object.ind$kurtosis_obs), kurtosis.ind.null=com.object.ind$Null_kurtosis)

	
	index.list.internal.sp<-list(CV.sp=t(com.object.sp$CV.NND_obs), CV.sp.null=com.object.sp$Null_CV,
							  kurtosis.sp=t(com.object.sp$kurtosis_obs), kurtosis.sp.null=com.object.sp$Null_kurtosis)

	rownames(tstats$T_IC.IR)<-rownames(tstats$T_IP.IC)
	index.list.external.ind<-list(T_IC.IR=tstats$T_IC.IR, T_IC.IR_nm=tstats$T_IC.IR_nm, 
								  range.ind=t(com.object.ind$Range_obs), range.ind.null=com.object.ind$Null_range,
								  mean.ind=t(com.object.ind$Mean_obs), mean.ind.null=com.object.ind$Null_mean)
	
	rownames(tstats$T_PC.PR)<-rownames(tstats$T_IP.IC)
	index.list.external.sp<-list(T_PC.PR=tstats$T_PC.PR, T_PC.PR_nm=tstats$T_PC.PR_nm,
								  range.sp=t(com.object.sp$Range_obs), range.sp.null=com.object.sp$Null_range,
								  mean.sp=t(com.object.sp$Mean_obs), mean.sp.null=com.object.sp$Null_mean)

    par(mfrow=c(2,2))
                 
	plot.filter(index.list.internal.ind, color.cond=color.cond, val.quant=val.quant, cex.text =cex.text, plot.ask=plot.ask, srt.text=srt.text, xlim=xlim, ylim=ylim, bysites=bysites)
	plot.filter(index.list.internal.sp, color.cond=color.cond, val.quant=val.quant, cex.text =cex.text, plot.ask=plot.ask, srt.text=srt.text, xlim=xlim, ylim=ylim, bysites=bysites)
	plot.filter(index.list.external.ind, color.cond=color.cond, val.quant=val.quant, cex.text =cex.text, plot.ask=plot.ask, srt.text=srt.text, xlim=xlim, ylim=ylim, bysites=bysites)
	plot.filter(index.list.external.sp, color.cond=color.cond, val.quant=val.quant, cex.text =cex.text, plot.ask=plot.ask, srt.text=srt.text, xlim=xlim, ylim=ylim, bysites=bysites)

	
    par(mfrow=c(1,1))
	
}

	

# plot ses of an index against an other variable wich correspond to plot. For example species richness or a gradient variable

plot.ses.var<-function(index.list, variable=NULL, color.traits=NULL, val.quant=c(0.025,0.975), resume=FALSE){

	y<-variable
	
	namesindex.all<-names(index.list)
	nindex<-length(names(index.list))/2
	namesindex<-names(index.list)[seq(1,nindex*2, by=2)]
	namestraits<-colnames(index.list[[1]])
	namescommunity<-rownames(index.list[[1]])
	
	ncom<-dim(index.list[[1]])[1]
	ntr<-dim(index.list[[1]])[2]
	
	if(is.null(color.traits)){
		color.traits<-palette(terrain.colors(ntr))
	}
	
	#________________________________________
	#Calcul of standardised effect size
	
	res<-list()
	for (i in seq(1,nindex*2, by=2)){
		res[[eval(namesindex.all[i])]] <- ses(obs=index.list[[i]], nullmodel=index.list[[i+1]], val.quant=val.quant)
	}

	old.par<-par()
	par(mfrow=c(ceiling(sqrt(nindex))-1, ceiling(sqrt(nindex))))
	
	ylim=c(min(y, na.rm=T), max(y, na.rm=T))
		
		
	for(i in seq(1,nindex*2, by=2)){
		if(resume==FALSE){xlim=c(min(c(-4,res[[eval(namesindex.all[i])]]$ses), na.rm=T), max(c(4,res[[eval(namesindex.all[i])]]$ses), na.rm=T))}
		else{xlim=c(min(c(-4,rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm=T)-apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm=T), rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm=T)+apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm=T)), na.rm=T), max(c(4,rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm=T)-apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm=T), rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm=T)+apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm=T)), na.rm=T))  }
		
		plot(0, 0 ,bty="n", cex.lab=0.8, xlab=paste("SES", namesindex.all[i]), ylim=ylim, xlim=xlim, pch=16, type="n")		
		abline(v=0, lty=1, col="black")	
		
		
		if(resume==TRUE){
			points(rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm=T), y, pch=16)
			segments(rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm=T)-apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm=T), y, rowMeans(res[[eval(namesindex.all[i])]]$ses, na.rm=T)+apply(res[[eval(namesindex.all[i])]]$ses, 1,sd,na.rm=T), y, pch=16)
		}
		
		else{
			for(nco in 1:ncom){
				abline(h=y[nco], lty=2, pch=0.5, col="lightgray")
			}
			
			for (t in 1: ntr){
				points(res[[eval(namesindex.all[i])]]$ses[,t] , y, pch=16, col=color.traits[t])		
			}
		}	
	}
	
	if(resume!=TRUE){
		plot(0, 0 ,bty="n", cex.lab=0.8, xlab=paste("SES", namesindex.all[i]), ylim=ylim, xlim=xlim, pch=16, type="n")
		legend("center", legend=namestraits, fill=color.traits, bty="n", ncol= round(sqrt(nlevels(as.factor(namestraits)))-1 ) )
	}
	
	par(old.par)
	
}


plot.dens<-function(traits=NULL, var.1=NULL, var.2=NULL, col.dens=NULL, plot.ask=TRUE, ylim.cex=1, cex.leg=0.8, polyg=TRUE, sort.by.traits=F)  {

	var.1<-as.factor(as.vector(var.1))
	var.2<-as.factor(as.vector(var.2))
	oldpar<-par(no.readonly = TRUE)
	par(ask=plot.ask)
	
	namestraits<-colnames(traits)
	namescommunity<-unique(var.1)
	ncom<-length(namescommunity)
	
	ntr<-dim(traits)[2]
	
	if(is.null(col.dens)) { col.dens<-rainbow(nlevels(as.factor(var.2))) }
	
	par(mfrow=c(2,2))
	
	if(!sort.by.traits){
		for(co in 1:ncom){
			for(t in 1:ntr){
				if(length(na.omit(traits[as.factor(var.1)==levels(as.factor(var.1))[co],t]))>1){
				
					plot(main= paste(namestraits[t],levels(as.factor(var.1))[co], " "), density(traits[as.factor(var.1)==levels(as.factor(var.1))[co],t], na.rm=T), ylim=c(0,max(density(traits[as.factor(var.1)==levels(as.factor(var.1))[co],t], na.rm=T)$y)*ylim.cex), xlim=c(min(traits[,t], na.rm=T),max(traits[,t], na.rm=T)), col="black")
					lines(density(traits[,t], na.rm=T), lty=2, col="grey")

					if (polyg==T) {
						x<-density(traits[as.factor(var.1)==levels(as.factor(var.1))[co],t], na.rm=T)$x
						y<-density(traits[as.factor(var.1)==levels(as.factor(var.1))[co],t], na.rm=T)$y
						polygon(c(x,rev(x)), c(rep(0,length(x)),rev(y)), border=NA, col=rgb(0.5,0.5,0.5,0.7))
					}
					
					for(s in 1:nlevels(as.factor(var.2))) {
						if(length(na.omit(traits[as.factor(var.1)==levels(as.factor(var.1))[co] & as.factor(var.2)==levels(as.factor(var.2))[s],t]))>1) 
						{lines(density(traits[as.factor(var.1)==levels(as.factor(var.1))[co] & as.factor(var.2)==levels(as.factor(var.2))[s],t], na.rm=T), col=col.dens[s])}
					
						else if(length(na.omit(traits[as.factor(var.1)==levels(as.factor(var.1))[co] & as.factor(var.2)==levels(as.factor(var.2))[s],t]))==1) 
						{points(0,na.omit(traits[as.factor(var.1)==levels(as.factor(var.1))[co] & as.factor(var.2)==levels(as.factor(var.2))[s],t]), col=col.dens[s])}
					}
				}
				
			}
		}
		par(mfrow=c(1,1))
		plot(0,0, type="n")
		legend("center", inset=0.05,levels(as.factor(var.2)), fill=col.dens, cex=cex.leg, bty="n", ncol= round(sqrt(nlevels(as.factor((var.2))))-1 ) )
	}	
	
	else if(sort.by.traits){
	}	
	par(oldpar)
}











