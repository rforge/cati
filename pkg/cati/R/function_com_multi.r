com.index.multi<-function(traits=NULL, index=NULL, by.factor=NULL, namesindex=NULL, nullmodels=NULL, ind.plot=NULL, sp=NULL, nperm=99, printprogress=TRUE){
	
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
	
	if(is.null(by.factor)) {  by.factor=rep(1,length(name_sp_sites)) }
	
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
		
			if(nullmodels[i]==1){nm<-array(unlist(traits.nm1),dim=c(ncol(traits), dim(traits)[1], nperm) )}
			else if(nullmodels[i]==2){nm<-array(unlist(traits.nm2),dim=c(ncol(traits), dim(traits)[1], nperm) )}
			else if(nullmodels[i]==3){nm<-array(unlist(traits.nm3),dim=c(ncol(traits), dim(traits)[1], nperm) )}
			else{print("nullmodels need 1, 2 or 3")}
			
			nm_n<-nm[,,n]
			colnames(nm_n)<-rownames(comm)
			rownames(nm_n)<-colnames(traits)
			
			functionindex= eval(index[i])
			
			dim2<-dim(by(t(nm_n), by.factor, function (x) eval(parse(text=functionindex))))[1]
			Null[[eval(namesindex[i])]] <- array(NA, dim=c(dim2, nperm) )
				
			if(is.null(dim2)) {
				Null[[eval(namesindex[i])]] <- array(NA, dim=c(1, nperm) )
			}
				
			for(n in 1:nperm){
				Null[[eval(namesindex[i])]][,n]  <- as.vector(by(t(nm[,,n]), by.factor, function (x) eval(parse(text=functionindex))))
			}
			
			if(printprogress==T){
				print(paste(eval(namesindex[i]), round(i/nindex*100,2),"%")) 
			} 
		}
		
	}
		  
	######################################## 
	####	Calcul of observed values	####
	######################################## 
	obs<-list()
	
	if(printprogress==T){print("calcul of observed values")}
	
	for(i in 1:nindex){
		if(nullmodels[i]==1){nm<-array(unlist(traits.nm1),dim=c(ncol(traits), dim(traits)[1], nperm) )}
		else if(nullmodels[i]==2){nm<-array(unlist(traits.nm2),dim=c(ncol(traits), dim(traits)[1], nperm) )}
		else if(nullmodels[i]==3){nm<-array(unlist(traits.nm3),dim=c(ncol(traits), dim(traits)[1], nperm) )}
		else{print("nullmodels need 1, 2 or 3")}
		
		nm_n<-nm[,,n]
		colnames(nm_n)<-rownames(comm)
		rownames(nm_n)<-colnames(traits)
		
		functionindex= eval(index[i])
		
		dim2<-dim(by(t(nm_n), by.factor, function (x) eval(parse(text=functionindex))))[1]
		obs[[eval(namesindex[i])]] <-rep(NA, times=dim2)
	
		if(nullmodels[i]==3) {
			traits.pop<-apply(traits, 2 , function (x) tapply(x, name_sp_sites, mean , na.rm=T))
			obs[[eval(namesindex[i])]] <- as.vector(by(t(traits.pop), by.factor, function (x) eval(parse(text=functionindex))))
		}
		
		else if(nullmodels[i]==1  |  nullmodels[i]==2) {
			obs[[eval(namesindex[i])]] <- as.vector(by(traits, by.factor, function (x) eval(parse(text=functionindex))))
		}
		
		obs[[eval(namesindex[i])]]<-as.vector(obs[[eval(namesindex[i])]])
		
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















plot.ses.multi<-function(index.list, col.index=c("red","purple","green"), type="normal", add.conf=TRUE, color.cond=TRUE, val.quant=c(0.025,0.975), grid.v=TRUE, grid.h=TRUE, xlim=NULL, ylim=NULL){
	#possible type = "simple",  "simple_range", "normal" and "barplot"	
	
	namesindex.all<-names(index.list)
	nindex<-length(names(index.list))/2
	namesindex<-names(index.list)[seq(1,nindex*2, by=2)]
	
	nfactor<-dim(as.matrix(index.list[[1]]))[1]
	
	if(length(col.index)<nindex){
		col.index<-palette(rainbow(nindex))
	}
	
	#________________________________________
	#Calcul of standardised effect size
	
	res<-list()
	for (i in seq(1,nindex*2, by=2)){
		res[[eval(namesindex.all[i])]] <- ses(obs=index.list[[i]], nullmodel=index.list[[i+1]], val.quant=val.quant)
	}
	
	if(is.null(ylim)) { ylim=c(0,5.5+(nindex+1))}
	if(is.null(xlim)) { xlim=c(min(c(unlist(res),-2), na.rm=T),max(c(unlist(res),2), na.rm=T))}
	
	par(mar=c(5, 7, 4, 2))
	plot(0, ylab="index",yaxt= "n", xlab="Standardized Effect Size", ylim=ylim, xlim=xlim, col="black", type="l")
	abline(v=0)
	
	if(grid.v==T) {
		range.<-max(c(unlist(res),2), na.rm=T)-min(c(unlist(res),-2), na.rm=T)
		
		vect.grid<-seq(min(unlist(res), na.rm=T),max(unlist(res), na.rm=T), by=round(range.,2)/9)
		for(j in vect.grid){
			abline(v=j, lty=2, col="lightgray")	
		}
	}
	
	if(grid.h==T) {
		for(j in seq(5.5,5.5+(nindex+1))  ){
			abline(h=j, lty=2, col="lightgray")
		}
	}
	
	#________________________________________
	#plot : possible type = "simple", "simple_range", "normal" and "barplot"	
	
	#__________
	if(type=="simple"){
			
		for(i in 1:nindex){
			
			if(color.cond==F){
				points(res[[eval(namesindex[i])]]$ses , rep(5.5+(nindex+1)-i, times=nfactor), pch=20, col=col.index[i])
			}
			
			if(color.cond==T){
				if(length(col.index)!=2*nindex) {col.index[(nindex+1):(nindex*2)]<-"grey"} 
				points(res[[eval(namesindex[i])]]$ses , rep(5.5+(nindex+1)-i, times=nfactor), pch=20, col=col.index[nindex+i])
				condition<-res[[eval(namesindex[i])]]$ses  > res[[eval(namesindex[i])]]$ses.sup   |   res[[eval(namesindex[i])]]$ses  < res[[eval(namesindex[i])]]$ses.inf 
				condition[is.na(condition)]<-FALSE					
				
				if(sum(condition)>0){
					points(res[[eval(namesindex[i])]]$ses [condition], rep(5.5+(nindex+1)-i, times=sum(condition)), pch=20, col=col.index[i])
				}
			}
							
			points(mean(res[[eval(namesindex[i])]]$ses ,na.rm=T), 5.5+(nindex+1)-i, pch=17, col=col.index[i])
											
			if(add.conf==T){
				points(mean(res[[eval(namesindex[i])]]$ses.sup ,na.rm=T), 5.5+(nindex+1)-i, pch="|", col=col.index[i])
				points(mean(res[[eval(namesindex[i])]]$ses.inf ,na.rm=T), 5.5+(nindex+1)-i, pch="|", col=col.index[i])
			}
	
			abline(seq(from=5.5, to=4.5+(nindex+1), by=nindex+1),b=0, lty=4, lwd=0.2)		
		}
	}
	
	#__________
	else if(type=="simple_range"){
		
		for(i in 1:nindex){

			if(color.cond==F){
				points(mean(res[[eval(namesindex[i])]]$ses ,na.rm=T), 5.5+(nindex+1)-i, pch=17, col=col.index[i])
			}
			
			if(color.cond==T){
				if(length(col.index)!=2*nindex) {col.index[(nindex+1):(nindex*2)]<-"grey"} 
				points(mean(res[[eval(namesindex[i])]]$ses ,na.rm=T), 5.5+(nindex+1)-i, pch=17, col=col.index[nindex+i])
				
				condition<-mean(res[[eval(namesindex[i])]]$ses , na.rm=T) > mean(res[[eval(namesindex[i])]]$ses.sup , na.rm=T)  |  mean( res[[eval(namesindex[i])]]$ses , na.rm=T) < mean(res[[eval(namesindex[i])]]$ses.inf , na.rm=T)
				if(sum(condition)>0){
					points(mean(res[[eval(namesindex[i])]]$ses ,na.rm=T), 5.5+(nindex+1)-i, pch=17, col=col.index[i])
				}
			}
			
			if(add.conf==T){
				points(mean(res[[eval(namesindex[i])]]$ses.sup ,na.rm=T), 5.5+(nindex+1)-i, pch="|", col=col.index[i])
				points(mean(res[[eval(namesindex[i])]]$ses.inf ,na.rm=T), 5.5+(nindex+1)-i, pch="|", col=col.index[i])
			}
			
			segments(mean(res[[eval(namesindex[i])]]$ses ,na.rm=T) + sd(res[[eval(namesindex[i])]]$ses ,na.rm=T), 5.5+(nindex+1)-i , mean(res[[eval(namesindex[i])]]$ses ,na.rm=T) - sd(res[[eval(namesindex[i])]]$ses ,na.rm=T), 5.5+(nindex+1)-i, col=col.index[i])
			
			points(min(res[[eval(namesindex[i])]]$ses ,na.rm=T), 5.5+(nindex+1)-i, pch="*", col=col.index[i])
			points(max(res[[eval(namesindex[i])]]$ses ,na.rm=T), 5.5+(nindex+1)-i, pch="*", col=col.index[i])
			
			abline(seq(from=5.5, to=4.5+(nindex+1), by=nindex+1),b=0, lty=4, lwd=0.2)		
		}
	}
	
	#__________
	else if(type=="normal"){
		
		
		for(i in 1:nindex){
				
			if(color.cond==F){
				points(res[[eval(namesindex[i])]]$ses , rep(5.5+(nindex+1)-i, times=nfactor), pch=20, col=col.index[i])
			}
			
			if(color.cond==T){
				if(length(col.index)!=2*nindex) {col.index[(nindex+1):(nindex*2)]<-"grey"} 
				points(res[[eval(namesindex[i])]]$ses , rep(5.5+(nindex+1)-i, times=nfactor), pch=20, col=col.index[nindex+i])
				condition<-res[[eval(namesindex[i])]]$ses  > res[[eval(namesindex[i])]]$ses.sup   |   res[[eval(namesindex[i])]]$ses  < res[[eval(namesindex[i])]]$ses.inf 
				condition[is.na(condition)]<-FALSE					
				
				if(sum(condition)>0){
					points(res[[eval(namesindex[i])]]$ses [condition], rep(5.5+(nindex+1)-i, times=sum(condition)), pch=20, col=col.index[i])
				}
			}
	
			if(add.conf==T){
				points(mean(res[[eval(namesindex[i])]]$ses.sup ,na.rm=T), 5.5+(nindex+1)-i, pch="|", col=col.index[i])
				points(mean(res[[eval(namesindex[i])]]$ses.inf ,na.rm=T), 5.5+(nindex+1)-i, pch="|", col=col.index[i])
			}
			
			segments(mean(res[[eval(namesindex[i])]]$ses ,na.rm=T) + sd(res[[eval(namesindex[i])]]$ses ,na.rm=T), 5.5+(nindex+1)-i , mean(res[[eval(namesindex[i])]]$ses ,na.rm=T) - sd(res[[eval(namesindex[i])]]$ses ,na.rm=T), 5.5+(nindex+1)-i, col=col.index[i])
			points(mean(res[[eval(namesindex[i])]]$ses ,na.rm=T), 5.5+(nindex+1)-i, pch=10, col=col.index[i])
			
			abline(seq(from=5.5, to=4.5+(nindex+1), by=nindex+1),b=0, lty=4, lwd=0.2)
	
			
		}
		
	}

	#__________
	else if(type=="barplot"){
		for(i in 1:nindex){
			
			segments(mean(res[[eval(namesindex[i])]]$ses ,na.rm=T), 5.5+(nindex+1)-i, 0, 5.5+(nindex+1)-i, pch=17, col=col.index[i], lwd=8)
			segments(mean(res[[eval(namesindex[i])]]$ses ,na.rm=T) + sd(res[[eval(namesindex[i])]]$ses ,na.rm=T), 5.5+(nindex+1)-i , mean(res[[eval(namesindex[i])]]$ses ,na.rm=T) - sd(res[[eval(namesindex[i])]]$ses ,na.rm=T), 5.5+(nindex+1)-i, col=col.index[i])
							
			if(add.conf==T){
				points(mean(res[[eval(namesindex[i])]]$ses.sup ,na.rm=T), 5.5+(nindex+1)-i, pch="|", col=col.index[i], cex=2)
				points(mean(res[[eval(namesindex[i])]]$ses.inf ,na.rm=T), 5.5+(nindex+1)-i, pch="|", col=col.index[i], cex=2)
			}
			abline(seq(from=5.5, to=4.5+(nindex+1), by=nindex+1),b=0, lty=4, lwd=0.2)
		}
		
	}
	else{print(paste("Error:",type,"is not a valid type of plot"))}
	
	legend("bottom", inset=.005, namesindex, fill=col.index, ncol=round(nindex/3+1) ,cex=0.6, bty="0")

	par(mar=c(5, 4, 4, 2) + 0.1) #return to default parameter
}


