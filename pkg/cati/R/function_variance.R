#### Partition of variance by nested anova (Messier et al, 2010)####
# traits is the traits matrix and factors are the nested factors to take into account in the partition of variance
partvar<-function(traits, factors, printprogress=TRUE){
	
	traits<-as.matrix(traits)
	factors<-as.matrix(factors)
	nfactors <- ncol(factors)
	ntraits  <- ncol(traits)
	res<-matrix(0, nrow=nfactors+1, ncol=ntraits)
	colnames(res)<-colnames(traits)
	
	if(!is.null(colnames(factors)))
		{rownames(res)<-c(colnames(factors), "whithin")} 
	else {rownames(res)<-c(paste("factor",1:(nfactors),sep=""),  "whithin") ; colnames(factors)<-c(paste("factor", 1:(nfactors),sep="")) } 
	
	attach(as.data.frame(factors))
	
	for (t in 1 : ntraits) {
		trait<-traits[,t]
		functionlme= paste('varcomp(lme(trait~1, random=~1|', paste(colnames(factors), collapse='/'), ",na.action=na.omit),1)", sep="")
		res[,t]<-as.vector(eval(parse(text=functionlme)))
		
		if(printprogress==TRUE) 
			{print(paste(round(t/ntraits*100, 2), "%", sep=" ")) }
		else{}
	}
	
	attach(as.data.frame(factors))
	
	print(res)
}

#Pie of variance partitionning
pie.partvar<-function(var.part, col.pie=NA){
  nfactors <- nrow(var.part)
  ntraits  <- ncol(var.part)
  old.par<-par()
  par(mfrow=c(round(sqrt(ntraits))+1, round(sqrt(ntraits))))
  for (t in 1 : ntraits) {
    pie(var.part[,t], main= colnames(var.part)[t], col=col.pie , labels=rownames(var.part))
  }
  par(old.par)
}

#barplot of variance partitionning
bar.partvar<-function(var.part,  col.bar=NA, leg=FALSE){
	
	old.par<-par()
	par(mar=c(5,6.5,4,2), cex=0.7)
	
	if(leg==FALSE)
		{barplot(var.part, col=col.bar, las=1, horiz=T, xlab="% of variance")}
	else if(leg==TRUE) 
		{barplot(var.part, col=col.bar, las=1, horiz=T, legend.text=rownames(var.part), args.legend= list(cex=0.6), xlab="% of variance")}
	
	else{stop("leg must be TRUE or FALSE")}
	 par(old.par)
}
