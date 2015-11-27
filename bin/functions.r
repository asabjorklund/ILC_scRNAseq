# help functions for analysis of single cell RNAseq data from Innate lymphoid cells. 
# Written by Asa K Bjorklund 2015

cv2.var.genes <- function(R,ERCC,plot=FALSE,cutP=0.1,main="CV2 vs mean",ret.coord=FALSE){
    # function adopted from code in Brenneke et al. script 	
    # with rpkms instead of counts 
      
    library(statmod)
    nSamp<-ncol(ERCC)

    # calculate for ERCC
    rem<-which(rowSums(ERCC)==0)
    if (length(rem)>0) {
        ERCC<-ERCC[-rem,]
    }
    mE<-rowMeans(ERCC)
    vE<-apply(ERCC,1,var)
    cv2E<-vE/mE^2

   # for remaining genes.
   # remove all genes with sum<1
    remR<-which(rowSums(R)<1)
    R<-R[-remR,]
    mR<-rowMeans(R)
    vR<-apply(R,1,var)
    cv2R<-vR/mR^2

    # creating the fit
    mForFit<- unname( quantile( mE[ which( cv2E > .3 ) ], .95 ) )
    useForFit <- mE >= mForFit
    fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/mE[useForFit] ), cv2E[useForFit] )
    a0 <- unname( fit$coefficients["a0"] )
    a1 <- unname( fit$coefficients["a1tilde"])

    # finding genes in R with variation over ercc
    minBiolDisp <- .5^2
    cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
    psia1theta <- 1+a1
    testDenom <- ( mR * psia1theta + mR^2 * cv2th ) / ( 1 + cv2th/nSamp )

    p <- 1 - pchisq( vR * (nSamp-1) / testDenom, nSamp-1 )
    padj <- p.adjust( p, "BH" )
    sig <- padj < cutP
    sig[is.na(sig)] <- FALSE
    if (plot) {
         # Prepare the plot (scales, grid, labels, etc.)
       suppressWarnings(  plot( NULL, xaxt="n", yaxt="n",
             log="xy", xlim = range(mR), ylim = range(cv2R),
             xlab = "mean rpkm", ylab = "CV2" ))
        axis( 1, 10^(-4:5), c("0.0001", "0.001","0.01","0.1", "1", "10", "100", "1000",expression(10^4), expression(10^5) ) )
        axis( 2, 10^(-4:5), c("0.0001", "0.001","0.01","0.1", "1", "10", "100", "1000",expression(10^4), expression(10^5) ) )
        abline( h=10^(-3:4), v=10^(-5:5), col="#D0D0D0", lwd=2 )
        # Add the data points
        colS<-rep("red",nrow(R))
        colS[sig]<-"blue"
        points(mR,cv2R,pch=20,cex=0.2,col=colS)
        # add points for ERCC
        points(mE,cv2E,pch=20,cex=0.7,col="black")
        # Plot the fitted curve
        xg <- 10^seq( -2, 6, length.out=1000 )
        lines( xg, (a1)/xg + a0, col="#FF000080", lwd=3 )
        df <- nSamp - 1
        lines( xg, ( (a1)/xg + a0 ) * qchisq( .975, df ) / df,
              col="#FF000080", lwd=2, lty="dashed" )
        lines( xg, ( (a1)/xg + a0 ) * qchisq( .025, df ) / df,
              col="#FF000080", lwd=2, lty="dashed" )
        legend("topright",c("ERCC","non variable","variable"),col=c("black","red","blue"),pch=20)
        title(sprintf("%s\n%d var genes",main,length(which(sig))))
    }
    # convert the indices of sig to the original matrix before removing genes.
    idx.sig<-rep(FALSE,nrow(R)+length(remR))
    idx.sig[-remR]<-sig
    if (ret.coord) { return(list(idx.sig=idx.sig,coord=cbind(mR,cv2R))) }
    else {  return(idx.sig) }
}

change.duplicate.names<-function(names){
        # will give each gene one unique name, if some are duplicated
        t<-table(names)
        sel<-which(t>1)
        for (n in names(sel)){
            m<-which(names==n)
            for (i in 1:length(m)){
                names[m[i]] <- sprintf("%s..%d",names[m[i]],i)
            }
        }
        return(names)
}

remove.uniquenames<-function(names){
        sub("\\.\\.\\d$","",names)
}

# take a vector and create a list with indexes of all the names.
make_sets<-function(x){
   names<-names(table(x))
   l<-list()
   for (n in names){
       l[[n]]<-which(x==n)
   }
   return(l)
}

make_colors<-function(names,legend,colors){
        # make a color vector for each group based on sample names, a list of legend and colors
        col<-mat.or.vec(1,length(names))
        for (i in 1:length(legend)){
            col[names==legend[i]]<-colors[i]
            }
            return(col)
}


make.vioplot<-function(l,name="gene",text.col="black",ylim=NULL,horizontal=FALSE){
    library(vioplot)
    lm<-unlist(lapply(l,mean))
    col.idx<-round(log2(lm+1)*10)+1
    col.idx[col.idx>length(crange.log)]<-length(crange.log)
    cols<-col.log[col.idx]
    l<-lapply(l, function(x) log2(x+1))
    lz<- which(unlist(lapply(l,sum))==0)
    # if one gene is non-expressed, need to add pseudo values around zero
    for (z in lz)  { l[[z]]<-rnorm(length(l[[z]]),mean=0,sd=1e-10) }
    vioplot.list.small(l,col=cols,drawRect=F,border=NA,ylim=ylim,horizontal=horizontal)
    text(-2.5,1,name,srt=0,adj=c(0,0.2),col=text.col)
}

plot.color.bar<-function(col=col.log,crange=crange.log,main="log2(rpkm+1)"){
  barplot(rep(1,length(col)),col=col,border=NA,main=main,space=0,axes=F,cex.main=0.6)
  # round to closest 1,10,100,1000 etc.
  by.exp<-floor(log10(max(na.omit(crange))-min(na.omit(crange))))
  if (by.exp==1){
    by = 1
  }else {
    by<-10^(by.exp-1)
  }
  # if crange is not at even numbers.
  if (max(na.omit(crange)) > 100) {
    crange<-round(crange/10)*10
  }
  rem<-crange  %% by
  evens<-which(rem==0)
  if (length(evens) > 15) {
    rem<-crange  %% (by*5)
    evens<-which(rem==0)
  }
  lab<-crange[evens]
  axis(1,at=evens,labels=lab,cex.axis=0.5,las=2)
}


# also remove axes.
vioplot.list.small<- function (datas, range = 1.5, h = NULL, ylim = NULL, names = NULL,
    # modification of viplot function from cran, to take a list as input and plot smaller violins with no axes, means etc.		     
    horizontal = FALSE, col = "magenta", border = "black", lty = 1,
    lwd = 1, rectCol = "black", colMed = "white", pchMed = 19,
    at, add = FALSE, wex = 1, drawRect = TRUE,las=2,ylab=NULL)
{
    n <- length(datas)
    if (missing(at))
        at <- 1:n
    upper <- vector(mode = "numeric", length = n)
    lower <- vector(mode = "numeric", length = n)
    q1 <- vector(mode = "numeric", length = n)
    q3 <- vector(mode = "numeric", length = n)
    med <- vector(mode = "numeric", length = n)
    base <- vector(mode = "list", length = n)
    height <- vector(mode = "list", length = n)
    baserange <- c(Inf, -Inf)
    args <- list(display = "none")
    if (!(is.null(h)))
        args <- c(args, h = h)
    for (i in 1:n) {
        data <- datas[[i]]
        data.min <- min(data)
        data.max <- max(data)
        q1[i] <- quantile(data, 0.25)
        q3[i] <- quantile(data, 0.75)
        med[i] <- median(data)
        iqd <- q3[i] - q1[i]
        upper[i] <- min(q3[i] + range * iqd, data.max)
        lower[i] <- max(q1[i] - range * iqd, data.min)
        est.xlim <- c(min(lower[i], data.min), max(upper[i],
            data.max))
        smout <- do.call("sm.density", c(list(data, xlim = est.xlim),
            args))
        hscale <- 0.4/max(smout$estimate) * wex
       base[[i]] <- smout$eval.points
        height[[i]] <- smout$estimate * hscale
        t <- range(base[[i]])
        baserange[1] <- min(baserange[1], t[1])
        baserange[2] <- max(baserange[2], t[2])
    }
    if (!add) {
        xlim <- if (n == 1)
            at + c(-0.5, 0.5)
        else range(at) + min(diff(at))/2 * c(-1, 1)
        if (is.null(ylim)) {
            ylim <- baserange
        }
    }
    if (is.null(names)) {
        label <- 1:n
    }
    else {
        label <- names
    }
    boxwidth <- 0.05 * wex
    if (!add)
        plot.new()
    if (!horizontal) {
        if (!add) {
            plot.window(xlim = xlim, ylim = ylim)
            #axis(2)
            axis(2, labels = FALSE,tick=FALSE,col="white")
            #axis(1)
            axis(1, at = at, labels = FALSE,tick=FALSE,col="white")
        }
        box(col="white")
        if(length(col)!=n){ col<-rep(col[1],n) }
        for (i in 1:n) {
            polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])),
                c(base[[i]], rev(base[[i]])), col = col[i], border = border,
                lty = lty, lwd = lwd)
            if (drawRect) {
                lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd,
                  lty = lty)
                rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2,
                  q3[i], col = rectCol)
                points(at[i], med[i], pch = pchMed, col = colMed)
            }
        }
    }
    else {
        if (!add) {
            plot.window(xlim = ylim, ylim = xlim)
            axis(1, labels = FALSE,tick=FALSE,col="white")
            axis(2, at = at, labels = FALSE,tick=FALSE,col="white")
        }
        box(col="white")
        if(length(col)!=n){ col<-rep(col[1],n) }
        for (i in 1:n) {
            polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]],
                rev(at[i] + height[[i]])), col = col[i], border = border,
                lty = lty, lwd = lwd)
            if (drawRect) {
                lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd,
                  lty = lty)
                rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] +
                  boxwidth/2, col = rectCol)
                points(med[i], at[i], pch = pchMed, col = colMed)
            }
        }
    }
    invisible(list(upper = upper, lower = lower, median = med,
        q1 = q1, q3 = q3))
}

convert.to.color <- function(x,colscale,col.range=NULL){
  x.range<-range(na.omit(x))
  by=0.1
  if (is.null(col.range)){
    by=10
    col.range<-seq(x.range[1],x.range[2],by=by)
  }
  col.def<-colscale(length(col.range))
  col.idx<-round((x-x.range[1])/by)+1
  col.idx[col.idx>length(col.range)]<-length(col.range)
  cols<-col.def[col.idx]
  return(list(cols=cols,col.def=col.def,col.range=col.range))
}


##########################################
# functions for PCA plotting
run.pca<-function(data,seln=0,log.vals=TRUE,samples.col=TRUE,center=TRUE){
        if (!samples.col){ data<-t(data) }

        # remove zero read genes
        z<-which(rowSums(data)==0)
        if (length(z>0)){
            data<-data[-z,]
        }

        if (log.vals) { data <- log2(data+0.1) }

        # select top varied genes
        if (seln>0){
            sD<-apply(data,1,sd)
            o<-order(sD,decreasing=T)
            data<-data[o[1:seln],]
        }


        myPca <- prcomp(t(data),center=center,scale.=FALSE)
        vars<- myPca$sdev^2
        vars<- vars/sum(vars)
        pcnames<-mat.or.vec(1,length(vars))
        for (i in 1:length(vars)) {
            pcnames[i]<-sprintf("PC%d %.5f",i,vars[i])
        }

        myPca$pc.names<-as.vector(pcnames)
        return(myPca)
}

pca.plot <- function (data, log.vals=TRUE,seln=0,selpc=1:2,center=TRUE,samples.col=TRUE, ...){
         # data is either a prcomp object or a data matrix
         if (class(data) != "prcomp"){
            data<-run.pca(data,seln=0,log.vals=log.vals,samples.col=samples.col,center=center)
         }
         tmpPca <- as.data.frame(data$x[,selpc])
         colnames(tmpPca)<-data$pc.names[selpc]
         plot(tmpPca, ...)
         invisible(data)
}

##########################################
# functions for protein/RNA plotting

# normalize each protein and rna within the plates
# set all values below zero to zero
norm_distr <- function(x,cut=0){
    x[x<=cut]<-0
    # if none is expressed return all zeros
    if (sum(x)==0){
      return(rep(0,length(x)))
    }else {
      n = (x-min(x))/(max(x)-min(x))
      return(n)
    }
}

plot_mean<-function(data,name){
  mP<-colMeans(data[,prot.idx])
  mR<-colMeans(data[,rna.idx])

  cP<-cor(mP,mR)
  cS<-cor(mP,mR,method="spearman")

  plot(mP,mR, xlab="Mean protein abundance",ylab="Mean RNA abuncance",main=sprintf("Mean, %s, Pearson=%.4f, Spearman=%.4f",name,cP,cS))
  offset<-(max(mR)-min(mR))/30
  text(mP,mR+offset,markers,cex=0.7)
}

calculate_simstats<-function(data1,data2,ncut,nSim){
        output<-mat.or.vec(nSim,2)
        nS<-nrow(data1)
        for (s in 1:nSim) {
            set<-sample(1:nS)[1:ncut]
            output[s,1]<-cor(colMeans(data1[set,]), colMeans(data2[set,]))
            output[s,2]<-cor(colMeans(data1[set,]), colMeans(data2[set,]),method="spearman")
        }
        m<-colMeans(output)
        s<-apply(output,2,sd)
        return(c(m,s))
}


run_simulation_plot<-function(data,nSim,noPar=FALSE){
  nM<-ncol(data)/2
  c1P<-apply(data, 1, function(x) cor(x[1:nM],x[(nM+1):(nM*2)]))
  c1S<-apply(data, 1, function(x) cor(x[1:nM],x[(nM+1):(nM*2)],method="spearman"))
  stat1<-c(mean(c1P),mean(c1S),sd(c1P),sd(c1S))

  data.markers<-data[,1:nM]
  data.genes<-data[,(1+nM):(nM*2)]

  cAllP<-cor(colMeans(data.markers), colMeans(data.genes))
  cAllS<-cor(colMeans(data.markers), colMeans(data.genes),method="spearman")
  statAll<-c(cAllP, cAllS,0,0)
  nS<-nrow(data)
  by<-5
  if (nS>150) { by <- 10 }
  if (nS>300) { by <- 50 }

  subsets<-c(10,25,seq(50,nS-1 ,by=by))
  nSet<-length(subsets)

  simdata<-mat.or.vec(nSet,4) #mean/sd with P or S,
  for (i in 1:nSet) {
    simdata[i,]<-calculate_simstats(data.markers,data.genes,subsets[i],nSim)
  }

  x<-c(1,subsets,nS)
  simdata2<-rbind(stat1,simdata,statAll)

  library(gplots)

  if (!noPar) { par(mfrow=c(1,2)) }
  plotCI(x,simdata2[,1],simdata2[,3],col="black",main="Pearson",ylab="Pearson correlation",xlab="Number of cells",ylim=c(-0.1,1),xaxt="n")
  lines(x,simdata2[,1],col="black")
  axis(1,at=x,las=2,cex.axis=0.6)

  plotCI(x,simdata2[,2],simdata2[,4],col="black",main="Spearman",ylab="Spearman correlation",xlab="Number of cells",ylim=c(-0.1,1),xaxt="n")
  lines(x,simdata2[,2],col="black")
  axis(1,at=x,las=2,cex.axis=0.6)
}

