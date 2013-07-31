counterPerChr <-
    function(seqname,fl,gnModel,mapq=10,lib.strand=c("none","sense","anti"),...){
        
        lib.strand <- match.arg(lib.strand)
        
        s.i <- seqinfo(BamFile(fl)) 
        seq.length <- seqlengths(s.i)[seqname]
        param <- ScanBamParam(which=GRanges(seqname, IRanges(1, seq.length)),
                              what='mapq',
                              tag='NH')

        ## Read the alignment file
        aln <- readGappedAlignments(fl,param=param)
        
        ## What type of RNA-Seq library are we dealing with?
        strand(aln) <- switch(lib.strand,
                              none = "*",
                              sense = strand(aln),
                              anti  = ifelse(strand(aln) == "+","-","+")
                              )
        ## Removing aligned reads with more than one aligments
        aln <- aln[!values(aln)$NH > 1]
        ## Removing the low map quality alignmentsw
        aln <- aln[!values(aln)$mapq <= mapq]
        ## Subseting our gnModel to only models on the curent chromsome
        if (class(gnModel) == "GRangesList"){
            gnModel <- gnModel[seqnames(unlist(range(gnModel))) == seqname]
        } else if (class(gnModel) == "GRanges"){
            gnModel <- gnModel[seqnames(gnModel) == seqname]
        } else {
            stop("Wrong gene model object")
        }

        ## Do not use aligments overlapping more than one genes
        aln <- aln[countOverlaps(aln, gnModel)==1 ]
        
        counts <- countOverlaps(gnModel, aln)
        names(counts) <- names(gnModel)
        return(counts)
    }


BamsGeneCount <- function(BFL,
                          gnModel,
                          lib.strand=c("none","sense","anti"),
                          nCores=16,
                          ...){
    lib.strand <- match.arg(lib.strand)
    ## Making sure there is an index
    parIndexBam(BFL,nCores)

    chrs <- sapply(BFL,function(BFL) seqnames(seqinfo(BFL)))
    files <- rep(names(BFL),each=length(seqnames(seqinfo(BFL[[1]]))))
    ## Reading the bam files, in parrels one chromosome at a time
    counts.raw <- mcmapply(counterPerChr
                           ,chrs
                           ,files
                           ,MoreArgs=list(gnModel=gnModel,lib.strand=lib.strand)
                           ,SIMPLIFY=FALSE
                           ,mc.cores=nCores
                           ,mc.preschedule=FALSE
                           )
    
    ## Reorganizing the lists of counts per chromosome into a list of counts per file
    counts <- lapply(split(counts.raw,f=files),function(cnts){
        to.ret <- do.call(c,cnts)
        names(to.ret) <- unlist(sapply(cnts,names))
        to.ret
    })
    
    names(counts) <- sub("\\.bam","",basename(names(BFL)))
    return(counts)
}

getReport <- function(DGE,ID2gene,FDR=0.01,logFC=log2(2)){
    require(edgeR)
    require(GenomicFeatures)
    
    DGE$isUp <- DGE$table$logFC >= logFC & DGE$table$FDR <= FDR
    DGE$isDown <- DGE$table$logFC <= -logFC & DGE$table$FDR <= FDR
        
    gene.symb <- ID2gene$symbol[match(rownames(DGE$table),ID2gene$id)]
    
    cntl.s <- DGE$design[,1]-rowSums(DGE$design) == 0
    exp.s <- DGE$samples$group == DGE$comparison
    
    read.counts <- round(DGE$fitted.values[,cntl.s|exp.s])
    DGE.stats <- t(apply(DGE$table,1,sprintf,fmt=c(rep("%0.3f",3),rep("%0.2E",2))))

    colnames(DGE.stats) <- colnames(DGE$table)
    
    fmt <- "[[http://flybase.org/reports/%s.html][%s]]"
    FB.links <- sprintf(fmt,rownames(DGE$table),gene.symb)

    RG <- unlist(range(gnModel))[rownames(DGE$table)]
    RG <- resize(RG,round(width(RG) * 1.1),fix='center')

    IGV.fmt <- '[[http://localhost:60151/goto?locus=%s:%s-%s][%s]]'
    IGVlinks <- sprintf(IGV.fmt,seqnames(RG),start(RG),end(RG),'toIGV')
    
    DGE$reports <- data.frame(symbol=FB.links,
                              read.counts,
                              DGE.stats,
                              "IGV links"=IGVlinks,
                              stringsAsFactors=FALSE)
    
    DGE$counts <- data.frame(up=sum(DGE$isUp),down=sum(DGE$isDown))

    return(DGE)
  }

uxonify <- function(txdb){ reduce(exonsBy(txdb,'gene')) }

ixonify <- function(txdb){
  uxons <- uxonfiy(txdb)
  genes <- unlist(range(uxons))
  psetdiff(genes,uxons)
}

parIndexBam <- function(BFL,cpu=12,...){
    require(Rsamtools)  
    require(parallel)
    
    BFL <- BFL[!file.exists(paste(path(BFL),".bai",sep=''))]
    if(length(BFL) == 0) return(NULL)
    
    ## No reason to be wastefull on ressources
    nClus <- ifelse(length(BFL) < cpu,length(BFL),cpu)
    
    ## Make sure the bam files are indexed
    ##mclapply(BFL,indexBam,mc.core=nClus,mc.preschedule=FALSE)
    mclapply(BFL,function(BF){
        system(paste('samtools index',path(BF),sep=" "))
    },mc.cores=nClus,mc.preschedule=FALSE)
    
    return(NULL)
}

topHatReport <- function(topHatDir,cpu=12,debug=FALSE){
  require(Rsamtools)  
  require(parallel)

  topHat.mapped <- BamFileList(list.files(path=topHatDir,pattern="^accepted_hits.bam$",recu=TRUE,full.names=TRUE))
  topHat.unmapped <- BamFileList(file.path(dirname(path(topHat.mapped)),'unmapped.bam'))

  stopifnot(length(topHat.mapped) == length(topHat.unmapped))
    
  if(debug){topHat.mapped <- topHat.mapped[1]
            topHat.unmapped <- topHat.unmapped[1]}
            
  parIndexBam(c(topHat.mapped,topHat.unmapped))
  
  toProcess <-
    do.call(rbind,lapply(topHat.mapped,function(file){data.frame(chr=seqnames(seqinfo(file)),
                                                                 file=path(file),
                                                                 stringsAsFactors=FALSE)
                                                    }))
  ## No reason to be wastefull on ressources
  nClus <- ifelse(nrow(toProcess) < cpu,nrow(toProcess),cpu)
  
  counter <-function(d){
    fl <- d$file
    seqname <- d$chr
    
    s.i <- seqinfo(BamFile(fl)) 
    seq.length <- seqlengths(s.i)[seqname]
    param <- ScanBamParam(which=GRanges(seqname, IRanges(1, seq.length)),
                          flag = scanBamFlag(isUnmappedQuery=NA)
                          ,tag='NH')
    
    ## Read the alignment file
    aln <- readGappedAlignments(fl,param=param)
    
    tab.cnt <- table(values(aln)$NH)
    
    ## normalize the counts per their number of hits
    tab.cnt <- tab.cnt/as.numeric(names(tab.cnt))
    
    counts <- c(unique=tab.cnt[as.numeric(names(tab.cnt))==1]/1e6,
                multi.match=sum(tab.cnt[as.numeric(names(tab.cnt))>1])/1e6)
  }
  ## map the counts per chromosome per file
  mapped.data <- mclapply(split(toProcess,1:nrow(toProcess)),
                          counter,
                          mc.cores=nClus,
                          mc.preschedule=FALSE
                          )
  
  ## Reduce to 1) chromosomes 2) files
  reduced.data <- lapply(split(mapped.data,toProcess$file),function(perFile){
    d <- do.call(rbind,perFile)
    colSums(d)
  })
  
  unmapped.reads <- mclapply(topHat.unmapped,function(BF) countBam(BF)$records/1e6,
                             mc.cores=nClus,
                             mc.preschedule=FALSE
                             )
  counts <- data.frame(t(mapply(c,reduced.data,unmapped.reads)))
  counts <- t(apply(counts,1,function(cnt)
                    c(sprintf("%0.1fM",c(sum(cnt),cnt)),sprintf("%0.1f%%",cnt/sum(cnt)*100))
                    ))
  
  counts <- cbind(basename(dirname(path(topHat.mapped))),counts)
  colnames(counts) <- c("Sample",
                        "total reads",
                        "singly mapped",
                        "multiply mapped",
                        "unmapped",
                        "singly mapped (%)",
                        "multiply mapped (%)",
                        "unmapped reads (%)")
  counts
}


reportFracSense <-
    function(BFL,txdbfile,min.lim=10,lib.strand=c("anti","sense"),CPU=16,...){
        lib.strand <- match.arg(lib.strand)
        txdb<-loadDb(txdbfile)
        genes <- unlist(range(exonsBy(txdb,"gene")))

        genes <- genes[seqnames(genes) %in% seqnames(seqinfo(BFL[[1]]))]
        
        fracSense <-function(seqname,BF,gnModel){
            s.i <- seqinfo(BF)

            seq.length <- seqlengths(s.i)[seqname]
            param <- ScanBamParam(which=GRanges(seqname, IRanges(1, seq.length)),
                                  tag='NH')
            
            ## Read the alignment file
            aln <- readGappedAlignments(BF,param=param)
            
            ## What type of RNA-Seq library are we dealing with?
            strand(aln) <- switch(lib.strand,
                                  sense = strand(aln),
                                  anti  = ifelse(strand(aln) == "+","-","+")
                                  )
            ## Removing aligned reads with more than one aligments
            aln <- aln[!values(aln)$NH > 1]
            
            gnModel <- gnModel[seqnames(gnModel) == seqname]
            anti.gnModel <- gnModel
            strand(anti.gnModel) <- ifelse(strand(gnModel)=="+","-","+")
            
            ## Only use aligments overlapping single genes
            aln <- aln[countOverlaps(aln, c(gnModel,anti.gnModel))==1 ]
            
            sense <- countOverlaps(gnModel, aln)
            anti <-  countOverlaps(anti.gnModel, aln)
            frac.sense <- sense/(sense+anti)
            
            frac.sense[(sense+anti) <= min.lim] <- NA
            names(frac.sense) <- names(gnModel)
            
            return(frac.sense)
        }
        
        chr <- unique(seqnames(genes))
        raw.res <- mcmapply(fracSense,
                            chr,
                            rep(BFL,each=length(chr)),
                            MoreArgs=list(genes),
                            mc.cores=CPU,
                            mc.preschedule=FALSE
                            )

        res <- lapply(split(raw.res,names(rep(BFL,each=length(chr)))),function(l) do.call(c,l))

        names(res) <- sub("\\.bam$","",basename(names(BFL)))
        return(res)
    }

getFeatCov <- function(feat,covs){
    chr.name <- as.character(unique(seqnames(feat)))
    strand <- as.character(strand(feat))
    covs <- covs[[strand]]
    
    if(start(feat) >= length(covs[[chr.name]])){
        return(Rle(NA,width(feat)))
    }
    start <- ifelse (start(feat)<1,1,start(feat))
    end <- ifelse(end(feat) > length(covs[[chr.name]])
                  ,length(covs[[chr.name]])
                  ,end(feat))
    
    cov <- seqselect(covs[[chr.name]],start=start,end=end)
    if(end(feat) > length(covs[[chr.name]])){
        pad <- end(feat) - sum(width(covs[[chr.name]]))
        cov <- c(cov,Rle(NA,pad))
    }
            if(start(feat)<1){
                pad <- 1-start(feat)
                cov <- c(Rle(NA,pad),cov)
            }
    ## If the transcript is on the other strand, reverse the order of coverage
    if(strand == '-'){cov <- rev(cov)}
            return(cov)
}


covAlongTx <-
    function(BFL,txdbfile,lib.strand=c("anti","sense"),min.lim=10,nCores=16,...){
        lib.strand <- match.arg(lib.strand)
        txdb<-loadDb(txdbfile)
        tx <- exonsBy(txdb,"tx")
        
        seqlevels(tx,force=TRUE) <- seqlevels(seqinfo(BFL[[1]]))
        seqlevels(tx,force=TRUE) <- seqlevels(tx)[!isCircular(tx)]
        
        txRelCov <-function(chr,BF,gnModel){
            gnModel <- gnModel[seqnames(unlist(range(gnModel))) == chr]
            
            s.i <- seqinfo(BF) 
            seq.length <- seqlengths(s.i)[chr]
            param <- ScanBamParam(which=GRanges(chr, IRanges(1, seq.length)),
                                  tag='NH')
            
            ## Read the alignment file
            aln <- readGappedAlignments(BF,param=param)
            ## fllip the aln strd if strandness is anti
            if (lib.strand == 'anti'){strand(aln) <- ifelse(strand(aln) == "+","-","+")}
            ## Removing aligned reads with more than one aligments
            aln <- aln[!values(aln)$NH > 1]
            
            ## Do not use aligments overlapping more than one genes
            aln <- aln[countOverlaps(aln, gnModel)==1 ]
            
            counts <- countOverlaps(gnModel, aln)

            gnModel <- gnModel[counts >= min.lim]
            exons <- unlist(gnModel)
            
            split.aln <- unlist(grglist(aln))
            split.aln <- split.aln[width(split.aln) > 3]
            ## Only keeps reads fully over
            exons.aln <- split.aln[unique(queryHits(findOverlaps(split.aln,exons,type='within')))]
            
            covs <- list('+'=coverage(exons.aln[strand(exons.aln) == '+']),
                         '-'=coverage(exons.aln[strand(exons.aln) == '-']))
            
            exon.cov <- lapply(seq_along(exons),function(i) getFeatCov(exons[i],covs))
            
            tx.cov <- sapply(split(exon.cov,rep(seq_along(gnModel),elementLengths(gnModel))),
                             function(cov) do.call(c,cov))

            tx.rel.cov <- mapply(function(x,y) x/y, tx.cov,counts[counts >= min.lim])
            names(tx.rel.cov) <- names(gnModel)
            tx.rel.cov
        }
        
        chr <- seqnames(seqinfo(tx))
        raw.res <- mcmapply(txRelCov,
                            chr,
                            rep(BFL,each=length(chr)),
                            MoreArgs=list(tx),
                            mc.cores=nCores,
                            mc.preschedule=FALSE,
                            SIMPLIFY=FALSE
                            )

        mat <- mclapply(do.call(c,raw.res),function(t) {
            start <- round(seq(1,length(t),length.out=101))[-101]
            end <- c(start[-1]-1,length(t))
            apply( data.frame(start,end),1,function(x) mean(t[x[1]:x[2]]))
        },mc.cores=nCores)
        

        skeleton <- sapply(split(raw.res,names(rep(BFL,each=length(chr)))),function(x) sum(sapply(x,length)))

        res <- lapply(split(mat,rep(names(skeleton),skeleton)),function(x) t(do.call(cbind,x)))

        names(res) <- sub("\\.bam$","",basename(names(BFL)))
        return(res)
    }

covsPerChr <-
    function(chr,bam.path,lib.strand=c("none","sense","anti"),...){
        
        lib.strand <- match.arg(lib.strand)
        
        s.i <- seqinfo(BamFile(bam.path)) 
        seq.length <- seqlengths(s.i)[chr]
        param <- ScanBamParam(which=GRanges(chr, IRanges(1, seq.length)),
                              what='mapq',
                              tag='NH')

        ## Read the alignment file
        aln <- unlist(grglist(readGappedAlignments(bam.path,param=param)))
        
        ## What type of RNA-Seq library are we dealing with?
        strand(aln) <- switch(lib.strand,
                              none = "*",
                              sense = strand(aln),
                              anti  = ifelse(strand(aln) == "+","-","+")
                              )

        ## Compute the coverage
        if (lib.strand == "none"){
            covs <- coverage(aln)[[chr]]
        } else {
            covs <- list('+'=coverage(aln[strand(aln)=='+'])[[chr]],
                         '-'=coverage(aln[strand(aln)=='-'])[[chr]]
                         )
        }
    }

bams2Covs <-
    function(BFL,lib.strand=c("none","sense","anti"),nCores=16,...){
        lib.strand <- match.arg(lib.strand)
        
        chrs <- sapply(BFL,function(BFL) seqnames(seqinfo(BFL)))
        bam.paths <- as.vector(sapply(BFL,function(BFL) rep(path(BFL),length(seqnames(seqinfo(BFL))))))
        
        covs <- mcmapply(covsPerChr,chrs,bam.paths,
                         MoreArgs=list(lib.strand),
                         SIMPLIFY=FALSE,
                         mc.cores=nCores,
                         mc.preschedule=FALSE
                         )
        
        covs <- lapply(split(covs,bam.paths),function(covs){
            if (lib.strand == 'none'){
                covs <- RleList(covs)
            } else {
                covs <- list('+' = RleList(lapply(covs,function(x) x$`+`)),
                             '-' = RleList(lapply(covs,function(x) x$`-`)))
                
            }
        })
        names(covs) <- gsub("\\.bam$","",basename(path(BFL)))
        return(covs)
    }

bams2bw <-
    function(BFL,destdir=c("bigwig"),lib.strand=c("none","sense","anti"),nCores=16,...){

        dir.create(destdir,FALSE,TRUE)
        
        covs <- bams2Covs(BFL,lib.strand,nCores)
        covs.GR <- mclapply(unlist(covs),as,'GRanges',
                            mc.cores=nCores,
                            mc.preschedule=FALSE)

        if(lib.strand == 'none'){
            bw <- paste(names(covs),"bw",sep=".")
        } else {
            bw <- paste(paste(rep(names(covs),sapply(covs,length)),
                              ifelse(as.vector(sapply(covs,names))=='+','p','m'),sep="_"),
                        "bw",sep=".")
        }

        bw.path <- file.path(destdir,bw,3)
        mclapply(seq_along(covs.GR),function(i) export(covs.GR[[i]],bw.path[[i]]),
                 mc.cores=nCores,
                 mc.preschedule=FALSE)
    }

                                               
