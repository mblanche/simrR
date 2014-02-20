chrToUCSC <- function(object){
    renameSeqlevels(object,structure(paste0("chr",seqlevels(object)),names=seqlevels(object)))
}

chrFromUCSC <- function(object){
    renameSeqlevels(object,structure(sub("^chr","",seqlevels(object)),names=seqlevels(object)))
}

counterPerChr <-
    function(i,chrs,BFLs,gnModel,mapq=10,allow.multi=FALSE,lib.strand=c("none","sense","anti"),...){
        lib.strand <- match.arg(lib.strand)

        BF <- BFLs[[i]]
        chr <- chrs[[i]]
                
        seq.length <- seqlengths(BF)[chr]
        param <- ScanBamParam(which=GRanges(chr, IRanges(1, seq.length)),
                              what='mapq',
                              tag='NH')

        ## Read the alignment file
        aln <- readGAlignments(BF,param=param)
        
        ## What type of RNA-Seq library are we dealing with?
        strand(aln) <- switch(lib.strand,
                              none = "*",
                              sense = strand(aln),
                              anti  = ifelse(strand(aln) == "+","-","+")
                              )
        ## Removing aligned reads with more than one aligments
        if(!allow.multi) aln <- aln[!values(aln)$NH > 1]
        ## Removing the low map quality alignmentsw
        if(!is.null(mapq)) aln <- aln[!values(aln)$mapq <= mapq]
        ## Subseting our gnModel to only models on the curent chromsome
        if (class(gnModel) == "GRangesList"){
            gnModel <- gnModel[seqnames(unlist(range(gnModel))) == chr]
        } else if (class(gnModel) == "GRanges"){
            gnModel <- gnModel[seqnames(gnModel) == chr]
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
                          mapq=10,
                          allow.multi=FALSE,
                          nCores=16,
                          as.list=FALSE,
                          ...){
    lib.strand <- match.arg(lib.strand)
    ## Making sure there is an index
    parIndexBam(BFL,nCores)
    ## Prepare the jobs to send to the multiCores
    chrs <- as.vector(sapply(BFL,seqlevels))
    BFL.chrs <- rep(BFL,elementLengths(lapply(BFL,seqlevels)))
    ## Reading the bam files, in parallel, one chromosome at a time
    counts.raw <- mclapply(seq_along(chrs)
                           ,counterPerChr
                           ,chrs
                           ,BFL.chrs
                           ,gnModel=gnModel
                           ,lib.strand=lib.strand
                           ,mc.cores=nCores
                           ,mc.preschedule=FALSE
                           ,mapq=mapq
                           ,allow.multi=allow.multi
                           ,...
                           )
    ## Reorganizing the lists of counts per chromosome, returning a SummearizedExperiment as default or a list
    ## As before
    counts <- mapply(function(cnts){to.ret <- do.call(c,cnts)
                                    to.add <- rep(0,sum(!names(gnModel) %in% names(to.ret)))
                                    names(to.add) <- names(gnModel)[!names(gnModel) %in% names(to.ret)]
                                    c(to.ret,to.add)
                                }
                     ,split(counts.raw,path(BFL.chrs))
                     ,SIMPLIFY=!as.list)

    ## package the object (or return a list)
    if (as.list){
        names(counts) <- sub("\\.bam$","",basename(names(BFL)))
    } else {
        colData <- DataFrame(files=file.path(normalizePath(dirname(names(BFL))),basename(names(BFL))),
                             row.names=sub("\\.bam$","",basename(names(BFL))))
        counts <- SummarizedExperiment(assays=SimpleList(counts=counts)
                                       ,colData=colData)
        if(!is.null(rownames(counts))) rowData(counts) <- gnModel[rownames(counts)]
    }
    return(counts)
}

uxonify <- function(txdb){ reduce(exonsBy(txdb,'gene')) }

ixonify <- function(txdb){
  uxons <- uxonfiy(txdb)
  genes <- unlist(range(uxons))
  psetdiff(genes,uxons)
}

parIndexBam <- function(BFL,nCores=16,...){
    BFL <- BFL[!file.exists(paste(path(BFL),".bai",sep=''))]
    if(length(BFL) == 0) return(NULL)
    ## Make sure the bam files are indexed
    mclapply(BFL,indexBam,mc.cores=nCores,mc.preschedule=FALSE)
    return(NULL)
}

topHatReport <- function(topHatDir,nCores=16){

  topHat.mapped <- BamFileList(list.files(path=topHatDir,pattern="^accepted_hits.bam$",recu=TRUE,full.names=TRUE))
  topHat.unmapped <- BamFileList(file.path(dirname(path(topHat.mapped)),'unmapped.bam'))

  stopifnot(length(topHat.mapped) == length(topHat.unmapped))
            
  parIndexBam(c(topHat.mapped,topHat.unmapped),nCores)
  
  toProcess <-
    do.call(rbind,lapply(topHat.mapped,function(file){data.frame(chr=seqnames(seqinfo(file)),
                                                                 file=path(file),
                                                                 stringsAsFactors=FALSE)
                                                    }))
  counter <-function(d){
      fl <- d$file
      seqname <- d$chr
      
      s.i <- seqinfo(BamFile(fl)) 
      seq.length <- seqlengths(s.i)[seqname]
      param <- ScanBamParam(which=GRanges(seqname, IRanges(1, seq.length)),
                            flag = scanBamFlag(isUnmappedQuery=NA)
                            ,tag='NH')
      
      ## Read the alignment file
      aln <- readGAlignments(fl,param=param)
      
      tab.cnt <- table(values(aln)$NH)
      
      ## normalize the counts per their number of hits
    tab.cnt <- tab.cnt/as.numeric(names(tab.cnt))
      
    counts <- c(unique=tab.cnt[as.numeric(names(tab.cnt))==1]/1e6,
                multi.match=sum(tab.cnt[as.numeric(names(tab.cnt))>1])/1e6)
  }
  ## map the counts per chromosome per file
  mapped.data <- mclapply(split(toProcess,1:nrow(toProcess)),
                          counter,
                          mc.cores=nCores,
                          mc.preschedule=FALSE
                          )
  
  ## Reduce to 1) chromosomes 2) files
  reduced.data <- lapply(split(mapped.data,toProcess$file),function(perFile){
      d <- do.call(rbind,perFile)
      colSums(d)
  })
  
  unmapped.reads <- mclapply(topHat.unmapped,function(BF) countBam(BF)$records/1e6,
                             mc.cores=nCores,
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
    function(BFL,txdb,min.reads=50,lib.strand=c("anti","sense","none"),nCores=16,...){
        lib.strand <- match.arg(lib.strand)
        genes <- unlist(range(exonsBy(txdb,"gene")))
        ##genes <- genes[seqnames(genes) %in% seqnames(seqinfo(BFL[[1]]))]
        
        fracSense <-function(id,chrs,BFL,gnModel){
            chr <- chrs[id]
            BF <- BFL[[id]]
            ## If reading a chr without genes, return imediately
            if(!chr %in% seqlevels(gnModel)) return()

            ## Otherwise, read the bam on that chr
            param <- ScanBamParam(which=GRanges(chr, IRanges(1, seqlengths(seqinfo(BF))[chr])),
                                  tag='NH')
            aln <- readGAlignments(BF,param=param)
            
            ## What type of RNA-Seq library are we dealing with?
            strand(aln) <- switch(lib.strand,
                                  none = strand(aln),
                                  sense = strand(aln),
                                  anti  = ifelse(strand(aln) == "+","-","+")
                                  )
            ## Removing aligned reads with more than one aligments
            aln <- aln[is.na(values(aln)$NH) | !values(aln)$NH > 1]
            
            gnModel <- gnModel[seqnames(gnModel) == chr]
            anti.gnModel <- gnModel
            strand(anti.gnModel) <- ifelse(strand(gnModel)=="+","-","+")
            
            ## Only use aligments overlapping single genes
            aln <- aln[countOverlaps(aln, c(gnModel,anti.gnModel))==1 ]
            
            sense <- countOverlaps(gnModel, aln)
            anti <-  countOverlaps(anti.gnModel, aln)
            frac.sense <- sense/(sense+anti)
            ## Only report for gene with at least min.lim reads overlapping
            frac.sense[(sense+anti) <= min.reads] <- NA
            names(frac.sense) <- names(gnModel)
            
            return(frac.sense)
        }
        
        chrs <- as.vector(sapply(BFL,seqlevels))
        BFL2chrs <- rep(BFL,sapply(BFL,function(BF) length(seqlevels(BF))))
        
        raw.res <- mclapply(seq_along(chrs),
                            fracSense,
                            chrs,
                            BFL2chrs,
                            genes,
                            mc.cores=nCores,
                            mc.preschedule=FALSE
                            )
        
        res <- lapply(split(raw.res,path(BFL2chrs)),function(l) do.call(c,l))
        
        names(res) <- sub("\\.bam$","",basename(path(BFL)))
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

getTxRelCov<- 
    function(BFL,tx,lib.strand=c("anti","sense","none"),min.lim=50,nCores=16,n=100,...){
        meta.covs <- featCovViews(BFL,tx,lib.strand,min.lim,nCores,...)
        ##Remove chromsomes with no elements
        meta.covs <- lapply(meta.covs,function(view){view[sapply(view,length) > 0]})
        ##################################################
        ## Helper functions
        ##################################################
        ## This breaks a region define by start and width in n fragmetns of almost
        ## equal size. Returns the broken parts as a IRanges object to be reatach to teh Views objet
        getBlocks <- Vectorize(function(start,width,n){
            IR.width <- rep(width %/% n,n) + c(rep(1,width %% n),rep(0,n-(width %% n)))
            starts <- cumsum(c(start,IR.width))
            IRanges(start = starts[-length(starts)],
                    width=IR.width)
        },c('start','width'))
        ##################################################
        ## Break the views in coverages of n parts
        ##################################################
        fragFeats <- function(view,n=100,nCores=16){
            ## First, fragment each IRanges in equal parts of almost equal size
            IRL <- IRangesList(mclapply(view,function(v){
                do.call(c,getBlocks(start(v),width(v),n))
            }
                                        ,mc.cores=nCores
                                        ,mc.preschedule=FALSE)
                               )
            Views(subject(view),IRL)
        }
        ##################################################
        ## Return the sum of the framgent coverage
        ##################################################
        getFragSum <- function(RVL,nCores=16){
            RleList(mclapply(RVL,function(view) Rle(viewSums(view)),mc.cores=nCores,mc.preschedule=FALSE))
        }
        ##################################################
        ## Compute the relative coveage
        ##################################################
        relCovs <- function(RLL.FRAG,RLL.FEAT){
            RLL.REL.COV <- RLL.FRAG/RleList(lapply(RLL.FEAT,rep,each=n))
            res <- lapply(RLL.REL.COV,function(RL){if(length(RL) == 0){return()}
                                                   matrix(as.vector(RL),ncol=n,byrow=TRUE)})
            res <- do.call(rbind,res[sapply(res,function(x) !is.null(dim(x)))])
        }
        ## This is where the meat is!
        frag.covs <- lapply(meta.covs,fragFeats)
        frag.sums <- lapply(frag.covs,getFragSum)
        feat.sums <- lapply(meta.covs,getFragSum)
        
        res <- mapply(relCovs,frag.sums,feat.sums,SIMPLIFY=FALSE)
    }

featCovViews2 <-
    function(BFL,gnModel,coverage.loading=c('bySeqname','byFile'),lib.strand=c("anti","sense","none"),nCores=16,...){
        lib.strand <- match.arg(lib.strand)
        coverage.loading <- match.arg(coverage.loading)
        
        ## Map the data for the coverages
        if(coverage.loading == 'byFile'){
            covs <- covByFile(BFL,lib.strand,nCores)
        }else{
            covs <- covByChr(BFL,lib.strand,nCores)
        }
        ## These two operations can be time consuming, let's compute them outside the inner loop
        gene.strand <- strand(unlist(range(gnModel)))
        gene.seqnames <- seqnames(unlist(range(gnModel)))

        ## Compute coverage views for spliced transcript in the gnModel
        views.raw <- mclapply(covs,function(cov){
            ## No need to go further if the visited seqnames is empty
            ## Might want to prefilter before... might be faster
            if(sum(cov$coverage) == 0) return(Views(Rle(),IRanges()))
            ## Extract the only the visited features
            sub.gn <- gnModel[gene.strand==cov$strand & gene.seqnames==cov$seqname]
            ## Extract the exons as GRanges
            exons <- unlist(sub.gn)
            ## Compute the Views for every exons
            exon.views <- Views(cov$coverage,ranges(exons))
            ## Keep only the transcripts with at least min.lim coverage
                                        #exon.views <- exon.views[sum(exon.views) >= min.lim]
            names(exon.views) <- NULL
            ## Return an RLE of the splice coverage of all each exons
            ## Flip the coverages if we are dealing with feature on the negative strand
            exon.covs <- do.call(c,viewApply(exon.views,function(x){if(cov$strand=='-'){rev(x)}else{x}}))
            ## compute the location of each transcripts (spliced exons) on the new spliced coverage
            tx.width <- sapply(split(width(exons),rep(seq_along(sub.gn),elementLengths(sub.gn))),sum)
            ## Finish by building a new view of the spliced coverage
            if (length(exon.covs)==0){
                ## If nothing to return, return an empty RleViews
                return(Views(Rle(),IRanges()))
            } else {
                ## Otherwise, return a views of the concat tx coverage
                tx.views <- Views(exon.covs,
                                  IRanges(start=c(1,cumsum(tx.width[-length(tx.width)])+1),width=tx.width))
                names(tx.views) <- names(sub.gn)
                return(tx.views)
            }
        },mc.cores=nCores,mc.preschedule=ifelse(length(covs)/nCores > 100,TRUE,FALSE))
        
        ## Reduce the data to the original BFL
        res <- mclapply(split(views.raw,sapply(covs,function(x) path(x$BF))),function(views){
            ## recover all the coverages into a single coverage
            names(views) <- NULL
            subject <- do.call(c,lapply(views,subject))
            ## Recover all the different ranges
            ranges <- lapply(views,ranges)
            ## recover the transcript names
            tx.names <- do.call(c,lapply(views,names))
            ## Compute the amount of shifts the second to the last ranges need to be incremented
            shifts <- cumsum(lapply(views,function(x) length(subject(x)))[-length(views)])
            ## Shifts the ranges and recover a single final range
            shifted <- do.call(c,c(list(ranges(views[[1]])),mapply(shift,ranges[-1],shifts,SIMPLIFY=FALSE))) 
            ## Assemble a single Views object
            res <- Views(subject,shifted)
            names(res) <- tx.names
            return(res)
        },mc.preschedule=FALSE,mc.cores=nCores)
        
    ## rename the final list of views
    names(res) <- sub("\\.bam","",basename(path(BFL)))
    ## return as a fancy Bioc object!
    return(RleViewsList(res))
}

covByChr <- function(BFL,lib.strand,nCores=16){
    if(lib.strand=='none'){
        isMinus <- NA
    }else{
        isMinus <- c(FALSE,TRUE)
    }
    chrs <- rep(as.vector(sapply(BFL,seqlevels)),each=length(isMinus))
    BFLs <- rep(BFL,sapply(BFL,function(BF) length(isMinus)*length(seqlevels(BF))))
    strands <- rep(isMinus,sum(sapply(BFL,function(BF) length(seqlevels(BF)))))
    covs.raw <- mclapply(seq_along(BFLs),function(i){
        BF <- BFLs[[i]]
        chr <- chrs[[i]]
            ## Setting up the chromosome region to read
        seq.length <- seqlengths(BF)[chr]
            param <- ScanBamParam(flag  = scanBamFlag(isMinusStrand=strands[[i]]),
                                  which = GRanges(chr, IRanges(1, seq.length)))
        ## Read the alignment file
        aln <- unlist(grglist(readGAlignmentsFromBam(BF,param=param)))
        ## Computing the coverages on the seqlevels
        cov <- coverage(aln)[[chr]]
        ## Return the results as a named list... (Could use S4 obj)
        list(coverage=cov,
             seqname=chr,
             strand =switch(lib.strand,
                 none = '*',
                 sense = ifelse(strands[[i]],'-','+'),
                 anti  = ifelse(strands[[i]],'+','-')),
             BF=BF)
    },mc.cores=nCores,mc.preschedule=FALSE)
}


covByFile <- function(BFL,lib.strand,nCores=16){
    ## Define the strand that need to be read
    if(lib.strand=='none'){
        isMinus <- rep(NA,length(BFL))
    }else{
        isMinus <- rep(c(FALSE,TRUE),length(BFL))
    }
    ## Assemble a list of bam file to visit, augmented with stranded
    BFLs <- rep(BFL,each=length(unique(isMinus)))
    ## Read the bam, compute the coverage
    cov.raw <- mclapply(seq_along(isMinus),function(i){
        param <- ScanBamParam(scanBamFlag(isMinusStrand=isMinus[[i]]))
        ## Read the alignment file
        cov <- coverage(unlist(grglist(readGAlignmentsFromBam(BFLs[[i]],param=param))))
        ## Computing the coverages on the seqlevels
    },mc.cores=nCores,mc.preschedule=FALSE)
    covs <- do.call(c,cov.raw)
    cov.strands <- switch(lib.strand,
                          none = rep('*',sum(sapply(cov.raw,length))),
                          sense = ifelse(rep(isMinus,sapply(cov.raw,length)),'-','+'),
                          anti  = ifelse(rep(isMinus,sapply(cov.raw,length)),'+','-'),
                          )
    cov.chrs <- unlist(lapply(cov.raw,names))
    cov.BFL <- rep(BFLs,sapply(cov.raw,length))
    ## Return a list of named lists
    mapply(list,coverage=covs,seqname=cov.chrs,strand=cov.strands,BF=cov.BFL,SIMPLIFY=FALSE)
}


getCovsMatrix <- function(cov.view,cut.number=100,type=c('relative','sum.coverage'),nCores=15,min.coverage=(50*50)){
    ## Filter tx with less than min.coverage coverage
    cov.view <- cov.view[sum(cov.view)>=min.coverage]
    ##Remove the transcripts shorter than the number of bins... Not sure what to do with them
    cov.view <- cov.view[width(cov.view)>=cut.number]
    ## Breaks every transcripts in 100 bins, find the cuts
    cuts <- viewApply(cov.view,function(x) cut(seq_along(x),cut.number,labels=FALSE),simplify=FALSE)
    ## Retrieve the coverages from teh Views as vectors
    cov.vecs <- viewApply(cov.view,as.vector,simplify=FALSE)
    ## Compute the coverage sums for each cuts
    res <- do.call(rbind,mclapply(seq_along(cov.vecs), function(i) sapply(split(cov.vecs[[i]],cuts[[i]]),sum),mc.cores=nCores))
    ## get teh gene names for each row
    row.names(res) <- names(res)
    ## return in relative space if asked
    if(type=='relative') res <- res/rowSums(res)
    ## Return
    return(res)
}

getCovsMatrix2 <- function(cov.views,cut.number=100,type=c('relative','sum.coverage'),nCores=15,min.coverage=(50*50)){
    type <- match.arg(type)
    ## Filter tx with less than min.coverage coverage
    cov.views <- mclapply(cov.views,function(cov.view) cov.view[sum(cov.view)>=min.coverage],mc.cores=nCores,mc.preschedule=FALSE)
    ##Remove the transcripts shorter than the number of bins... Not sure what to do with them
    cov.views <- mclapply(cov.views,function(cov.view) cov.view[width(cov.view)>=cut.number],mc.cores=nCores,mc.preschedule=FALSE)
    ## if coverage is very sparse, some Views my by empty after the filtering. Keep track of them
    is.empty <- sapply(cov.views,length)==0
    ## Breaks every transcripts in 100 bins, find the cuts
    cuts <- mclapply(cov.views[!is.empty],function(cov.view){
        viewApply(cov.view,function(x) cut(seq_along(x),cut.number,labels=FALSE),simplify=FALSE)
    },mc.cores=nCores,mc.preschedule=FALSE)
    names(cuts) <- NULL
    cuts <- do.call(c,cuts)
    ## Retrieve the coverages from teh Views as vectors
    cov.vecs <- mclapply(cov.views[!is.empty],function(cov.view){
        viewApply(cov.view,as.vector,simplify=FALSE)
    },mc.cores=nCores,mc.preschedule=FALSE)
    names(cov.vecs) <- NULL
    cov.vecs <- do.call(c,cov.vecs)
    ## Compute the coverage sums for each cuts
    res.raw <- mclapply(seq_along(cov.vecs), function(i) sapply(split(cov.vecs[[i]],cuts[[i]]),sum),mc.cores=nCores)
    ## Need to reduce to the original cov.views
    res.t <- lapply(split(res.raw,rep(seq_along(cov.views[!is.empty]),sapply(cov.views[!is.empty],length))),do.call,what=rbind)
    ## get the gene names for each row
    res.t <- mapply(function(r,c.v){
        row.names(r) <- names(c.v)
        return(r)
    },res.t,cov.views[!is.empty])
    ## return in relative space if asked
    if(type=='relative') res.t <- lapply(res.t,function(res) res/rowSums(res))
    ## Create the list of resutls
    res <- list()
    res[is.empty] <- NA
    res[!is.empty] <- res.t
    ## Name the list with the original cov.views names
    names(res) <- names(cov.views)
    ## Return
    return(res)
}

featCovViews <-
    function(BFL,features,lib.strand=c("anti","sense","none"),min.lim=50,nCores=16,...){
        if(!is(BFL,'BamFileList')) stop("A BamFileList is required")
        if(!length(BFL)>0) stop("Need to pass at least one BamFile",call.=FALSE)
        lib.strand <- match.arg(lib.strand)

        BFL.seqlevs <- unique(as.vector(sapply(BFL,seqlevels)))
        
        ## is there at least one chr name in BFL found in features? if not, bail out...
        if (sum(BFL.seqlevs %in% seqlevels(features)) == 0){
            stop("No compatible seqnames between the bam files and features were found",call.=FALSE)
        }
        ## Make sure I will not stumble of problematic tx
        if (!all(seqlevels(features) %in% BFL.seqlevs)){
            seqlevs <- BFL.seqlevs[BFL.seqlevs %in% seqlevels(features)]
            seqlevels(features,force=TRUE) <- seqlevs
        }
        if(any(isCircular(features)[!is.na(isCircular(features))])) {
            seqlevels(features,force=TRUE) <- seqlevels(features)[!isCircular(features)]
        }
        
        computeCovs <-function(id,chrs,BFL.chrs,gnModel){
            chr <- chrs[id]
            BF <- BFL.chrs[[id]]
            gnModel <- gnModel[seqnames(unlist(range(gnModel))) == chr]
            s.i <- seqinfo(BF) 
            seq.length <- seqlengths(s.i)[chr]
            param <- ScanBamParam(which=GRanges(chr, IRanges(1, seq.length)),
                                  tag='NH')
            ## Read the alignment file
            aln <- readGAlignments(BF,param=param)
            ## fllip the aln strd if strandness is anti
            ## What type of RNA-Seq library are we dealing with?
            strand(aln) <- switch(lib.strand,
                                  none = "*",
                                  sense = strand(aln),
                                  anti  = ifelse(strand(aln) == "+","-","+")
                                  )
            ## Removing aligned reads with more than one aligments
            aln <- aln[!values(aln)$NH > 1]
            ## Split gapped alignmets
            split.aln <- unlist(grglist(aln))
            
            ##gnModel <- gnModel[countOverlaps(unlist(range(gnModel)),aln) > min.lim]
            gnModel <- gnModel[countOverlaps(gnModel,aln) > min.lim]
            gn.strand <- as.vector(strand(unlist(range(gnModel))))
            
            ## There is a ~3X gain in speed in calculting the coverage if only
            ## Subseting on reads overlaping the features (including the findOverlaps overhead)
            comps.aln <- split.aln[unique(queryHits(findOverlaps(split.aln,unlist(gnModel))))]
            if (lib.strand == 'none'){
                covs <- list('+'=(all.covs <- coverage(comps.aln)),
                             '-'=all.covs)
            }else {
                covs <- lapply(c('+','-'),function(s) coverage(comps.aln[strand(comps.aln) == s]))
            }
            names(covs) <- c('+','-')
            
            cov.views <- sapply(c('+','-'),function(s){
                cov <- covs[[s]][[as.character(chr)]]
                comps <- unlist(gnModel[gn.strand==s])
                comps.view <- Views(cov,ranges(comps))
                names(comps.view) <- NULL
                ## Return an RLE of the coverage for every feature parts in gnModel
                ## Flip the coverages if we are dealing with feature on the negative strand
                comps.cov <- do.call(c,viewApply(comps.view,function(cov){if(s=='-'){rev(cov)}else{cov}}))
                feat.width <- sapply(split(width(comps),rep(seq_along(gnModel[gn.strand==s]),
                                                            elementLengths(gnModel[gn.strand==s]))),sum)
                if (length(comps.cov)==0){
                    ## If nothing to return, return an empty RleViews
                    return(Views(Rle(),IRanges()))
                } else {
                    ## Otherwise, return a views of the concat tx coverage
                    feat.views <- Views(comps.cov,
                                        IRanges(start=c(1,cumsum(feat.width[-length(feat.width)])+1),width=feat.width))
                    names(feat.views) <- names(gnModel[gn.strand==s])
                    return(feat.views)
                }
            })
            Views(c(subject(cov.views$'+'),subject(cov.views$'-')),
                  c(ranges(cov.views$'+'),shift(ranges(cov.views$'-'),length(subject(cov.views$'+')))))
        }
        
        ## Simillarly, no need to visit seqnames in BFL not found in features
        ## No need to be wastefull, runs only on seqnames that are in the features DB
        chrs <- as.vector(sapply(BFL,function(BF) seqlevels(BF)[seqlevels(BF) %in% seqlevels(features)]))
        BFL.chrs <- rep(BFL,sapply(BFL, function(BF) length(seqlevels(BF)[seqlevels(BF) %in% seqlevels(features)])))
        
        raw.covs <- mclapply(seq_along(chrs),
                             computeCovs,
                             chrs,
                             BFL.chrs,
                             features,
                             mc.cores=nCores,
                             mc.preschedule=FALSE
                             )
        
        res <- mapply(function(v,names){ r <- RleViewsList(v)
                                         names(r) <- names
                                         return(r)
                                     }
                      ,split(raw.covs,path(BFL.chrs))
                      ,split(chrs,path(BFL.chrs))
                      )

        names(res) <- sub("\\.bam$","",basename(names(BFL)))
        return(res)
    }

covsPerChr <-
    function(id,chrs,BFL,lib.strand=c("none","sense","anti")){
        chr <- chrs[id]
        BF <- BFL[[id]]
        
        lib.strand <- match.arg(lib.strand)
        
        s.i <- seqinfo(BF)
        seq.length <- seqlengths(s.i)[chr]
        param <- ScanBamParam(which=GRanges(chr, IRanges(1, seq.length)),
                              what='mapq',
                              tag='NH')

        ## Read the alignment file
        aln <- unlist(grglist(readGAlignments(BF,param=param)))
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
    function(BFL,lib.strand=c("none","sense","anti"),nCores=16){
        if(!is(BFL,'BamFileList')) stop("A BamFileList is required")
        if(!length(BFL)>0) stop("Need to pass at least one BamFile")
        lib.strand <- match.arg(lib.strand)
        
        chrs <- as.vector(sapply(BFL,seqlevels))
        BFL2chrs <- rep(BFL,sapply(BFL,function(BF) length(seqlevels(BF))))
                
        covs <- mclapply(seq_along(chrs),
                         covsPerChr,
                         chrs,
                         BFL2chrs,
                         lib.strand,
                         mc.cores=nCores,
                         mc.preschedule=FALSE
                         )
        names(covs) <- chrs
        
        covs <- lapply(split(covs,path(BFL2chrs)),function(covs){
            if (lib.strand == 'none'){
                covs <- RleList(covs)
            } else {
                covs <- list('+' = RleList(lapply(covs,function(x) x$`+`),compress=FALSE),
                             '-' = RleList(lapply(covs,function(x) x$`-`),compress=FALSE))
            }
        })
        names(covs) <- gsub("\\.bam$","",basename(path(BFL)))
        return(covs)
    }

bams2bw <-
    function(BFL,destdir=c("bigwig"),lib.strand=c("none","sense","anti"),lib.norm=TRUE,nCores=16){
        if(!is(BFL,'BamFileList')) stop("A BamFileList is required")
        if(!length(BFL)>0) stop("Need to pass at least one BamFile")
        lib.strand <- match.arg(lib.strand)
                

        if(!is(BFL,'BamFileList')) stop("A BamFileList is required")
        
        if(!length(BFL)>0) stop("Need to pass at least one BamFile")

        ## Create bam index if they don't exits
        parIndexBam(BFL,nCores)
        
        dir.create(destdir,FALSE,TRUE)
        
        covs <- bams2Covs(BFL,lib.strand,nCores)
        
        covs.GR <- mclapply(unlist(covs),as,'GRanges',
                            mc.cores=nCores,
                            mc.preschedule=FALSE)
        ## Ok, every even slot in the covs.GR list is of the negative strand (ie 1:lenght(list)%%2)
        ## Flip the sign on the coverage value if the library is stranded
        if(lib.strand!='none'){
            anti.cov <- covs.GR[seq_along(covs.GR)%%2 == 0]
            anti.cov <- lapply(anti.cov,function(anti.cov){
                values(anti.cov)$score <- -values(anti.cov)$score
                anti.cov
            })
            covs.GR[seq_along(covs.GR)%%2 == 0] <- anti.cov
        }
        ## Normalize the coverage over the library size
        if(lib.norm){
            covs.GR <- lapply(covs.GR,function(c){
                values(c)$score <- values(c)$score/(abs(sum(as.numeric(values(c)$score)))/1e8)
                return(c)
            })
        }
        ## Creating file names and file paths        
        if(lib.strand == 'none'){
            bw <- paste(names(covs),"bw",sep=".")
        } else {
            suffix <- paste0("_",c('p','m'),'.bw')
            bw <- as.vector(sapply(sub("\\.bam$","",basename(names(BFL))),paste0,suffix))
        }
        bw.path <- file.path(destdir,bw)
        ## Exporting as bigwigs
        mclapply(seq_along(covs.GR),function(i) export(covs.GR[[i]],bw.path[[i]]),
                 mc.cores=nCores,
                 mc.preschedule=FALSE)
    }
