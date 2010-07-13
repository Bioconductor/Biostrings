### =========================================================================
### MultipleAlignment objects
### -------------------------------------------------------------------------
###

## setClass("MultipleAlignment",
##     representation(set="BStringSet",
##                    masks="NormalIRanges",
##                    rowMasks="integer"))

## virtual class for metadata
setClass("AlignmentMetadata",
    representation=representation(
        "VIRTUAL"))

## TODO: ... classes for each aligner


## virtual class for multiple alignments
setClass("MultipleAlignment",
         representation("VIRTUAL",
## TODO:                       metadata="AlignmentMetadata",
                        unmasked="XStringSet",
                        alignmentMask="NormalIRanges",
                        sampleMask="NormalIRanges"))

## concrete classes for multiple alignments
setClass("DNAMultipleAlignment",
    contains="MultipleAlignment",
    representation=representation(
        unmasked="DNAStringSet"))

## setClass("RNAMultipleAlignment",
##     contains="MultipleAlignment",
##     representation=representation(
##         unmasked="RNAStringSet"))

## setClass("AAMultipleAlignment",
##     contains="MultipleAlignment",
##     representation=representation(
##         unmasked="AAStringSet"))


## DNAMultipleAlignment constructor
DNAMultipleAlignment <- function(x=character(), start=1, end=nchar(x[[1]]),use.names=TRUE)
{
  ## Test that all lengths are == to each other
  if(length(x)>0){
   if(!all(unlist(lapply(x,nchar)) %in% nchar(x[[1]]))){
      stop("All the Strings in an MAS Set have to be the same length")}
  }
  new("DNAMultipleAlignment",
## TODO:     metadata=
      unmasked=DNAStringSet(x=x,start=start,end=end,width=NA,
        use.names=use.names),
      alignmentMask=as(IRanges(), "NormalIRanges"),
      sampleMask=as(IRanges(), "NormalIRanges"))
}


##Accessor methods:
setGeneric("alignmentMask", signature="x",
           function(x) standardGeneric("alignmentMask"))
setMethod("alignmentMask", "MultipleAlignment",
          function(x) x@alignmentMask)  
setGeneric("alignmentMask<-", function(x,value)
           standardGeneric("alignmentMask<-"))
setReplaceMethod("alignmentMask", signature(x="MultipleAlignment",
  value="NULL"), function(x, value) {x@alignmentMask<-as(IRanges(),
                                       "NormalIRanges"); x} )
setReplaceMethod("alignmentMask", signature(x="MultipleAlignment",
  value="NormalIRanges"), function(x, value) {x@alignmentMask<-value; x})


setGeneric("sampleMask", signature="x",
           function(x) standardGeneric("sampleMask"))
setMethod("sampleMask", "MultipleAlignment",
          function(x) x@sampleMask) 
setGeneric("sampleMask<-", function(x,value) standardGeneric("sampleMask<-"))
setReplaceMethod("sampleMask", signature(x="MultipleAlignment",
  value="NULL"), function(x, value) {x@sampleMask<-as(IRanges(),
                                       "NormalIRanges"); x} )
setReplaceMethod("sampleMask", signature(x="MultipleAlignment",
  value="NormalIRanges"), function(x, value) {x@sampleMask<-value; x})

setGeneric("unmask",signature="x", function(x) standardGeneric("unmask"))
setMethod("unmask", signature="MultipleAlignment",
          function(x) {alignmentMask(x) <- NULL;
                       sampleMask(x) <- NULL; x})



setMethod("as.character","MultipleAlignment",function(x){as.character(x@unmasked)})
setMethod("rownames","MultipleAlignment",function(x){names(x@unmasked)})
setReplaceMethod("rownames", "MultipleAlignment",
                 function(x,value){names(x@unmasked)<-value; x})

setMethod("ncol","MultipleAlignment",function(x){nchar(x@unmasked)})
setMethod("nrow","MultipleAlignment",function(x){length(x@unmasked)})




######################################################################
## Below here we will implement versions of methods that are "mask aware"
######################################################################


## I definitely want to OL these to consider the mask slot
setMethod("consensusMatrix","MultipleAlignment",
          function(x){consensusMatrix(x@unmasked)})
setMethod("consensusString","MultipleAlignment",
          function(x){consensusString(x@unmasked)})
setMethod("alphabetFrequency","MultipleAlignment",
          function(x){alphabetFrequency(x@unmasked)})






########################################################################## 
## stuff for the show method

## TODO: For the show method, We need to make sure that we display a special
## message when the entire thing is masked...

##msk is a normalIRanges object which *IS* a single mask from a maskCollection.
.charStrFormat = function(str, msk){
    msk <- Mask(mask.width=nchar(str),
                start=start(msk),end=end(msk),width=width(msk))
    masks(str) <- msk
    ans <- as.character(unmasked(str))
    nir0 <- as(str, "NormalIRanges")
    for (i in seq_len(length(nir0))) {
      strip <- paste(rep.int("#", width(nir0)[i]), collapse="")
      substr(ans,  start(nir0)[i], end(nir0)[i]) <- strip
    }
    ans = gsub("#","",ans)
    ##ans = gsub("-","*",ans)##Temp test
    ans
}


SeqSnippetSlice <- function(x, width, index, mask)
{
    if(index<=1){start=1}else{start=width*(index-1)}
    if (width < 7L)
        width <- 7L
    seqlen <- length(x)
    truncEnd <- start+width
    if ((truncEnd) > seqlen && seqlen > start) {
      paste(subseq(.charStrFormat(x, mask), start=start, end=seqlen),
            sep="")      
    } else if ((truncEnd) <= seqlen) {
      paste(subseq(.charStrFormat(x, mask), start=start, width=width),
            sep="")            
    }else{}## print nothing
}


.namesMASW <- 20
.XStringMAS.show_frame_header <- function(iW, widthW, with.names)
{
    cat(format("", width=iW+1),
        ## format("width", width=widthW, justify="right"),
        sep="")
    if (with.names) {
        cat(format(" seq", width=getOption("width")-iW-widthW-.namesMASW-1),
            format("names", width=.namesMASW, justify="left"),
            sep="")
    } else {
        cat(" seq")
    }
    cat("\n")
}

.XStringMAS.show_frame_line <- function(x, i, iW, widthW, index, hiddenRows)
{
   if(!i %in% hiddenRows){## only draw rows that are not hidden
    width <- ncol(x)[i]
    snippetWidth <- getOption("width") - 2 - iW - widthW
    if (!is.null(rownames(x)))
        snippetWidth <- snippetWidth - .namesMASW - 1
    seq_snippet <- SeqSnippetSlice(x@unmasked[[i]], snippetWidth, index, mask=alignmentMask(x))
    
    if (!is.null(rownames(x)) && !is.null(seq_snippet))
      seq_snippet <- format(seq_snippet, width=snippetWidth)
    cat(format(paste("[", i,"]", sep=""), width=iW, justify="right"), " ",
        " ",
        seq_snippet,
        sep="")
    if (!is.null(rownames(x)) && !is.null(seq_snippet)) {
      snippet_name <- rownames(x)[i]
      if (is.na(snippet_name))
        snippet_name <- "<NA>"
      else if (nchar(snippet_name) > .namesMASW)
        snippet_name <- paste(substr(snippet_name, 1, .namesMASW-3),
                              "...", sep="")
      cat(" ", snippet_name, sep="")
    }
    cat("\n")
  }
}


.XStringMAS.show_frame <- function(x, half_nrow=9L)
{
  ##hard codes the number of sets of rows we plan to display
  .frameLines <- 3
  ##gather hidden rows
  .hiddenRows <- as.integer(sampleMask(x))

    lx <- nrow(x)
    iW <- nchar(as.character(lx)) + 2 # 2 for the brackets
    ncharMax <- max(ncol(x))
    widthW <- max(nchar(ncharMax), nchar("width"))
    .XStringMAS.show_frame_header(iW, widthW, !is.null(rownames(x)))

  for(index in seq_len(.frameLines)){   
    if (lx <= 2*half_nrow+1) {
      for (i in seq_len(lx))
        .XStringMAS.show_frame_line(x, i, iW, widthW, index, .hiddenRows)
        cat("\n")
    } else {
      for (i in 1:half_nrow)
        .XStringMAS.show_frame_line(x, i, iW, widthW, index, .hiddenRows)
      cat(format("...", width=iW, justify="right"),
          format("...", width=widthW, justify="right"),
          "...\n")
      for (i in (lx-half_nrow+1L):lx)
        .XStringMAS.show_frame_line(x, i, iW, widthW, index, .hiddenRows)
        cat("\n")
    }
  }
}


setMethod("show", "MultipleAlignment",
    function(object)
    {
        cat("  A ", class(object), " instance with ", nrow(object), " sequences represented.","\n", sep="")
        if (nrow(object) != 0)
            .XStringMAS.show_frame(object)
    }
)




##########################################################################
## more specialized functions

## function to generate this from a clustal file.

clustReader = function(file){
  con = file(file)
  data = readLines(con)
  close(con)  
  ## drop the header and empty lines
  data = data[grep("^CLUSTAL", data, invert=TRUE)]
  data = data[data!=""]
  ## get the index of the 1st line to be a complete blank line
  count = grep("^(\\s|\\*)+$", data)[1] - 1
  data = data[grep("^(\\s|\\*)+$",data, invert=TRUE)]
  ## Therefore we shall gather and then drop the IDs
  ids =  gsub("(^gi\\S+)\\s+.+","\\1",data)[1:count]
  data = gsub("^gi\\S+\\s+","",data)
  ## And also the positions from the end
  data = gsub("\\s+\\d+$","",data)
  ## And our count tells us how to split these up
  data = split(data, factor(rep(1:count, length(data)/count)))
  data = unlist(lapply(data, paste, collapse=""))
  names(data) = ids
  ##TODO: there needs to be a more elegant way than hard coding this choice!
  ##So plan to add a parameter or else discuss a more elegant way for the
  ##constructor to "figure it out"
  DNAMultipleAlignment(data) 
}



