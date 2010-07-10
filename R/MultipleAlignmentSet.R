### =========================================================================
### MultipleAlignmentSet objects
### -------------------------------------------------------------------------
###

setClass("MultipleAlignmentSet",
    representation(set="BStringSet",
                   masks="MaskCollection",
                   rowMasks="integer"))

## MultipleAlignmentSet constructor
MultipleAlignmentSet <- function(x=character(), start=1, end=nchar(x[[1]]),use.names=TRUE)
{
  ## Test that all lengths are == to each other
  if(length(x)>0){
   if(!all(unlist(lapply(x,nchar)) %in% nchar(x[[1]]))){
      stop("All the Strings in an MAS Set have to be the same length")}
  }
  new("MultipleAlignmentSet", set=BStringSet(x=x,start=start,end=end,width=NA,
                         use.names=use.names),
                       masks=Mask(mask.width=nchar(x[[1]]),start=0,
                         end=0),
                       rowMasks=integer(length=length(x))
      )
}

setMethod("[","MultipleAlignmentSet",
          function(x, i, j, ..., drop){
            if (!missing(j) || length(list(...)) > 0)
              stop("invalid subsetting")
            if (missing(i) || (is.logical(i) && all(i)))
              return(x)
            if (is.logical(i))
              i <- which(i)
            if (!is.numeric(i) || any(is.na(i)))
              stop("invalid subsetting")
            if (any(i < 1) || any(i > length(x)))
              stop("subscript out of bounds")
            x@set <- x@set[i]; x})

setMethod("[[","MultipleAlignmentSet",
          function(x, i, j, ..., value){
            if (!missing(j) || length(list(...)) > 0)
              stop("invalid subsetting")
            if (missing(i) || (is.logical(i) && all(i)))
              return(x)
            if (is.logical(i))
              i <- which(i)
            if (!is.numeric(i) || any(is.na(i)))
              stop("invalid subsetting")
            if (any(i < 1) || any(i > length(x)))
              stop("subscript out of bounds")            
            x@set[[i]]})

setMethod("names","MultipleAlignmentSet",function(x){names(x@set)})
setReplaceMethod("names", "MultipleAlignmentSet",
                 function(x,value){names(x@set)<-value; x})

setMethod("as.list","MultipleAlignmentSet",function(x){as.list(x@set)})
setMethod("as.character","MultipleAlignmentSet",function(x){as.character(x@set)})
setMethod("nchar","MultipleAlignmentSet",function(x){nchar(x@set)})
setMethod("length","MultipleAlignmentSet",function(x){length(x@set)})



######################################################################
## Below here we implement versions of methods that are "mask aware"
######################################################################

## TODO: change narrow so that it does what it says it does (instead of the
## opposite) OR just jettison narrow (we don't really need it here, I am just
## using it for some quick test.
setMethod("narrow","MultipleAlignmentSet",
          function(x,start,end,width,use.names){
            x@masks <- Mask(mask.width=nchar(x[[1]]),start,end,width);x})



## TODO: remove this after we refactor subsetColumns()
setGeneric("xscat", signature="...",
    function(...) standardGeneric("xscat"))
setMethod("xscat","MultipleAlignmentSet",
          function(...){
            ans_set <- do.call(xscat, lapply(list(...), function(x)x@set))
            new("MultipleAlignmentSet", set=ans_set)
          })



## I definitely want to OL these to consider the mask slot
setMethod("consensusMatrix","MultipleAlignmentSet",
          function(x){consensusMatrix(x@set)})
setMethod("consensusString","MultipleAlignmentSet",
          function(x){consensusString(x@set)})
setMethod("alphabetFrequency","MultipleAlignmentSet",
          function(x){alphabetFrequency(x@set)})








## I might need to do something like this (but modded to remove chars instead
## of # them) I need to be able to call this on the Bstrings, BUT, I also only
## need to be able to do it in the context of these local show methods, so
## probably this is just a local function instead of a method...
## setMethod("as.character", "MaskedXString",
##     function(x)
##     {
##         ans <- as.character(unmasked(x))
##         nir0 <- as(x, "NormalIRanges")
##         for (i in seq_len(length(nir0))) {
##             strip <- paste(rep.int("#", width(nir0)[i]), collapse="")
##             substr(ans,  start(nir0)[i], end(nir0)[i]) <- strip
##         }
##         ans = gsub("#","",ans)
##         ans
##     }
## )


.charStrFormat = function(str, msk){
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





## TODO: For the show method, We need to make sure that we display a special
## message when the entire thing is masked...

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

.XStringMAS.show_frame_line <- function(x, i, iW, widthW, index)
{
    width <- nchar(x)[i]
    snippetWidth <- getOption("width") - 2 - iW - widthW
    if (!is.null(names(x)))
        snippetWidth <- snippetWidth - .namesMASW - 1   
    seq_snippet <- SeqSnippetSlice(x[[i]], snippetWidth, index, mask=x@masks)
    
    if (!is.null(names(x)) && !is.null(seq_snippet))
      seq_snippet <- format(seq_snippet, width=snippetWidth)
    cat(format(paste("[", i,"]", sep=""), width=iW, justify="right"), " ",
        " ",
        seq_snippet,
        sep="")
    if (!is.null(names(x)) && !is.null(seq_snippet)) {
      snippet_name <- names(x)[i]
      if (is.na(snippet_name))
        snippet_name <- "<NA>"
      else if (nchar(snippet_name) > .namesMASW)
        snippet_name <- paste(substr(snippet_name, 1, .namesMASW-3),
                              "...", sep="")
      cat(" ", snippet_name, sep="")
    }
    cat("\n")
}


.XStringMAS.show_frame <- function(x, half_nrow=9L)
{
  ##hard codes the number of sets of rows we plan to display
  .frameLines <- 3
  
    lx <- length(x)
    iW <- nchar(as.character(lx)) + 2 # 2 for the brackets
    ncharMax <- max(nchar(x))
    widthW <- max(nchar(ncharMax), nchar("width"))
    .XStringMAS.show_frame_header(iW, widthW, !is.null(names(x)))

  for(index in seq_len(.frameLines)){   
    if (lx <= 2*half_nrow+1) {
      for (i in seq_len(lx))
        .XStringMAS.show_frame_line(x, i, iW, widthW, index)
        cat("\n")
    } else {
      for (i in 1:half_nrow)
        .XStringMAS.show_frame_line(x, i, iW, widthW, index)
      cat(format("...", width=iW, justify="right"),
          format("...", width=widthW, justify="right"),
          "...\n")
      for (i in (lx-half_nrow+1L):lx)
        .XStringMAS.show_frame_line(x, i, iW, widthW, index)
        cat("\n")
    }
  }
}


setMethod("show", "MultipleAlignmentSet",
    function(object)
    {
        cat("  A ", class(object), " instance with ", length(object), " sequences represented.","\n", sep="")
        if (length(object) != 0)
            .XStringMAS.show_frame(object)
    }
)




## function to make this from a clustal file.
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
  MultipleAlignmentSet(data) 
}


## helper function for subsetting the columns
subsetColumns <- function(x, start, end){
  names = names(x)
  if(length(start) != length(end)){
    stop("You must supply an equal number of starts and ends.")
  }
  x@masks <- Mask(mask.width=nchar(x[[1]]),start,end);
  x
}


