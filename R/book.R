setClass("book",
         representation(title="character",
	                author="character",
			chapters="list"))
# How do I specify I want chapters to be a matrix of mode "character"?

MyBook <- function()
{
	b <- new("book")
	b@title <- "Hello world"
	b@author <- "bibi"
	b@chapters <- list(
		list(title="Hello", text=c("aaaa","eeeeeeeee"), start=1, end=24),
		list(title="world", text=c("zzz","ssssss","iii"), start=25, end=39)
		)
	b
}

setMethod("length",
          signature(x = "book"),
          function (x)
      {
          length(x@chapters)
      })