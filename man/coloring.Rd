\name{coloring}

\alias{update_X_palette}
\alias{update_DNA_palette}
\alias{update_RNA_palette}
\alias{update_AA_palette}
\alias{update_B_palette}

\title{Update XString Display Colors}

\description{
  Functions to set custom user color palettes for displaying \code{XString} objects.
}

\usage{
update_DNA_palette(colors=NULL)
update_RNA_palette(colors=NULL)
update_AA_palette(colors=NULL)
update_B_palette(colors=NULL)
}

\arguments{
  \item{colors}{
    A named list of colors to update, with entries \code{fg} and \code{bg} specifying the foreground and background colors, respectively. Colors can be specified in any way compatible with \code{\link[crayon]{make_style}} from the \code{crayon} package. Defaults to \code{NULL}, which resets the color palette to the default color scheme. See Details and Examples for more information.
  }
}

\details{
  These functions allow users to change the default color scheme of Biostrings. Each function expects a \code{list} with named entries corresponding to the values to update. Each entry can specify \code{'fg'} and \code{'bg'} values, corresponding to the foreground and background colors (respectively). If \code{'fg'} is not specified, it defaults to \code{rgb(1,1,1)} (white). If \code{'bg'} is not specified, it defaults to transparent.

  These functions will only update the values passed, leaving the rest of the colors as-is. For example, calling \code{update_AA_palette(list(A=list(fg="green")))} would update the coloring for \code{A} while leaving all other colors as the default schema.

  To reset all colors to the default palette, call the function with no arguments (\code{NULL}).

  To remove a coloring for a specific value, provide a named entry with value \code{NULL}. For example, \code{update_AA_palette(list(A=NULL))} will remove the coloring for \code{A}.

  \code{update_DNA_palette} and \code{update_RNA_palette} are identical internally, so either function can be used to update colorings for \code{T,U}.

  See the Examples section for more examples of custom colorings.
}

\value{
  Invisibly returns the new color mapping, consisting of a named character vector. Calling \code{cat} on the return value will print out all letters with their respective coloring.
}

\author{Aidan Lakshman <AHL27@pitt.edu>}

\seealso{
  \link{XString-class}
}

\examples{
## display default colors
DNAString(paste(DNA_ALPHABET, collapse=''))
AAString(paste(AA_ALPHABET, collapse=''))
BString(paste(LETTERS, collapse=''))

## create new palettes
DNA_palette <- list(
  A=list(fg="blue",bg="black"),
  T=list(fg="red",bg='black'),
  G=list(fg='green',bg='black'),
  C=list(fg='yellow',bg='black')
)
update_DNA_palette(DNA_palette)
DNAString(paste(DNA_ALPHABET, collapse=''))

## reset to default palette
update_DNA_palette()
DNAString(paste(DNA_ALPHABET, collapse=''))

## colors can also be specified with `rgb()`
AA_palette <- list(
  A=list(fg="white", bg="purple"),
  B=list(fg=rgb(1,1,1), bg='orange')
)
update_AA_palette(AA_palette)
AAString(paste(AA_ALPHABET, collapse=''))

## remove all coloring for QEG
update_AA_palette(list(Q=NULL, E=NULL, G=NULL))
AAString(paste(AA_ALPHABET, collapse=''))

## reset to default
update_AA_palette()
AAString(paste(AA_ALPHABET, collapse=''))

## We can also add colors to BStrings,
## which are normally not colored

## if 'fg' is not specified, defaults to rgb(1,1,1)
## if 'bg' is not specified, background is transparent
B_palette <- list(
  A=list(bg='green'),
  B=list(bg="red"),
  C=list(bg='blue'),
  D=list(fg="orange"),
  E=list(fg="yellow")
)
update_B_palette(B_palette)
BString(paste(LETTERS, collapse=''))

## can also directly view the changes with cat
cat(update_B_palette(B_palette), '\n')

## reset to default
update_B_palette()
BString(paste(LETTERS, collapse=''))
}
