\name{readWIG and writeWIG}
\alias{readWIG}
\alias{writeWIG}
\title{Read and Write genomic WIG (fixed step) files}

\description{
Functions dedicated for reading and writing WIG fixed step from their named
list representation in R.
}

\usage{
readWIG(fileName)
writeWIG(wigData, fileName, folder=NULL, fixedStep=50, addExtension=TRUE)
}

\arguments{
  \item{fileName}{Path to the file on which to write the results or from which
  to read the data}
  \item{wigData}{A named list of numeric vectors representing each chromosome
  (or track) piled data scores to be written}
  \item{folder}{Eventual folder path (for writeWIG) that will be concatenated
  to the filename (can be NULL if the user prefers to specify the full path in
  fileName)}
  \item{fixedStep}{Size of the bins used for the representation of piled data
  in the named list (the function can not detect it since this information is
  not carried by wigData)}
  \item{addExtension}{Single logical. Specify if the extension string "wig"
  should be added to the filename if 'fileName' argument does not end with it
  (case is ignored for search).} 
} 

\details{
IMPORTANT : These functions can only handle fixed step WIGs with a limited
format.
They assume a very UNflexible format with two lines of description for each
track and can't be used in other contexts.
Binsize specified in the file is IGNORED while reading the file and the user is
assumed to know it.
}

\value{
readWIG returns a named list of numeric vectors (each track in the wig gives a
new list element), each numeric value correspond to a bin in the corresponding
track.
}

\author{Romain Fenouil, Nicolas Descostes}
\seealso{
  \code{\link[Pasha]{mergeWigs}}
}

\examples{
## Define an artificial wig list
generateWigScores <- function(scores=1:20, lengthTrack=1000)
{
    return(sample(scores, lengthTrack, replace=TRUE))
}

myWig <- list("chr1"=generateWigScores(), "chr2"=generateWigScores())

# Write a wig file that can be read in genome browsers
writeWIG(myWig, "myWigFixedSteps", folder=tempdir(), fixedStep=200)

# Read the file generated
myWigFromFile <- readWIG(file.path(tempdir(), "myWigFixedSteps.wig"))

# Check that scores are equal for all chromosomes
if(!all(mapply(all.equal,myWig,myWigFromFile))) stop ("Error, the original 
generated data look different from what has been read");
}
\keyword{ wig }
