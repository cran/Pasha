\name{pileupDouble}
\alias{pileupDouble}
\alias{C_pileupDouble}
\title{Piling vector of double}

\description{
This function is inspired from ShortRead deprecated 'pileup' function (Martin
Morgan & Simon Anders). Deprecated function usage is now tipically replaced by
'coverage' from GenomicRanges.
A readaptation of this strategy had to be done in order to allow double
values as weight (particularly useful for multireads).

Briefly, this function computes a numeric vector of pileup scores (coverage)
emcompassing all the coordinates covered by blocks defined by arguments 'start',
'fragLenth', 'dir', and 'readLength'.
}

\usage{
pileupDouble(start, fragLength, dir, readLength, weight)
}

\arguments{
  \item{start}{Coordinates of leftmost end of reads on the chromosome}
  \item{fragLength}{Size of elongation in bp (DNA fragment length)}
  \item{dir}{Strand on which the read has been aligned on the reference
  sequence. IMPORTANT : see remark in description}
  \item{readLength}{Length of the sequenced read in bp}
  \item{weight}{Double value specifying the weight for each read (can be useful
  to put it <1 in some cases for reads aligned in multiple places)} 
}

\details{
Arguments 'fragLength' 'dir' 'readLength' and 'weight' must have the same
length than 'start'.
If one of these argument has length==1, the value will be repeated for all reads.

This function re-encodes internally the strand (dir) to a factor with levels
c('+', '-') to ensure consistency in results. All other values will be ignored
with a warning.

Internally, the function iterates on each defined block/read. 
For blocks on positive strand ('dir'=='+'), the score incrementing region is
defined from the 'start' coordinate to 'start'+'fragLength'.
For blocks on negative strand ('dir'=='-'), the score incrementing region is 
defined from 'start'+'readLength'-'fragLength' to 'start'+'readLength' (it is
important to keep in mind that 'start' always represent the leftmost coordinate
of the read, whatever the strand is.) 


Adaptations were made as compared to ShortRead (Martin Morgan & Simon Anders)
pileup function :

\describe{
    \item{}{R wrapper is in charge of generating vectors in memory}
    \item{}{R wrapper converts the data types to C compatible formats (avoids
    R's conversion code in C function)}
    \item{}{R wrapper ensures the arguments have the appropriate length }
    \item{}{Resulting values are written in the vector passed as argument
    (pointer) to C (previously created in the wrapper)}
    \item{}{The chromosome length is not needed anymore, the function returns a
    vector large enough to cover all reads}
    \item{}{If a read extension generates out of bound coordinates, a warning
    message appears but the computation goes on}
    }
}

\value{
A (large) numeric vector.
}

\author{Romain Fenouil}
\seealso{
  \code{\link[IRanges]{coverage}}
}

\examples{
# generate artificial reads coordinates  on a chromosome
nbReads <- 1000000
chrLength <- 10000000

startPositions <- sample(1:chrLength, nbReads, replace=TRUE)
fragmentsLength <- trunc(rnorm(n=nbReads, mean=146, sd= 40))
strandAligned <- factor(sample(c('+','-'), 
                        nbReads,
                        replace=TRUE),
                        levels=c('+', '-'))
readLength <- 36
weight <- 1

# compute the piled vector
res <- pileupDouble(start=startPositions, 
                    fragLength=fragmentsLength,
                    dir=strandAligned, 
                    readLength=readLength, 
                    weight=weight)

# plot distribution of scores on the chromosome
hist(res)
}
