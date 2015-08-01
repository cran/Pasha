\name{generatePiled}
\alias{generatePiled}
\title{Pileup reads in a numeric vector with several piling and elongation
options} 

\description{
This function converts the reads information (position, length, and eventually
weight) from an 'AlignedData' object into enrichment scores along the
concerned region (represented as a numeric vector).
Its behaviour depends on several arguments described below (see also pdf
document attached to the package). In brief, it is able to handle single or
paired -ends reads, classic or midpoint piling, and elongation or not,
depending on the users needs.
It is one of the three major functions used in the pipeline.
}

\usage{
generatePiled(alignedDataObject, 
              elongationSize, 
              averageReadSize, 
              midPoint=FALSE)
}

\arguments{
  \item{alignedDataObject}{An instance of the class AlignedData containing the
  reads. Typically reads are selected to be part of a single chromosome
  (seqnames).
  Alternatively, this argument can be an environment containing the object
  (useful for pass-by-reference in case of big objects). In this case, the
  environment must contain the object in a variable named 'value'}
  \item{elongationSize}{A positive numeric value or NA. See details}
  \item{averageReadSize}{An atomic strictly positive integer. For computation
  purpose, this function needs to know the reads length. If not available in
  the object (qwidth is an optional field in the class representation), the
  user have to specify an average size as argument of the function}
  \item{midPoint}{An atomic logical. Defines the type of piling that will be
  applied to reads. See details} }

\details{
This function will produce very different results depending on the parameters
 provided by the user (elongationSize, midpoint) and the type of experiment
(paired-ends or single-ends).

First, the specified elongationSize can induce different piling strategies. 

The value '0' is special and interpreted as 'no elongation'. The piling result
will only represent the sequenced reads region. Internally the elongation value
is replaced by each read size.
If the user asks for midpoint piling strategy with elongationSize 0, the result
will only represent scores corresponding to the first leftmost or rightmost
base for reads respectively on positive and negative strand.

A strictly positive value is only valid for single-ends experiments
(paired-ends objects are assumed to carry their own insert size based on the
difference of coordinates between a read and its mate, hence NA is expected for
paired-ends experiments, if not, a warning is raised and the value is ignored).
This value represents the average DNA fragments (inserts) size in the
experiment (eventually computed with estimateElongationSize). For classic
elongation process (no midpoint) this value will be used to elongate the read
prior to pileup in order to represent the full DNA fragment in the results.
A value lower than the reads size is technically possible but scientifically
improbable, and can induce misleading interpretations, specially when compared
to elongation 0.
For midpoint piling, only the the first leftmost or rightmost base (for reads
respectively on positive and negative strand) are taken in account during
piling. However the elongationSize (coming from corresponding argument for
singleEnds or from object for paired-ends) is used to shift this score to the
expected 'center' or DNA fragment (ie. DNA fragment size size divided by 2).

As written earlier, NA is expected for paired-ends experiments but NOT for
single-ends ones which need an elongation value. Otherwise an error will stop
the execution.

All these cases are documented and represented as graphics in the pdf file
attached to this package.
}

\value{
Returns a Rle object representing an integer vector, each score corresponding
to the enrichment at a chromosomic coordinate. The vector is long enough to
represent the full region covered by reads in the object.
Rle objects allow to carry metadata and this function stores the chromosome
name in it.
Note that this function can also deal with an object describing reads for
several  (seqnames). In this case the function will return a list of Rle
objects, each element named as the concerned chromosome.
}

\author{Romain Fenouil}
\seealso{
  \code{\link[Pasha]{AlignedData-class}}
  \code{\link[Pasha]{processPipeline}}
  \code{\link[Pasha]{estimateElongationSize}}
  \code{\link[S4Vectors]{Rle}}
}

\examples{
# Get the path to the example BAM file
exampleBAM_fileName <- system.file("extdata",
                                   "embedDataTest.bam",
                                   package="Pasha")

# Reading aligned file
myData <- readAlignedData(folderName="",
                          fileName=exampleBAM_fileName, 
                          fileType="BAM",
                          pairedEnds=FALSE)

# Split the dataset in a list by chromosomes (not mandatory in this case but
# done like this in the pipeline for consistency with others functions which
# need objects with single chromosomes) 
myData_ChrList <- split(myData, seqnames(myData))

# Compute the piling of reads on each chromosome
piledUpVectors <- sapply(myData_ChrList, 
                         generatePiled, 
                         elongationSize=146,
                         averageReadSize=NA,
                         midPoint=TRUE) 
}
