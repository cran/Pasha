\name{estimateElongationSize}
\alias{estimateElongationSize}
\title{Insert size estimation in-silico}

\description{ This function tries to infer the sequenced original DNA fragments
size from the reads positions on both strands. It is one of the three main
functions used in the pipeline.
} 

\usage{ estimateElongationSize(alignedDataObject, 
                               expName, 
                               stepMin=50,
                               stepMax=450, 
                               stepBy=10, 
                               averageReadSize=NA, 
                               resultFolder=".") 
      } 

\arguments{
  \item{alignedDataObject}{An instance of the class AlignedData containing the
   reads. Typically reads are selected to be part of a single chromosome. 
   Alternatively, this argument can be an environment containing the object 
   (useful for pass-by-reference in case of big objects). In this case, the 
   environment must contain the object in a variable named 'value'}

  \item{expName}{An atomic string. A name or ID of the experiment, used for 
  creating filenames for figures generated during the estimation}

  \item{stepMin}{An atomic strictly positive integer. Defines the minimum DNA 
  fragment size assumed to be sequenced, the algorithm will start evaluation 
  from this value up to 'stepMax'}

  \item{stepMax}{An atomic strictly positive integer. Defines the maximum DNA 
  fragment size assumed to be sequenced, the algorithm will estimate the 
  fragment size up to this value}

  \item{stepBy}{An atomic strictly positive integer. Represents the resolution 
  (in basepairs) of the DNA fragment size estimation. Note that the 
  computation time is proportional to the number of values to be tested (ie. 
  (stepMax-stepMin)/stepBy)}

  \item{averageReadSize}{An atomic strictly positive integer. For computation
  purpose, this function needs to know the reads length. If not available in
  the object (it's an optional field in the class representation), the user
  have to specify an average size as argument of the function}

  \item{resultFolder}{An atomic string. Path to a valid folder where the
  figures will be created} } 

\details{ The goal of this function is, for single-end experiments, to
estimate the average size of DNA fragments that were sequenced (as opposed to
reads size which only represents one of the ends of fragments).
In brief, the function computes separately the piled-up score of both strands
and estimate the shifting between enriched regions.
As opposed to the classic cross-correlation typically used for this concern,
this function is highly sensistive to enrichment and thus favors the estimation
for most enriched regions.
This function lies on a critic code written in C, see 'elongationEstimation' for
more details.
} 

\value{ A graphic per chromosome is plotted to illustrate the score for each
shifting. Based on the local maxima and the reads size, a decision is made about
the size of sequenced DNA fragments.

Returns a numeric being the most probable elongation size for the provided
chromosomes in 'alignedDataObject' Note that this function can also deal with an
object describing reads for several chromosomes. In this case the function will
return a list of integers, each element named as the concerned chromosome.
}
\author{Romain Fenouil}
\seealso{
  \code{\link[Pasha]{AlignedData-class}}
  \code{\link[Pasha]{processPipeline}}
}

\examples{
# Get the path to the example BAM file (and the index)
exampleBAM_fileName <- system.file("extdata",
                                   "embedDataTest.bam",
                                   package="Pasha")

# Reading aligned file 
myData=readAlignedData(folderName="",
                       fileName <- exampleBAM_fileName,
                       fileType="BAM",
                       pairedEnds=FALSE)

# Split the dataset in a list by chromosomes (not mandatory in this case but
# done like this in the pipeline for consistency with others functions which
# need objects with single chromosomes) 
myData_ChrList <- split(myData, seqnames(myData))

# Compute the fragment size estimation on all chromosomes
readsSize <- sapply(myData_ChrList, estimateElongationSize, "myExperiment",
resultFolder=tempdir())
}
