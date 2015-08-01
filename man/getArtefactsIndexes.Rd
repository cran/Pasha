\name{getArtefactsIndexes}
\alias{getArtefactsIndexes}
\title{Identify reads implied in piles potentially considered as PCR artefacts}
\description{ This function uses an 'AlignedData' object (containing reads from
a single chromosome) to give several informations about the reads piles
distribution and identify the reads that can be considered as artefactuals.
It is one of the three major functions used in the pipeline.
} 
\usage{ getArtefactsIndexes(alignedDataObject, 
                            expName, 
                            thresholdToUse=1,
                            thresholdForStats=c(1:5,10,20,50,100),
                            resultFolder=".") } 
                            
\arguments{
  \item{alignedDataObject}{An instance of the class AlignedData containing the
  reads. Reads in the object MUST be part of a unique chromosome (use split on
  seqnames if necessary, see example).
  Alternatively, this argument can be an environment containing the object
  (useful for pass-by-reference in case of big objects). In this case, the
  environment must contain the object in a variable named 'value'}

  \item{expName}{An atomic string. A name or ID of the experiment, used for
  creating filenames for figures generated during the processing}
  \item{thresholdToUse}{An atomic positive integer. Defines the threshold that
  will be used to call artefactual regions and identify reads associated to
  it.}

  \item{thresholdForStats}{A positive integer vector. Defines several
  thresholds that will be plotted on graphics together with the
  'thresholdToUse' value.}

  \item{resultFolder}{An atomic string. Path to a valid folder where the
  figures will be created} }

\details{ Briefly three steps are performed by this function.
First, on the selected reads, a distribution of piles (reads with the exact
same coordinates) will be computed, and a summary plot will be created.
Second, based on the thresholds provided, the algorithm will establish the list
of reads that lies in regions where the piles are higher than each threshold.
At this step, a plot will be appended to previous figure and will summarize the
proportion of reads in the object that would be considered as artefacts for each
threshold.
Finally, in paired-end experiments, each read will be examined and checked
whether its mate is also implied in another artefactual pile or not. A plot is
also produced at this step to graphically summarize the number of reads that will
be reported for each threshold.
For paired-ends experiments, this plot will also give information about orphan
(only one of both reads is in an 'artefactual' pile) or paired artefacts.

For single-ends experiments, all reads in an 'artefactual' pile will be
reported as 'to be removed' except one.
For paired-ends experiments however all reads correponding to a pile will be
removed since the previous strategy would require an arbitrary pair to be
selected among the reads in the pile.
}

\value{
This function returns a 3-level nested list.

First level correspond to the strand of the reads since both strands are
treated separately.

Second level describe the different categories of information provided. It
includes "indexesToRemoveFromTotal" and
"indexesToRemoveFromTotal_isPairedArtefact".

'indexesToRemoveFromTotal', is an integer vector of indexes pointing to the
elements in the object that are considered as "reads implied in an artefactual
pile".

'indexesToRemoveFromTotal_isPairedArtefact' is however a logical vector, each
element corresponding to an element of 'indexesToRemoveFromTotal' (a read) and
specify if this read had its mate in an artefactual pile too.
Information returned only concerns the threshold defined in 'thresholdToUse'.
}

\author{Romain Fenouil}
\seealso{
  \code{\link[Pasha]{AlignedData-class}}
  \code{\link[Pasha]{processPipeline}}
  \code{\link[S4Vectors]{Rle}}
}

\examples{
# Get the path to the example BAM file (and the index)
exampleBAM_fileName <- system.file("extdata",
                                   "embedDataTest.bam",
                                   package="Pasha")
                                
                                

# Reading aligned file
myData <- readAlignedData(folderName="", 
                          fileName=exampleBAM_fileName, 
                          fileType="BAM",
                          pairedEnds=FALSE)

# Split the dataset in a list by chromosomes (not mandatory in this case
# but done like this in the pipeline for consistency with others functions
# which need objects with single chromosomes)
myData_ChrList <- split(myData, seqnames(myData))

# Get the indexes of reads in artefactual piles for each chromosome
indexesArtefacts <- sapply(myData_ChrList, 
                           getArtefactsIndexes, 
                           expName="myExperiment",
                           thresholdToUse=1, 
                           thresholdForStats=c(1:5), 
                           resultFolder=tempdir())
}
