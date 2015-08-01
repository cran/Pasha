\name{mergeWigs}
\alias{mergeWigs}
\title{Merge several different wig files in a unique result file}
\description{ This function computes the average for each genomic position in
several WIG files (fixed step).
} 

\usage{ mergeWigs(files, binSize, outputFolder) } 
\arguments{
  \item{files}{A named list of character vectors, each one containing the paths
  to the wig files that will be merged with other ones, one can specify several
  elements (as much as expected merge results) that will be treated
  sequentially}
  \item{binSize}{The size of the bins (steps) in the wig files, must be the
  same for all files}
  \item{outputFolder}{The folder where the results will be written, created
  recursively if doesn't exist} } 

\details{ For each element of the files list, this function reads all the wigs
files specified and compute an average of all wigs for each position of each
chromosome.
This function only works for fixedStep wigs and is based on the unflexible
format defined by readWIG and writeWIG functions. Thus it might not be compliant
with all the WIG fixed step specification.
The result is directly written to disk as wig in the specified folder.
} 

\value{ NULL }
\author{Romain Fenouil}
\seealso{
  \code{\link[Pasha]{readWIG}}
  \code{\link[Pasha]{writeWIG}}
} 

\examples{
# Generate artificial wig fixed step files with random counts
exampleFolder <- tempdir(); nbFiles=6; nbChromosomes=4

for(fileIndex in 1:nbFiles)
{
    wigScoresContent <- lapply(as.list(1:nbChromosomes), 
        function(x) {sample(x=1:20, size=1000, replace=TRUE, prob=NULL)})
    names(wigScoresContent)<- paste("chr", 1:nbChromosomes, sep="");
    writeWIG(wigScoresContent, 
                paste("wigFileExample", fileIndex, sep=""),
                folder=exampleFolder, fixedStep=200) 
}

# Build a file list containing wig files to merge together, and a
# reference to the file to be used as input for score subtraction
wigFileList <- paste(exampleFolder, 
                     "/",
                     "wigFileExample",
                     1:nbFiles,
                     ".wig", sep="")

# Merge files together (average the scores at each coordinates),
# producing two WIG files
filesToMerge <- list()
filesToMerge[["MyEXP1"]] <- wigFileList[1:(nbFiles/2)]
filesToMerge[["MyEXP2"]] <- wigFileList[((nbFiles/2)+1):nbFiles]

mergeWigs(files=filesToMerge, binSize=200, outputFolder=exampleFolder)
}
\keyword{ merge }
