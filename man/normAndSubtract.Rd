\name{normAndSubtractWIG}
\alias{normAndSubtractWIG}
\title{Normalization and subtraction of WIG files}
\description{ Function that linearly rescales scores in provided WIG files and
subtract Input WIG file if asked } 

\usage{ 
normAndSubtractWIG(wigFileList,
                   inputFile=NA, 
                   rescaleInput=FALSE,
                   meanBGLevelRescale=c(0,10), 
                   subtractInput=TRUE, 
                   binSize=50) 
      } 

\arguments{
  \item{wigFileList}{A character vector, each element describing the full path
  to the wig file that has to be treated, when several elements, each one will
  be treated sequentially}
  \item{inputFile}{An atomic character vector describing the full path to a wig
  file from which the genomic scores will eventually be subtracted to
  wigFileList wig files scores, NA otherwise}
  \item{rescaleInput}{An atomic logical deciding if the eventual inputFile
  scores have to be rescaled (with the same strategy employed for wigFileList
  files) before using them for subtraction}
  \item{meanBGLevelRescale}{A numeric vector of size 2, the first number is the
  lower limit of the resulting scores after scaling, the second number
  represent the score that will represent the average enrichment on the genome
  after rescaling, all scores being rescales linearly around this value}
  \item{subtractInput}{An atomic logical specifyijng that the input scores
  (eventually rescaled) from inputFile, must be subtracted to resulting files
  or not}
  \item{binSize}{An atomic integer value specifying the step size in the wig
  files} } 

\details{ For each wig file in the character vector, this functions reads all
the scores for each chromosome, rescale the range of scores, eventually subtract
the scores from an input experiment, and saves the results in the same folder as
the original file.
The rescaling linearly spreads the scores from the lower score to the average
score (of the whole original file) to respectively the first and second element
of meanBGLevelRescale argument.
If an input file is specified, the scores stored for each genomic location are
eventually (depending on the arguments) subtracted from corresponding scores in
the original wigs.
If the user decides to rescale the input before subtraction (recommended), the
rescaled inputFile will also be saved in the the original file folder with the
_Rescaled extension and can then further be used directly for further
computation, with no need to rescale it again (gain of time).
} 

\value{ Writes the resulting rescaled and input subtracted files in the same
folders as originals, with suffix in filenames _rescaled and _BGsub depending on
the parameters }

\author{Romain Fenouil}
\seealso{
  \code{\link[Pasha]{readWIG}}
  \code{\link[Pasha]{writeWIG}}
  \code{\link[Pasha]{mergeWigs}} 
} 

\examples{
# Generate artificial wig fixed step files with random counts
exampleFolder <- tempdir() 
nbFiles <- 3
nbChromosomes <- 4

for(fileIndex in 1:nbFiles)
{
    wigScoresContent <- lapply(as.list(1:nbChromosomes), 
        function(x) {sample(x=1:20, size=1000, replace=TRUE, prob=NULL)})
    names(wigScoresContent) <- paste("chr", 1:nbChromosomes, sep="")
    writeWIG(wigScoresContent, 
             paste("wigFileExample", fileIndex, sep=""),
             folder=exampleFolder, 
             fixedStep=200)
}

# Build a file list containing wig files to merge together, and a reference
# to the file to be used as input for score subtraction
wigFileList <- file.path(exampleFolder, 
                         paste("wigFileExample",
                               1:(nbFiles-1), 
                               ".wig", 
                               sep=""))
                                
inputFile <- file.path(exampleFolder, 
                       paste("wigFileExample", 
                             nbFiles,
                             ".wig", 
                             sep=""))

# Run the merge function on 'wigFileList' which will rescales the scores and
# then subtract the wig scores used as input experiment (rescaled too). 
normAndSubtractWIG(wigFileList, 
                   inputFile, 
                   rescaleInput=TRUE,
                   meanBGLevelRescale=c(0,10), 
                   subtractInput=TRUE, 
                   binSize=200)
}
\keyword{ normalize }
