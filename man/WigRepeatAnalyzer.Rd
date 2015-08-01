\name{WigRepeatAnalyzer}
\alias{WigRepeatAnalyzer}
\alias{AnalyzeRepeat}
\title{WIG Repeat Analyzer}

\description{
This script aims at analyzing the distribution (in coverage and score) of a
signal (stored in a WIG file) according to genome repeat annotations.
}

\usage{
WigRepeatAnalyzer(filename, 
                  inputFolder, 
                  outputFolder="./",
                  repeatMaskerFilePath,
                  isRegex=FALSE) 
}

\arguments{
  \item{filename}{ An atomic character string. The file name of the wig file
  (fixed step WIG). }
  \item{inputFolder}{ An atomic character string. The path to the wig file
  (fixed step WIG). }
  \item{outputFolder}{ An atomic character string. The path to the folder where
  analysis results must be stored. }
  \item{repeatMaskerFilePath}{ An atomic character string. The path to the file
  containing the repeat annotations (Repeat Masker file). }
  \item{isRegex}{ An atomic logical. If TRUE, the filename parameter is
  interpreted as a regular expression (LC_SYNTAX) and the script will search
  for a unique file corresponding to the provided regular expression. If no or
  several file are found, the scripts ends with error. } 
} 
\details{ The script consider the part of the signal contained in the WIG file
that fall into repeat annotations provided by user.

This annotation file must formatted as the Repeat Masker annotations file
provided by UCSC. The required information are:

\itemize{
  \item{Column 6 : Chromosome name}
  \item{Column 7 : Start position}
  \item{Column 8 : End position  }
  \item{Column 12 : Repeat Class}
  \item{Column 13 : Repeat Family}
}


The script compute the coverage of each repeat class and family (i.e. the
percentage of positions falling into each annotations) and the weight of each
class and family (i.e. the percentage of score falling into each annotations).
All results are provided as barplots figures and text files.
}
\value{
Barplot figures and text files of raw data for coverage and weight of repeat
classes and families.
}

\author{Lionel Spinelli}
\seealso{
  \code{\link[Pasha]{processPipeline}}
}

\examples{
# Build temp dir
exampleFolder <- tempdir()

# Generate artificial wig fixed step files with random counts
nbChromosomes <- 1
bin_size <- 200
wigScoresContent <- lapply(as.list(1:nbChromosomes), 
     function(y) {sample(x=1:20,
                         size=floor(4333587/bin_size), 
                         replace=TRUE,
                         prob=NULL)})
names(wigScoresContent) <- 1:nbChromosomes

writeWIG(wigScoresContent,
         "wigFileExample",
         folder=exampleFolder, 
         fixedStep=bin_size)

# Define input file
my_wig_file <- "wigFileExample.wig"

# Define the RepeatMasker file
my_repeat_masker_file <- system.file("extdata",
                                     "mm9_RepeatMasker_2000l.txt",
                                     package="Pasha")

\dontrun{
# Launch the script
result <- WigRepeatAnalyzer(filename=my_wig_file, 
                            inputFolder=exampleFolder,
                            outputFolder=exampleFolder, 
                            repeatMaskerFilePath=my_repeat_masker_file)
} 
}
