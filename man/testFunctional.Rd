\name{.testFunctional}
\alias{.testFunctional}
\alias{.testFunctional_generateResults}
\title{Functions dedicated to the functional validation of Pasha package.
} 
\description{ The goal of these two functions is to generate a reference
dataset with the currently installed Pasha package in the user's environment
and/or to compare it to the precomputed results in the developpers environment.
} 
\usage{ 
.testFunctional(folderOutput=tempdir(), testType="regular",
verbose=FALSE) 

.testFunctional_generateResults(testFileBAM_folderName,
                                testFileBAM_fileName, 
                                testFileMulti_fileName,
                                testType="complete") 
      }

\arguments{
  \item{folderOutput}{An existing folder (this is tested and can eventually
  raise an error) in which the generated results will be written and the
  precomputed results will be inflated.}
  \item{testFileBAM_folderName}{The folder in which the BAM file for results
  generation stands (.testFunctional passes the path to the
  package-embedded test BAM file)}
  \item{testFileBAM_fileName}{The filename of the BAM file used to generate the
  results (.testFunctional passes the name of the package-embedded test
  BAM file)}
  \item{testFileMulti_fileName}{The filename of the file that is used as test
  to generate ouptut results for comparison. It simulates reads that were
  treated with the multiread processing functions before injecting it in the
  pipeline.}
  \item{testType}{A value in c("regular", "complete"). A regular test will
  generate ~15 WIG and GFF files, a complete test however generates ~750 result
  files from the test BAM file and aims at passing in all possible loops and
  blocks of the pipeline.}
  \item{verbose}{Logical. If verbose is FALSE, the pipeline stdout will be
  redirected to a log file to avoid overflow of test messages.} 
} 

\details{
Once installed, Pasha gives access to a functional test with the function
'.testFunctional'.

The package includes a SAM/BAM file designed to cover most of the cases for
paired-ends or single-ends datasets (see publication for details).
This test consists in using the provided pipeline and the embedded BAM file to
generate results and compare them to a previously generated (precomputed)
dataset that has been functionally validated by developpers.
The precomputed dataset is embedded in the package in a zip file and contains
all the WIG and GFF files that could be necessary for validation. This dataset
is generated and compiled by package developpers using
.testFunctional_generateResults.
Two kind of tests can be computed:
A regular test will go through the most common options for single-end and
paired-ends, it will generate 15 to 20 results file and compare them to the
references, the generation and comparisons should take less than a minute to
complete.
A complete test will launch the generation of results by using as much as
possible different options and parameters to the pipeline function. It will
typically generate ~1000 wig and gff files and compare them to the validated
reference, this test can be long (15 to 30 minutes) but guarantees to obtain the
same results as the developpers for all kind of analysis.
} 

\value{ .testFunctional prints to stdout a summary of file generation and
comparison with reference dataset.
}

\author{Romain Fenouil}
\seealso{
  \code{\link[Pasha]{AlignedData-class}}
  \code{\link[Pasha]{processPipeline}}
 }

