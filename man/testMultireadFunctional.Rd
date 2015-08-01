\name{.testFunctionalMultiread}
\alias{.testFunctionalMultiread}
\alias{.testFunctionalMultiread_generateResults}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions dedicated to the functional validation of multiread feature of Pasha package.}

\description{
The goal of these two functions is to generate a reference dataset with the
currently installed Pasha package in the user's environment and/or to compare
it to the precomputed results in the developpers environment.}

\usage{
.testFunctionalMultiread(folderOutput=tempdir(), 
                         testType="regular",
                         verbose=FALSE) 

.testFunctionalMultiread_generateResults(testFileBOW_folderName,
                                         testFileBOW_fileName,
                                         folderOutputGenerated, 
                                         verbosity=0) 
      }

\arguments{
  \item{folderOutput}{An existing folder (this is tested and can eventually
  raise an error) in which the generated results will be written and the
  precomputed results will be inflated.}
  \item{testFileBOW_folderName}{The folder in which the .bow file for results
  generation stands (.testFunctionalMultiread passes the path to the
  package-embedded test .bow file)}
  \item{testFileBOW_fileName}{The filename of the ?bow file used to generate
  the results (.testFunctionalMultiread passes the name of the
  package-embedded test .bow file)}
  \item{folderOutputGenerated}{The folder in which the intermediary generated
  data files stands (.testFunctionalMultiread passes this path)}
  \item{testType}{A value in c("regular", "complete"). A regular test will
  generate multiread specific results but skips files for which comparison can
  be problematic because of writing floating points with C++ routines (which
  formatting depends on the machine/platform floating numbers implementation).
  Because this test function is called for package validation, a regular test
  exclude these files which would otherwise prevent the package validation. A
  complete test can still be done manually for testing purposes.}
  \item{verbose}{Logical. If verbose is FALSE, the test stdout will be
  redirected to a log file to avoid overflow of test messages.}
  \item{verbosity}{ An integer. The verbose level : 0 = no message, 1 = trace
  level } 
}

\details{
Once installed, Pasha gives access to a functional test of features linked to
multiread management with the function '.testFunctionalMultiread'.

The package includes a .bow file (obtained throufh bowtie with --concise
option) designed to cover most of the cases for multiread datasets (see
publication for details).
This test consists in using the provided R funcntions and the embedded .bow
file to generate results and compare them to a previously generated
(precomputed) dataset that has been functionally validated by developpers.
The precomputed dataset is embedded in the package in a zip file and contains
all the .bow and .txt files that could be necessary for validation. This
dataset is generated and compiled by package developpers using
.testFunctionalMultiread_generateResults.
}

\value{
.testFunctionalMultiread prints to stdout a summary of file generation and
comparison with reference dataset.
}

\author{Lionel Spinelli}
\seealso{
  \code{\link[Pasha]{processPipeline}}
 }