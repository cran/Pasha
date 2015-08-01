\name{multiread_RemoveArtifact}
\alias{multiread_RemoveArtifact}
\alias{RemoveArtifact}
\title{Artifact Remover}
\description{ This script aims to remove the reads issued from experiement
artifact detected in a bowtie aligned read file.} 
\usage{
multiread_RemoveArtifact(alignedFile, 
                         outputPath, 
                         referenceFile,
                         incrArtefactThrEvery=10000000, 
                         verbosity=0) 
      } 

\arguments{
  \item{alignedFile}{ An atomic character string. The full path to the file
  containing the reads aligned by bowtie with the --concise option. }
  \item{outputPath}{ An atomic character string. Path where the resulting file
  will be written. }
  \item{referenceFile}{ An atomic character string. Either a full path to a
  reference file (see details for format specification), or the ID of one
  reference included in the package (see details for available ones). }
  \item{incrArtefactThrEvery}{ A complex parameter (see details). A numeric
  value or NA. A positive numeric value defines a threshold to consider piles
  like 'artefacts' as 'number of reads in the experiment devided by
  incrArtefactThrEvery'. A negative value allows to directly set the threshold
  manually, independently of the number of reads. A NA will ignore the eventual
  artefactual piles. }
  \item{verbosity}{ An integer. The verbose level : 0 = no message, 1 = trace
  level } } 

\details{ The script looks at the reads aligned at the exact same position.
If the number of such reads exceed a limit, these reads are considered as due to
experiment artefact. In that case, only one read is kept and the rest is
removed.
The limit defining artefact is controled by the parameter 'incrArtefactThrEvery'
such as:

limit = Total number of reads / incrArtefactThrEvery

The script returns a file formatted as the input file (text file formated as
the bowtie output with --concise option).
The file is named according to the name of the input file suffixed with
"_noA.bow".
}

\value{
A text file formated as the bowtie output with --concise option.
}

\author{Lionel Spinelli}
\seealso{
  \code{\link[Pasha]{processPipeline}}
  \code{\link[Pasha]{multiread_UniformDispatch}}
  \code{\link[Pasha]{multiread_CSEMDispatch}} 
}

\examples{
# Define input aligned file
# Define input aligned file
my_aligned_file <- system.file("extdata",
                               "embededDataTest_MultiSignal.bow",
                               package="Pasha")

# Define the genome reference file
genome_reference_file <- system.file("resources", "mm9.ref",package="Pasha")

# Define the output folder
my_output_folder <- tempdir()

# Launch the script
multiread_RemoveArtifact(my_aligned_file,
                         my_output_folder,
                         genome_reference_file,
                         7000000,
                         0) 
}
