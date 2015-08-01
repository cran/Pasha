\name{multiread_CSEMDispatch}
\alias{multiread_CSEMDispatch}
\alias{CSEMDispatch}
\title{Multiread scoring CSEM dispatch}

\description{ This script aims to dispatch the scoring of multi-read aligned
reads according to the CSEM algorithm developped by Chung et al. (see
"Discovering Transcription Factor Binding Sites of Genomes with Multi-Read
Analysis of ChIP-seq Data" (2011) PLoS Computational Biology).
} 

\usage{ multiread_CSEMDispatch(alignedFile, 
                               outputFolder, 
                               referenceFile,
                               window_size=101, 
                               iteration_number=200, 
                               incrArtefactThrEvery=NA, 
                               verbosity=0) 
      }

\arguments{
  \item{alignedFile}{ An atomic character string. The full path to the file
  containing the reads aligned by bowtie with the --concise option. }
  \item{outputFolder}{ An atomic character string. The path to the folder where
  the file output by the script must be stored. }
  \item{referenceFile}{ An atomic character string. Either a full path to a
  reference file (see details for format specification), or the ID of one
  reference included in the package (see details for available ones). }
  \item{window_size}{ A positive integer. The size of the window used by the
  algorithm (see algorithm details). Default value is 101.}
  \item{iteration_number}{ A positive integer. The number of iteration executed
  by the algorithm (see algorithm details). Default value is 200. }
  \item{incrArtefactThrEvery}{ A complex parameter (see details). A numeric
  value or NA. A strictly positive numeric value activate the option that allow
  to remove the 'artifacts', defining a threshold to consider piles like
  'artifacts' as 'number of reads in the experiment devided by
  incrArtefactThrEvery'. A NA will ignore the eventual artifactual piles.}
  \item{verbosity}{ An integer. The verbose level : 0 = no message, 1 = trace
  level } } 

\details{ The script consider the reads that have been aligned in several
location by bowtie (multi-reads). At each read, it assign a score determined by
the CSEM algorithm (Chung et al. "Discovering Transcription Factor Binding Sites
of Genomes with Multi-Read Analysis of ChIP-seq Data" (2011) PLoS Computational
Biology).
The script output a tab separated value text file formated as below:
\itemize{
  \item{Column 1 : Chromosome name}
  \item{Column 2 : Strand}
  \item{Column 3 : Position}
  \item{Column 4 : Score}
  } 

The output file is named according to the name of the input file suffixed with
"_csemDispatch.txt".

If the parameter 'incrArtefactThrEvery' is set to a strictly positive value,
the input file is first passed to a script removing the 'artefacts'. This
script looks at the reads aligned at the exact same position.
If the number of such reads exceed a limit, these reads are considered as due to
experiment artifact.
In that case, only one read is kept and the rest is removed.
The limit defining artefact is controled by the parameter
'incrArtefactThrEvery' such as:

limit = Total number of reads / incrArtefactThrEvery
}

\value{
A tab separated value text file formated as below:
\itemize{
  \item{Column 1 : Chromosome name}
  \item{Column 2 : Strand}
  \item{Column 3 : Position}
  \item{Column 4 : Score}
  }
}

\author{Lionel Spinelli}
\seealso{
  \code{\link[Pasha]{processPipeline}}
  \code{\link[Pasha]{multiread_RemoveArtifact}} 
  \code{\link[Pasha]{multiread_UniformDispatch}}
}

\examples{
# Define input aligned file
my_aligned_file <- system.file("extdata",
                               "embededDataTest_MultiSignal.bow",
                               package="Pasha")

# Define the output folder
my_output_folder <- tempdir()

# Define the genome reference file
genome_reference_file <- system.file("resources",
                                     "mm9.ref",
                                     package="Pasha")

# Launch the script
multiread_CSEMDispatch(my_aligned_file, 
                       my_output_folder, 
                       genome_reference_file,
                       incrArtefactThrEvery=7000000, 
                       verbosity=1) 
}
