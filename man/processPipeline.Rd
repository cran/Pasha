\name{processPipeline}
\alias{processPipeline}
\title{A pipeline to accurately transform aligned reads to genomic enrichment
scores}
\description{ 
This function automatizes the processing of chromatin sequecing data analysis 
by chaining all required steps from the reading of aligned reads to the writing
of enrichment scores in WIG/GFF files.

}

\usage{
processPipeline(
        # I/O GENERAL PARAMETERS
        INPUTFilesList,
        resultSubFolder             = "Results_Pasha",
        reportFilesSubFolder        = ifelse(length(resultSubFolder)>1,
                                                resultSubFolder[2],
                                                "ReportFiles"),
        WIGfs                       = TRUE,
        WIGvs                       = FALSE,
        GFF                         = FALSE,
        # COMPLEX PARAMETERS (SINGLE OR VECTORS OR LIST OF IT)
        incrArtefactThrEvery        = 7000000,
        binSize                     = 50,
        elongationSize              = NA,
        rangeSelection              = IRanges(0,-1), 
        annotationFilesGFF          = NA, # GFF files
        annotationGenomeFiles       = NA, # path to file or "mm9", "hg19"... 
        # SINGLE PARAMETERS
        elongationEstimationRange   = c(mini=150, maxi=400, by=10),
        rehabilitationStep          = c("orphans","orphansFromArtefacts"),
        removeChrNamesContaining    = "random|hap",
        ignoreInsertsOver           = 500,
        nbCPUs                      = 1,
        keepTemp                    = TRUE, # Keep the intermediary files
                                            # that led to the final ones
                                            # (rehab and multireads) 
        logTofile                   = "./log.txt", 
        eraseLog                    = FALSE, 
        # LIST PARAMETERS (one element per expName)
        multiLocFilesList           = list()) # A list with experiments
                                              # names and associated
                                              # filenames to treat 
}

\arguments{
  \item{INPUTFilesList}{A named list of list. Each element name will be used as
  experiment ID to generate output files and messages. Each sublist must
  provide the following elements as described (see details for more information
  on individual elements) : 'folderName', 'fileName', 'fileType', 'chrPrefix',
  'chrSuffix', 'pairedEnds', 'midPoint'}
  
  \item{resultSubFolder}{A single character string. Name of the subfolder to
  be created in the folder where each experiment file is. The results files
  will be stored in this subfolder}
  
  \item{reportFilesSubFolder}{A single character string. Name of the
  subfolder in which log files and QC reports will be stored. This name is
  used to create a filder under resultSubFolder}
  
  \item{WIGfs}{A single boolean. If TRUE, a result file in WIG fixed step
  format (see binSize parameter) will be generated}
  
  \item{WIGvs}{A single boolean. If TRUE, a result file in WIG variable step
  format will be generated}
  
  \item{GFF}{A single boolean. If TRUE, a result file in GFF format (see
  binSize parameter) will be generated. Typically useless except for
  compatibility with some peak detection algorithms such as ones used for
  chIP-on-chip experiments}
  
  \item{incrArtefactThrEvery}{A complex parameter (see details). A numeric
  value or NA. A positive numeric value defines a threshold to consider piles
  like 'artefacts' as 'number of reads in the experiment devided by the
  value'. A negative value allows to directly set the threshold manually,
  independently of the number of reads. A NA will ignore the eventual
  artefactual piles}
  
  \item{binSize}{A complex parameter (see details). A striclty positive value
  (>0) defining the size of the steps in the resulting WIG file. When 1 there
  is one score per genomic coordinate (no binning)}
  
  \item{elongationSize}{A complex parameter (see details). A numeric value
  that tells how much each read must be extended to fit the actual insert
  size. If NA, the elongation estimation module will be used to automatically
  determine the overrepresented insert size in the experiment. If 0, the
  reads will not be extended. For paired-end experiments, this value will
  override the real insert size, except if NA value is used}
  
  \item{rangeSelection}{A complex parameter (see details). An 'IRanges'
  object that defines the different groups of reads to be made. In case of
  single-end experiments, these ranges are applied to reads size whereas in
  case of paired-ends experiments, the groups are made based on the insert
  size. An empty range (as defined for default value) will not perform any
  selection, using all the reads independently of their size.}
  
  \item{annotationFilesGFF}{A complex parameter (see details). A named vector
  of gff file paths. If this argument and annotationGenomeFiles are provided,
  the pipeline will generate or each rangeSelection (and total) a pdf file
  summarizing reads occupancy among annotations in gff files.}
  
  \item{annotationGenomeFiles}{A complex parameter (see details). Needed for
  plotting reads statistics on annotations (see argument 'annotationFilesGFF').
  A single path to a genome reference file or the ID of a genome for which
  a reference is bundled in the package (to see bundled files use the command :
  'dir(system.file("resources", package="Pasha"))'.
  IMPORTANT : in order to use a bundled reference file, one must NOT specify
  the file extension (examples : hg18, hg19, mm9...)}
  
  \item{elongationEstimationRange}{A numeric vector with 3 named elements :
  'mini', 'maxi', 'by'. These values define the range and the granularity
  that will be used for elongation estimation}
  
  \item{rehabilitationStep}{A character vector. Can contain 'orphans',
  'orphansFromArtefacts'. See 'Rehabilitation steps' in 'details' section}
  
  \item{removeChrNamesContaining}{A single regular expression. Defines a
  pattern that will match chromosome (seqnames) names to be removed from the
  experiment}
  
  \item{ignoreInsertsOver}{A single strictly positive integer, or NA. In
  case of paired-ends experiments, one might want to ignore inserts above a
  certain size (which are probably the result of misalignment)}
  
  \item{nbCPUs}{A single strictly positive integer.If several cores are
  available, the program can work on several chromosomes in parallel. This
  decreases the time needed for processing but uses more memory}
  
  \item{keepTemp}{A single logical. If TRUE, after the merge, the
  intermediary files will not be erased. It concerns the wig files of
  subgroups such as 'orphans', 'orphansFromArtefacts', 'multireads'. If FALSE
  only the merged result will be kept}
  
  \item{logTofile}{A single string, or NULL. If the string defines a valid
  filename, the log messages for all experiments will be written in it. Note
  that there is a local copy of the log for each experiment in the result
  folder}
  
  \item{eraseLog}{A single logical. In its default behavior, the pipeline
  stops if a log file with the same name already exists. This is a security
  to avoid deleting previously computed results. One can disable this
  security by putting this parameter as TRUE}
  
  \item{multiLocFilesList}{A named (or empty) list. Names of elements must
  match the experiments names defined in 'INPUTFilesList'. Each element must
  be the full path to a text file resulting from multiread processing (see
  multiread specific commands in the package)}
} 
    
\details{ 
The 'pipeline' covers most of the required steps to convert chromatin aligned
reads information to piled-up enrichments scores (ie. WIG files).

It makes use of the package functions to :
\itemize{
    \item{Read and store in memory dataset(s) provided by the user}
    \item{Generate eventual subgroups of interest based on the insert size (in
    case of paired-end experiment) or on the reads size, providing statistics
    about the regions covered and the number of reads concerned}
    \item{Help to identify and remove artefactual enrichments}
    \item{Estimate in-silico the size of original DNA fragments (inserts, in
    case of single-end experiments), or plotting the insert size distribution
    otherwise}
    \item{Extension and piling of reads with several options to obtain a score
  per genomic coordinates}
    \item{Writing results as WIG variable, WIG fixed-step, or GFF files}
}


The flexible implementation allows users to specify a bunch of parameters values  
for each experiment (using an advanced options system, see details). The
pipeline will then loop over the different values and provide results in
separate folders.

It also provides mechanisms to handle orphan reads (incomplete pairs) all along
the process and deal with reads that aligned in multiple locations (see specific
functions below).

The work can also be spread automatically on several processors/cores (at the
expense of memory) by taking advantage of the parallelization options offered by
'parallel'.



The pipeline (started by processPipeline) defines a workflow using the main
functions provided in the package with an additional layer dedicated to
automatization, graphics/QC, summary of experiments, and organization of
resulting files.



How to specify experiment files parameters 'INPUTFilesList'

It is possible to define a list of experiment that will be processed
sequentially.
This argument is a list of list. Each element of this list is freely named by
the user in order to recognize the experiment. These names will be used to name
the experiment in logs and results.
Each element contains a nested list with fixed names that must match the
 followings : 'folderName', 'fileName', 'fileType', 'chrPrefix', 'chrSuffix',
'pairedEnds', 'midPoint' 

 \describe{
    \item{folderName}{is the full path to the folder containing the file
    containing aligned reads to be loaded.}
    
    \item{fileName}{is the filename containing the aligned data, typically the
    output from a program of read alignment.}
    
    \item{fileType}{defines the format in which the data is stored in the file.
    The pipelinehandle preferentially standard BAM files ('BAM' value), but can
    also read several proprietary format by using the ShortRead library (see
    readAligned function from ShortRead library for a complete list of
    supported formats).}
    
    \item{chrPrefix}{is a regular expression string defining the prefix
    concatenated to each chromosome name. Typically this value is 'chr' but can
    be different depending on how the reference genjome for alignment has been
    created and named.}
    
    \item{chrSuffix}{see 'chrPrefix' argument. Typically, this is defined as an
    empty string ''.}
    
    \item{pairedEnds}{a logical value specifying is the experiment should be
    considered as single or paired -end. In the latter case, only BAM files are
    supported and the format has to be cin conformity with BAM policy, such as
    bowtie output.}
    
    \item{midPoint}{a logical value. If TRUE a specific piling method will be
    applied to the reads and the output files will have '_midpoint' suffix. See
    below for details on midpoint method.}
    }

Graphics

Depending on the type of experiment, several graphics can be produced in the
output log folder.
In case of paired-end experiments, a plot will be produced to give information
about the distribution of inserts size.
For single-end experiments, if no manual elongation size is provided
(recommended), the elongation estimation module will produce a plot summarizing
the estimation for each chromosome.
If several subgroups (based on insert size or reads size) are asked, a global
distribution plot will be produced on which the concerned subpopulation is
highlighted.
The 'artefacts' module will also produce a summary sheet for each chromosome
with several graphics (for each strand separately). First, a distribution of
the piles (reads sharing the exact same coordinates) heights on the chromosome.
Second, a cumulative plot describing the proportion of reads in the total
population that are lost for several choices of threshold (see
getArtefactIndexes function). Finally, a summary of the number of reads removed
for different thresholds and the one choosen by the user (plus more information
about orphan reads in case of paired-ends experiments).
Finally if arguments 'annotationFilesGFF' and  'annotationGenomeFiles' are
specified and valid, a report will be generated to summarize the distribution
of the reads in the annotation files (for each range in case rangeSelection
argument is provided).


Experiments summary

During the processing, the pipeline provide a lot of information about the
experiment and the current operations.
All these informations are saved and stored in the log files. Among these
information, one will find the number of reads in the experiment (aligned and
unaligned if provided in the original file), a summary of the reads size,
alignment scores, and all parameters that were asked.


Results File Organization

Because the pipeline allow to specify a lot of different parameters, a
hierarchical organization of resulting files is produced.
First, the 'resultSubFolder' will be produced in the folder where the input
file is after appending a suffix (_PE for paired-ends experiments and _SE for
single-end ones).
In the latter one, another subfolder is created for each subgroup (based on
insert size or reads size). all final results (merged WIGs and GFFs) will be
generated in these folders. A subfolder named as specified in
'reportFilesSubFolder' will also be created and filled with figures and log
files. Another subfolder will be created to store temporary output files (see
keepTemp parameter).


Complex parameters

The pipeline automatization lies on a specific mechanism for some parameters
('incrArtefactThrEvery', 'binSize', 'elongationSize' and 'rangeSelection').
First, these arguments can be used as single values. In this case, the same
value will be used for all experiments described in 'INPUTFilesList' argument.
Second, one can specify a different value for each experiment defined in
'INPUTFilesList' argument. For this concern, the user has to create a named
list with as many elements as number of experiments, and the same names as
'INPUTFilesList' argument.
Finally, it is possible to give several values for each or all experiments. In
this case, the pipeline will loop over these values, recompute what has to be
recomputed, and write resulting files with differentiated filenames.
A specific function is dedicated for checking the consistency of parameters
provided for one or all experiments, and raise an error in case of conflictual
or invalid values.


Rehabilitation steps

For paired-ends experiments, the pipeline tries to save the 'half-pairs'
(called 'orphans'). In brief, reads from incomplete pairs (mate missing because
of misalignment) are separated from the other ones.
Later, when the average insert size has been estimated, it is used for
rehabilitation of incomplete pairs by a step of elongation prior to piling of
these reads. These pileups subcategories are then merged with the other reads
to produce the final result files.
Because pairs can also be broken when eliminating PCR artefacts, a second
category called 'orphansFromArtefacts' are treated the same way, producing
another temporary output that is merged with the other reads.
The user is asked to consider or ignore these reads with the argument
'rehabilitationStep'. Finally the user can also look at the separate categories
result files (WIGs and GFFs) prior merging by using the keepTemp parameter.


Midpoint

the midpoint piling strategy has mainly been thought for MNase experiments
which allows to track for nucleosome positionning. These experiments assume
that each sequenced DNA fragment originate from DNA wrapped around nucleosome
and that non-protected DNA is digested. Thus extremities of fragments are
supposed to represent the boudaries of nucleosomes.
Instead of a classic elongation and piling (which represents the nucleosome
density), the midpoint strategy will only consider the extremities or the exact
center of each DNA fragment (depending on 'elongation parameter', see
package-attached pdf document for a summary), giving a more accurate view of
nucleosome positionning as opposed to the nucleosome density typically
observed.


Range selection

When using the argument 'rangeSelection' with either a numeric value or an
IRanges object, the pipeline will consider separately the reads (or pair of
reads) with different size.
If a single number 'n' is specified, the program will attempt to cut the size
distribution in 'n' intervals with an equal number of reads/pairs in each
(quantile function is used), and process each group separately. When there is
more than one group or if the group does not include all the reads/pairs, the
pipeline will automatically add another group containing all reads (referred as
'allReads').
When the user defines IRanges objects for 'rangeSelection' argument, the
pipeline will process separately each group of size defined by the intervals in
the object. NOTE : the lower value of each interval is EXCLUDED from the group
selection, whereas the upper one is INCLUDED. In this case, only the groups
defines by user are processed (no 'allReads' is added).
Results from each group will be written in separate folders identified by the
concerned range of size.


------------------------------------------

Multireads

The functions dedicated to multiread produce output that can then be injected
in the main pipeline. This information will be integrated in the signal to
produce the final results (merged WIG/GFF files) taking in account the
multireads.


Step 1 - Align the reads
------------------------
To align the reads, you will run Bowtie command with the options requested to
get read aligned in several location reported.
In order to get an output file not too big (since multiread alignement can
bring very large files), the "--concise" option of Bowtie must be used. Other
options are standard.
For instance, if you want to keep the read aligning on the mm9 genomes 100 times
or less, authorizing 2 mismatches and run on 8 processors

  bowtie -q -v 2 -a -m 100 -p 8 --concise mm9 input_file.fastq > output_file.bow

Step 2 - Prepare the genome reference file
------------------------------------------
The scripts used in the next steps requires a reference file containing
information on the chosen genome.
This reference file is a simple text file with ".ref" extension that must
contain:

\itemize{
  \item{first line the number of chromosomes in the considered genome}
  \item{second line, the list of chromosomes names (separated by spaces)}
  \item{third line, the list of chromosomes sizes (separated by spaces)}
}

You can use the 'bowtie-inspect' command on the reference genome to get the
information requested to build this reference file.

Step 3 - Dispatching signal along multiread positions
-----------------------------------------------------
You have now to decide how each multiread read score will be dispatched along
its mutiple positions. Two options are offered here:

* Uniform method: The signal is equally dispatched along the read position i.e.
for instance, if a read has 10 aligned positions, each position will receive a
score of 1/10. To apply this method, execute the R command:

multiread_UniformDispatch(alignedFile, outputFolder, referenceFile,
incrArtefactThrEvery, verbosity)

where :

\describe{
    \item{alignedFile}{the .bow file outputed by bowtie in step 3 (or in step 1
    if you decided not to remove artifacts)}   
    \item{outputFolder}{the folder where you want to see the output file fall
    in.}
    \item{referenceFile}{the .ref file created in step 2}
    \item{incrArtefactThrEvery}{the ratio used to detect artifact. Default value
    is NA.}
    \item{verbosity}{the verbose level (0 = no message, 1 = trace level).
    Default value is 0.}
    }


This command will output a file named with the alignedFile name suffixed by
"_uniformDispatch.txt".

* Chung et al. method: The score on read position is allocated according to an
algorithm developped by Chung et al. (see "Discovering Transcription Factor
Binding Sites of Genomes with Multi-Read Analysis of ChIP-seq Data" (2011) PLoS
Computational Biology).  To apply this method, execute the R command:

multiread_CSEMDispatch( alignedFile, outputFolder, referenceFile, window_size,
iteration_number, incrArtefactThrEvery, verbosity)

where :

\describe{
    \item{alignedFile}{the .bow file outputed by bowtie in step 3 (or in step 1
    if you decided not to remove artifacts)}
    \item{outputFolder}{the folder where you want to see the output file fall
    in.}
    \item{referenceFile}{the .ref file created in step 2}
    \item{window_size}{the size of the window used by the algorithm (see
    algorithm details). Suggested value is 101.}
    \item{iteration_number}{the number of iteration executed by the algorithm.
    Suggested value is 200.}
    \item{incrArtefactThrEvery}{the ratio used to detect artifact. Default value
    is NA.}
    \item{verbosity}{the verbose level (0 = no message, 1 = trace level).
    Default value is 0.}
    }

This command will output a file named with the alignedFile name suffixed by
"_csemDispatch.txt".

Important note:
For both commands, if the parameter 'incrArtefactThrEvery' is set to a strictly
positive value, the input file is first passed to a script removing the
'artifacts'. This script looks at the reads aligned at the exact same position.
If the number of such reads exceed a limit, these reads are considered as due to
experiment artifact.
In that case, only one read is kept and the rest is removed.
The limit defining artefact is controled by the parameter
'incrArtefactThrEvery' such as:

limit = Total number of reads / incrArtefactThrEvery

Step 4 - Launch Pasha for multiread
----------------------------------
One you have obtained a file containing the information of the dispatched
scores over the read positions, you can use it in Pasha by passing the obtained
file in the variable 'multiLocFilesList'. See details above on
'multiLocFilesList' parameter.

Step 5 - Analyze repeat distribution
------------------------------------
Once Pasha has been launched, you have obtained WIG files. Running the
following command will permit you to obtain information and statistics about
the dispertion of your signal on annotations identified as repeated regions in
the genome.

Use the following command to launch the script:

WigRepeatAnalyzer(filename, inputFolder, outputFolder, repeatMaskerFilePath, 
isRegex)

where:


\describe{
    \item{filename}{the file name of the wig file (fixed step WIG).}
    \item{inputFolder}{the path to the wig file.}
    \item{outputFolder}{the path to the folder where analysis results must be
    stored.}
    \item{repeatMaskerFilePath}{the path to the file containing the repeat
    annotations (Repeat Masker file).}
    \item{isRegex}{If TRUE, the filename parameter is interpreted as a regular
    expression (LC_SYNTAX) and the script will search for a unique file
    corresponding to the provided regular expression. If no or several file are
    found, the scripts ends with error.} 
    }

The script compute the coverage of each repeat class and family (i.e. the
percentage of positions falling into each annotations) and the weight of each
class and family (i.e. the percentage of score falling into each annotations).
All results are provided as barplots figures and text files.
}

\value{
A list containing information referring to all processed experiments. Each
element contains at least one nested list element called "execTime", which
summarize the time spent on important data processing steps (see log file).
}
\author{Romain Fenouil}
\seealso{
  \code{\link[Pasha]{getArtefactsIndexes}}
  \code{\link[Pasha]{estimateElongationSize}}
  \code{\link[Pasha]{generatePiled}} 
}
\examples{
  \dontrun{

  # This first (don't run) part of example aims at presenting all
  # available parameters and some variations, see at the bottom of
  # this section for a running example
  
  # Define an experiment description list with a classic epigenetic
  # mark and one MNase experiment that will be seen as nucleosome
  # 'density' or nucleosome 'positionning' (see midpoint parameter)
  
  myExps <- list()
  myExps[["mES_H3K4me3"]] <- list('folderName'="/home/exp",
                                  'fileName'="SRR432543.BAM",
                                  'fileType'="BAM",
                                  'chrPrefix'="chr",
                                  'chrSuffix'="",
                                  'pairedEnds'=FALSE,
                                  'midPoint'=FALSE)
                                             
  myExps[["mES_MNase"]] <- list('folderName'="/home/exp",
                                'fileName'="SPT543426.BAM",
                                'fileType'="BAM", 
                                'chrPrefix'="chr", 
                                'chrSuffix'="", 
                                'pairedEnds'=TRUE, 
                                'midPoint'=FALSE)
                                              
  myExps[["mES_MNase_MIDPOINT"]] <- list('folderName'="/home/exp",
                                         'fileName'="SPT543426.BAM",
                                         'fileType'="BAM", 
                                         'chrPrefix'="chr", 
                                         'chrSuffix'="",
                                         'pairedEnds'=TRUE, 
                                         'midPoint'=TRUE)
  
  # Call the pipeline for the three experiments with basic parameters
      
  processPipeline(
  
  #### I/O GENERAL PARAMETERS
  
  # Experiments description list
  INPUTFilesList              = myExps,
  
  # name of the folder that will contain results
  resultSubFolder             = "Results",
          
  # name of the folder that wil contain logs and figures
  reportFilesSubFolder        = "ReportFiles",
      
  # generate results as WIG fixed steps
  WIGfs                       = TRUE,
               
  # generate results as WIG variable steps
  WIGvs                       = TRUE,
               
  # generate results as GFF files
  GFF                         = FALSE,
              
  #### COMPLEX PARAMETERS (SINGLE OR VECTORS OR LIST OF IT)
  
  # The threshold to detect artefactual piles will be incremented
  # by one every 10Million reads aligned for each experiment
  incrArtefactThrEvery        = 10000000,
           
  # Along the genome one score every 50 basepairs will be computed
  binSize                     = 50,
                 
  # The reads will be extended according to the in-silico estimation
  # algorithm or based on the pairs alignments (insert size) 
  elongationSize              = NA,
   
  # No subgroups selection for specific inserts or reads size
  rangeSelection              = IRanges(0,-1),
  
  # no GFF files given, the module plotting statistics on reads and
  # annotations will not be loaded
  annotationFilesGFF          = NA,
      
  # path to file or "mm9", "hg19"... This argument is needed only if 
  # gff files are specified in 'annotationFilesGFF' argument
  annotationGenomeFiles       = NA,
  
  #### SINGLE PARAMETERS
  
  # For single-end experiments, the fragment size will be estimated
  # between 50 and 400 with a resolution of 10bp
  elongationEstimationRange   = c(mini=50, maxi=400, by=10),
   
  # The pipeline will try to save half-pairs from alignment and the 
  # ones broken during 'artefact' removal
  rehabilitationStep          = c("orphans","orphansFromArtefacts"),
  
  # Remove chromosomes with names containing "random" or "un" 
  removeChrNamesContaining    = "random|un",
        
  # For paired-ends ignore inserts > 500bp according to alignment
  ignoreInsertsOver           = 500,
  
  # Use 1 cpu (recommended as first try to estimate the memory usage)              
  nbCPUs                      = 1,
                  
  # Do not erase pileups and result files from subcategories prior to
  # merging (orphans etc...)
  keepTemp                    = TRUE,
               
  # make a copy of the log for all experiments
  logTofile                   = "./log.txt",
        
  # In case the same computation is restarted, do not warn the user
  # and erase previous results
  eraseLog                    = TRUE,
               
  #### LIST PARAMETERS (one element per expName)
  
  # An eventual list of multiread repartition results
  multiLocFilesList           = "");
                
  ########
  # The four "complex parameters" could have been declared like this
  # to generate more results
  
  # as a vector of values, each value will be used sequentially for
  # all experiments
  # incrArtefactThrEvery <- c(10000000,NA, -10)
  
  # as a list for specifying arguments for individual experiments
  # binSize              <- list("mES_H3K4me3"=200, 
  #                                   "mES_MNase"=50,
                                  "mES_MNase_MIDPOINT"=50)
  
  # mixed, some experiment have one value, others have several
  # elongationSize       <- list("mES_H3K4me3"=c(NA,0),
  #                                   "mES_MNase"=c(146,NA),
  #                                   "mES_MNase_MIDPOINT"=NA) 
                                     
  # Compute without elongating reads (0), a fixed numeric value (not
  # recommended), or estimate in-silico (or based on pairs) the
  # optimal elongation (NA) 
  # rangeSelection <- list("mES_H3K4me3" =IRanges(0,-1), 
  #                        "mES_MNase"=c(IRanges(0,-1),
  #                                      IRanges(0,100), 
  #                                      IRanges(100,1000)),
  #                        "mES_MNase_MIDPOINT"=c(IRanges(0,-1),
  #                                               IRanges(0,100), 
  #                                               IRanges(100,1000)))
  
  
}  


#############################################
#### Actual runnable example on BAM file ####
#############################################

# Define a temporary directory where the example will run
exampleFolder <- tempdir()

# Get the path to the example BAM file and copy it (with the index)
testFileBAM_fileName <- "embedDataTest.bam"
testFileBAM_fullPath <- system.file("extdata",
                                    testFileBAM_fileName,
                                    package="Pasha")
file.copy(testFileBAM_fullPath, exampleFolder)

testFileBAI_fileName <- "embedDataTest.bam.bai"
testFileBAI_fullPath <- system.file("extdata",
                                    testFileBAI_fileName,
                                    package="Pasha")
file.copy(testFileBAI_fullPath, exampleFolder)
# Create the data structure containing information on the experiments
INPUTFilesList <- list()
INPUTFilesList[["testBAM"]] <- list(folderName=exampleFolder, 
                                    fileName=testFileBAM_fileName,
                                    fileType="BAM", 
                                    chrPrefix="chr", 
                                    chrSuffix="", 
                                    pairedEnds=TRUE, 
                                    midPoint=FALSE)
\dontrun{
# Start the pipeline using default parameters
processPipeline(INPUTFilesList)
}
}