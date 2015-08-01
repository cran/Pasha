\name{AlignedData-class}
\docType{class}

% Class
\alias{AlignedData-class}
\alias{AlignedData}

% Constructors
\alias{readAlignedData}
\alias{readMultipleAlignedData}

% Methods
\alias{mapq<-,AlignedData-method}
\alias{mapq,AlignedData-method}
\alias{[,AlignedData,ANY,ANY-method}
\alias{checkPairsOK,AlignedData-method}
\alias{seqnames<-,AlignedData-method}
\alias{seqnames,AlignedData-method}
\alias{dropChromosomePattern,AlignedData-method}
\alias{cleanNonSimplePairs,AlignedData-method}
\alias{filterInsertSize,AlignedData-method}
\alias{filterReadSize,AlignedData-method}
\alias{flag<-,AlignedData-method}
\alias{flag,AlignedData-method}
\alias{getOrphansIndexes,AlignedData-method}
\alias{getStats,AlignedData-method}
\alias{readID<-,AlignedData-method}
\alias{readID,AlignedData-method}
\alias{isize<-,AlignedData-method}
\alias{isize,AlignedData-method}
\alias{length,AlignedData-method}
\alias{mrnm<-,AlignedData-method}
\alias{mrnm,AlignedData-method}
\alias{normalizeChrNames,AlignedData-method}
\alias{notAligned<-,AlignedData-method}
\alias{notAligned,AlignedData-method}
\alias{pairedEnds<-,AlignedData-method}
\alias{pairedEnds,AlignedData-method}
\alias{position<-,AlignedData-method}
\alias{position,AlignedData-method}
\alias{qwidth<-,AlignedData-method}
\alias{qwidth,AlignedData-method}
\alias{dropUndefinedStrand,AlignedData-method}
\alias{sortByPairs,AlignedData-method}
\alias{strand<-,AlignedData-method}
\alias{strand,AlignedData-method}
\alias{weight<-,AlignedData-method}
\alias{weight,AlignedData-method}

% Functions
\alias{mapq}
\alias{mapq<-}
\alias{checkPairsOK}
\alias{seqnames}
\alias{seqnames<-}
\alias{dropChromosomePattern}
\alias{cleanNonSimplePairs}
\alias{filterInsertSize}
\alias{filterReadSize}
\alias{flag}
\alias{flag<-}
\alias{getOrphansIndexes}
\alias{readID}
\alias{readID<-}
\alias{initialize}
\alias{isize}
\alias{isize<-}
\alias{length}
\alias{mrnm}
\alias{mrnm<-}
\alias{normalizeChrNames}
\alias{notAligned}
\alias{notAligned<-}
\alias{pairedEnds}
\alias{pairedEnds<-}
\alias{position}
\alias{position<-}
\alias{qwidth}
\alias{qwidth<-}
\alias{readAlignedData}
\alias{dropUndefinedStrand}
\alias{sortByPairs}
\alias{strand}
\alias{strand<-}
\alias{weight}
\alias{weight<-}


\title{Class \code{AlignedData}}

\description{ This class gives tools to read, store and manage simple data from
HTS experiments. This class is specially optimized for chromatin-oriented
sequencing experiments (ChIP-Seq, MNase-Seq, Faire-Seq) since it is mainly
focused on reads position rather than sequences, ignores alignments with
non-simple CIGAR strings (non canonical alignments, ie. insertions/deletions)
and secondary alignments, and is not designed to represent features with
undefined strands.
It directly uses Rsamtools to read data from BAM files. It has a similar global
structure to ShortRead and GAlignments classes but have some specific
features.
}

\details{The specific features of this class are that (i) read sequences are not
read from aligned files, (ii) only R native types are used to store data in
memory, and (iii) it can represent both paired-end or single-end data by using a
similar internal structure as what has been defined for SAM/BAM formats. For
reading data from 'non-BAM' formats (aligner-specific), the constructor relies
on ShortRead library (see ShortRead for details).

This class is a simplistic representation of aligned reads for users aiming at
using piling options offered by Pasha package. Secondary alignments, reads
aligned with insertions/deletions (non simple CIGAR), and not aligned reads are
discarded by the constructor, and reads with a strand value other than -/+ will
raise a warning (strand representation keeps levels c('+', '-', '*') for
consistency with bioconductor packages). For other applications requiring
exhaustive handling of CIGAR strings (insertions/deletions), undefined strand
values, or to keep reads that could not be aligned, the use of GAlignments class
family is recommended.

Specific functions are provided for a seamless conversion to/from GAlignments
class. It is however important to note that eventual reads with non-simple CIGAR
string will be dropped (with a warning) during conversions from GAlignments
objects. 

Internally, this class defines different kinds of slots. First, slots where a
value is mandatory for each read. Second, slots with atomic values (one
for the whole object : "notAligned", "pairedEnds"), and finally optionnal slots
which can be empty (no information) or should provide information for each read
("mapq", "isize", "flag", "mrnm", "weight", "qwidth").


\describe{
  \item{readID:}{integer - Read identifier. This value is inferred from
  'groupid' field when reading paired-end BAM files (see
  \link[Rsamtools]{BamFile} and 'asMates'). A simple integer if inferred from
  'qname' field otherwise.)}
  \item{seqnames:}{factor - Chromosome sequence name}
  \item{strand:}{factor - Strand value (accessors and subsetters methods will
  try to maintain consistent values for factor levels ('+','-'))}
  \item{position:}{integer - The 1-based leftmost position/coordinate (POS field
  in a SAM/BAM record)}
  \item{mapq:}{integer - Mapping quality}
  \item{isize:}{integer - Insert size (paired-ends only)}
  \item{qwidth:}{integer - Query width (read sequence length)}
  \item{flag:}{integer - Flag (see SAM format specification)}
  \item{mrnm:}{factor - Mate-pair chromosome sequence name}
  \item{notAligned:}{integer (single value) - Number of reads that wouldn't be
  aligned during mapping (if information availabe)}
  \item{pairedEnds:}{logical (single value) - Object declared as paired-ends}
  \item{weight:}{numeric - Weight attributed to read for piling (multireads
  only)} }


Paired-end specificities

The limitations imposed by the constructor when reading paired-end datasets (no
unaligned reads, and no secondary alignments) allow to represent the mate reads
association by their shared 'qname' (reported as 'readID'). Altering such IDs
can thus break the pairs relation as seen by methods of the class.

Several methods are dedicated for paired-ends datasets (see 'checkPairsOK',
'filterInsertSize', 'getOrphansIndexes', 'cleanNonSimplePairs', 'sortByPairs').
These methods will emit a warning and not process data when called for a
single-end object.

}


\section{Objects from the Class}{Objects can be created by a call of the form
\code{new("AlignedData", readID, seqnames, strand, position, mapq,
isize, qwidth, flag, mrnm, notAligned, pairedEnds, weight)}, or by using the
dedicated constructor that builds an object from a file.}

\section{Constructors}{
  \describe{
    \item{\code{readAlignedData(folderName, fileName, fileType="BAM",
    pairedEnds=FALSE)}}{This function reads data from a file and creates an
    AlignedData object from it. It handles BAM files using the Rsamtools
    library (ignores reads with non-simple CIGAR string) and uses ShortRead
    library to handle several other proprietary formats (see ShortRead
    readAligned function).}
    } 
}

\section{Subsetting}{
  \describe{
    \item{[}{Allow to select a subpopulation of reads, all slots are affected
    except the atomic ones for which the value is copied in resulting object.
    This function can also drop the unused levels for chromosome (seqnames) and
    mrnm slots only if the argument drop=TRUE is used.} 
    } 
}

\section{Coercion}{
  \describe{
  \item{as(object, 'GAlignments')}{Converts an AlignedData object to a
  GAlignment object. This conversion drops the eventual pairs related
  information, as well as mapping quality, number of unaligned read in the
  dataset, and flags.}
  \item{as(object, 'AlignedData')}{Converts a GAlignments object to a
  AlignedData object. This conversion will drop eventual reads that aligned
  with a non-simple CIGAR string.} 
  } 
}
 

\section{Accessors}{
  \describe{
    \item{\code{readID(object)}, \code{readID(object) <- value}}{Get or set the
    numeric ID of reads (typically inferred from read name by constructor, see
    details). Altering these values can break the reads mating for paired-end
    data.}
    
    \item{\code{seqnames(object))}, \code{seqnames(object) <- value}}{Get or
    set a factor describing the seqnames (chromosome name) to which each read
    has been aligned.}
    
    \item{\code{strand(object)}, \code{strand(object) <- value}}{Get
    or set the factor describing the strand to which reads are aligned.}
    
    \item{\code{position(object)}, \code{position(object) <- value}}{Get or set
    the 1-base leftmost genomic coordinate of the reads.}
    
    \item{\code{mapq(object)}, \code{mapq(object) <- value}}{Get
    or set the numeric quality of alignment (if available). See Rsamtools
    library for details.} 
    
    \item{\code{isize(object)}, \code{isize(object) <- value}}{Get or set the
    insert size for paired reads (if available).}
    
    \item{\code{qwidth(object)}, \code{qwidth(object) <- value}}{Get or set the
    length of each read sequence (if available).} 
    
    \item{\code{flag(object)}, \code{flag(object) <- value}}{Get or set the
    numeric (integer) flag values of reads (if available). See BAM/SAM
    documentation for flag interpretation. Modification of this slot can lead to
    inconsistencies in further processing (specially for paired-end
    experiments).}
    
    \item{\code{mrnm(object)}, \code{mrnm(object) <- value}}{Get or set the mate
    read chromosome name (if available).}
    
    \item{\code{notAligned(object)}, \code{notAligned(object) <- value}}{Get or
    set the number of reads reported as not aligned by aligner (only available
    from BAM format when object is created via the constructor).}
    
    \item{\code{pairedEnds(object)}, \code{pairedEnds(object) <- value}}{Get or
    set the declared status of the object concerning paired/single ends reads.}
    
    \item{\code{weight(object)}, \code{weight(object) <- value}}{Get or set the
    weight that has been attributed to each read for piling (if available).} 
    }
}

\section{Methods}{
  \describe{
    \item{\code{show(object)}}{Writes a quick summary of the object content to
    the standard output (uses cat).}  
  
    \item{\code{checkPairsOK(object)}}{Perform a validity test on paired
    objects. Checks that the object only contain 'full-pairs' (no orhpans) and
    that all reads are correctly sorted (eg. interlaced reads and mates).}
    
    \item{\code{length(object)}}{Get the number of elements in the object.}

    \item{\code{dropChromosomePattern(object, pattern, quiet=TRUE)}}{Filter
    out reads for which the chromosome name (seqnames) match a certain pattern.
    Helps to get rid of undesired reads (align on mitochondrial chromosomes for
    instance).}

    \item{\code{filterInsertSize(object, rangeMin, rangeMax, includeLower=FALSE,
    quiet=TRUE)}}{Select pairs of reads for which insert size is defined in a
    specific range. Returns a new object with the selected reads. By default
    (includeLower=FALSE), the lower range value ('rangeMin') is EXCLUDED from
    selection whereas the upper one is INCLUDED in the selection. Insert size
    is defined as a symmetric (positive and negative numeric) value for reads of
    a pair. This function is designed for selecting both reads of pairs and
    therefore ignores the numeric sign for filtering.}

    \item{\code{filterReadSize(object, rangeMin, rangeMax, includeLower=FALSE,
    quiet=TRUE)}}{Select reads for which the sequence length is in a
    specified range. Returns a new object with the selected reads. By default
    (includeLower=FALSE), the lower range value ('rangeMin') is EXCLUDED from
    selection whereas the upper one is INCLUDED in the selection.}

    \item{\code{getOrphansIndexes(object, quiet=TRUE)}}{Get the indexes of
    reads for which a mate can't be identified (based on their ID).}

    \item{\code{cleanNonSimplePairs(object, quiet=TRUE)}}{For paired-ends
    objects, return a filtered object from which the pairs with flag 'unmapped
    mate' and/or the pairs with both reads on the same strand are "orphanized"}

    \item{\code{normalizeChrNames(object, chrPrefix, chrSuffix)}}{The main
    functions of this library (see processPipeline) assumes that chromosomes are
    always represented as chr*, * being the actual chromosome name. This
    functions helps to fit this requirement by removing eventual prefix/suffix
    around the actual chromosome name and then adding 'chr' as prefix. This
    function returns an object with modified seqnames (chromosome names).}

    \item{\code{sortByPairs(object, quiet=TRUE)}}{In case of pairedEnds
    experiments, this function sorts the reads. It tries to group the index of
    reads from pairs and uses the flag to order the reads by 'first in pair' and
    'second in pair' bits as firsts and seconds respectively.}
    
    \item{\code{dropUndefinedStrand, quiet=TRUE}}{Return a AlignedData object
    after removing eventual reads with unsupported strand values (supported strand
    values are '+' and '-').}

    } 
}


\section{Arguments}{
\describe{
    \item{folderName}{Path to the directory containing the file to read.}
    
    \item{fileName}{File containing aligned reads. Several formats are handled.}
    
    \item{fileType}{File format for reads description. This format can vary
    depending on the aligner and/or the parameters used. Typical format will be
    BAM files from which the maximum of information can be read from. Other
    proprietary formats are handled thanks to the ShortRead library (see
    corresponding documentation for details)}
    
    \item{pairedEnds}{A logical defining if the data read from file should be
    considered as paired-end or not.}
    
    \item{object}{Object of the class 'AlignedData' to which the function call
    is applied.}
    
    \item{pattern}{A regular expression describing the pattern to be recognized
    and replaced}
    
    \item{quiet}{An atomic logical. If TRUE, the function will not write
    messages to stdout}
    
    \item{rangeMin}{An atomic integer, defining the minimal size to be selected}
    
    \item{rangeMax}{An atomic integer, defining the maximal size to be selected}
    
    \item{includeLower}{An tomic logical defining if the lower value should be
    included when filtering reads/inserts size based on a range of values
    (default FALSE)}
    
    \item{chrPrefix}{An atomic regular expression describing the pattern to be
    recognized and replaced}
    
    \item{chrSuffix}{An atomic regular expression describing the pattern to be
    recognized and replaced}
    
    \item{breaks}{An integer, or a vector of integer. This parameters behaves
    like the 'breaks' parameter in function 'cut'} 
    }
}


\author{Romain Fenouil}
\seealso{
  \code{\link[ShortRead]{readAligned}}
  \code{\link[GenomicAlignments]{GAlignments}}
  \code{\link[Rsamtools]{scanBam}}
  \code{\link[Rsamtools]{BamFile}}
}
\examples{

# Get the BAM file used in RSamTools as example
BAMfile <- system.file("extdata","ex1.bam",package="Rsamtools")
# Testing the public constructor of AlignedData class (from Pasha)
alignedDataObject <- readAlignedData(folderName="" ,
                                     fileName=BAMfile, 
                                     fileType="BAM", 
                                     pairedEnds=TRUE)
print(alignedDataObject)
                             
# Rename the chromosomes from seqX to chrX
alignedDataObject <- normalizeChrNames(alignedDataObject,
                                       chrPrefix="seq",
                                       chrSuffix="")
print(alignedDataObject)
                               
                               
#### Splitting object (indirect use of '[') 
# Split the object and the associated data by chromosomes 
splitListalignedDataObject <- split(alignedDataObject, 
                                    seqnames(alignedDataObject))
# check that the object has correctly been split
print(names(splitListalignedDataObject))
print(sapply(lapply(splitListalignedDataObject, seqnames), unique))


#### Orphans and pairs consistency
print(length(alignedDataObject))
# Identify the reads with no mate and remove them
orphansReadsIndex <- getOrphansIndexes(alignedDataObject)
print(paste("Nb of orphan reads :", sum(orphansReadsIndex)))
alignedDataObject <- alignedDataObject[!orphansReadsIndex]
print(length(alignedDataObject))

# Check that pairs are correctly represented in the object
# ie. mate pairs are consecutive
print(checkPairsOK(alignedDataObject))
# Sort by pairs
alignedDataObject <- sortByPairs(alignedDataObject)
# Check pairs consistency again
print(checkPairsOK(alignedDataObject))


#### Filtering

## Chromosome patterns
# (drop reads on chromosomes containing the character '1')
print(alignedDataObject)
resultFiltering <- dropChromosomePattern(alignedDataObject, 
                                           pattern="1")
print(resultFiltering)

## Range of insert size
print(table(isize(alignedDataObject)))
resultFiltering <- filterInsertSize(alignedDataObject, 
                                    rangeMin=190, 
                                    rangeMax=195)
print(table(isize(resultFiltering)))
print(resultFiltering)

## Reads sequence size
print(table(qwidth(alignedDataObject)))
resultFiltering <- filterReadSize(alignedDataObject, 
                                  rangeMin=35, 
                                  rangeMax=50)
print(table(qwidth(resultFiltering)))
print(resultFiltering)
                          
## Reads with atypic characteristics for chromatin alignments
# Mate unmapped flag set, or mates aligned on same strand
resultFiltering <- cleanNonSimplePairs(alignedDataObject, quiet=FALSE)
# Reads with a strand value different from expected '+' or '-'
resultFiltering <- dropUndefinedStrand(resultFiltering, quiet=FALSE)

}
\keyword{classes}
