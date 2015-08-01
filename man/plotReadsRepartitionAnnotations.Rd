\name{plotReadsRepartitionAnnotations}
\alias{plotReadsRepartitionAnnotations}
\title{Graphical information about reads and annotations repartition}

\description{ This function plots from an internal alignedData object the number
of annotations per chromosome, the number of annotations by file, the number of
reads per chromosome, the length of chromosomes, the percentage of coverage of
each annotation according to the size of the reference genome and the number of
reads of the reference experiment that aligned to a given annotation.
} 
\usage{ 
plotReadsRepartitionAnnotations(alignedData,  
                                gff_names_vec, 
                                expName,
                                pdfFileName, 
                                genomeReferenceFile) 
      } 

\arguments{
  \item{alignedData}{An internal R object containing information about the
  aligned reads}
  \item{gff_names_vec}{A vector containing all file pathes to annotations in
  gff format}
  \item{expName}{A character string giving the name of the reference experiment
  from which information about reads alignment were retrieved}
  \item{pdfFileName}{The file path to the pdf output}
  \item{genomeReferenceFile}{File path to the reference genome in .ref format
  (see details)} } 

\details{ The pdf output is divided into three pages:

- The first page contains 4 barplots giving respectively the number of each
annotation by chromosomes (seqnames) (the correspondence between colors and
annotations is given at the bottom), the number of annotations contained in each file of
gff_names_vec argument, the number of reads that aligned to each chromosome and
the length of each chromosome.

- The second page contains two pie charts: The first one gives the percentage of
coverage of the reference genome by each annotation, correspondence of colors is
given at the right of the pie chart. The second pie chart gives the percentage
of reads that aligned to a particular annotation, correspondence of colors is
also given at the right of the pie chart.

- The third page is another representation of the pie chart enabling to compare
directly the difference of occupancy of each annotation with the number of reads
aligned to it.

The 'genomeReferenceFile' is the path to a file containing information about
the genome that has been used to align the reads. It is a 3 lines text files.

The same files are used for multiread (see vignette).

Some files are provided in the package for some of the most commonly used
genomes. For a list of them type: 

dir(pattern="*.ref", path=system.file("resources", package="Pasha")) 

First line describes the number of chromosomes in the genomes. 
Second line contains the size (in bp) of all chromosomes (space separated). 
Last lines containes the names of the chromosomes (space separated).

Example :
3
43255 10345 2456
chrI chrII chrIII
}

\value{
A pdf file with different graphics representing information described in the details section.
}
\author{Nicolas Descostes}
\keyword{ reads annotations overlapping repartition }
