\name{binVector}
\alias{binVector}
\alias{C_binVector}
\title{Divide a numeric vector in equal bins of signal}
\description{
Reduce the size of a numeric vector (typically large) to a more compact one by
averaging the numeric values of 'binSize' (constant) consecutive cells.
}
\usage{
binVector(piledValues, binSize)
}
\arguments{
  \item{piledValues}{A numeric vector}
  \item{binSize}{Size of the bins}
}
\details{
This function is typically used to represent piling scores for large
chromosomes, specially if the scores have to be written in a file (WIG fixed
step).
When the last bin has fewer number of items, an average is also computed
but only on the items available.
}
\value{
A numeric vector filled by the C function
}
\author{Romain Fenouil}
\examples{
# Create a vector of values and plot it
myValues <- sin(seq(1,10, length.out=300))*10
plot(x=1:length(myValues), myValues, type="b", pch=".", cex=2)

# Bin the vector by chunks of 20 points (compute the average value)
binnedValues <- binVector(myValues,20)

# Plot the binned value on top (with centered coordinates)
points(((x=1:length(binnedValues))*20)-10, 
         y=binnedValues, 
         col="blue",
         type="b",
         pch=5)
         
}
