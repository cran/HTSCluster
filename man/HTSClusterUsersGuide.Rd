\name{HTSClusterUsersGuide}
\alias{HTSClusterUsersGuide}

\title{View HTSCluster User's Guide}
\description{Finds the location of the HTSCluster User's Guide and optionally opens it.}

\usage{
HTSClusterUsersGuide(view=TRUE)
}

\arguments{
  \item{view}{logical, should the document be opened using the default PDF document reader?}
}

\value{Character string giving the file location.
If \code{view=TRUE}, the PDF document reader is started and the User's Guide is opened, as a side effect.}

\details{
The function \code{vignette("HTSCluster")} will find the short HTSCluster Vignette which describes how to obtain the HTSCluster User's Guide.
The User's Guide is not itself a true vignette because it is not automatically generated using \code{\link{Sweave}} during the package build process.
This means that it cannot be found using \code{vignette}, hence the need for this special function.

If the operating system is other than Windows, then the PDF viewer used is that given by \code{Sys.getenv("R_PDFVIEWER")}.
The PDF viewer can be changed using \code{Sys.putenv(R_PDFVIEWER=)}.

Note that this function was adapted from that defined by Gordon Smyth in the edgeR package.
}

\seealso{
\code{\link{system}}
}

\author{Gordon Smyth}

\examples{
# To get the location:
HTSClusterUsersGuide(view=FALSE)
# To open in pdf viewer:
\dontrun{HTSClusterUsersGuide()}
}

\keyword{documentation}
