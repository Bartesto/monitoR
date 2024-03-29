% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_prep.r
\name{just_nums}
\alias{just_nums}
\title{Creates correct suffix for day of month}
\usage{
just_nums(n)
}
\arguments{
\item{n}{numeric value for day of month.}
}
\value{
A character value with correct suffix.
}
\description{
\code{just_nums} takes a numeric value and returns a character equivalent
    with the correct suffix, e.g. 'th', 'nd' etc.
}
\examples{
just_nums(n = 2)

}
\author{
Bart Huntley, \email{bart.huntley@dbca.wa.gov.au}
}
