% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/display-stars.r
\name{display_stars}
\alias{display_stars}
\alias{animate_stars}
\title{Star glyph tour path animation.}
\usage{
display_stars(...)

animate_stars(data, tour_path = grand_tour(3), ...)
}
\arguments{
\item{...}{other arguments passed on to \code{\link{stars}}}

\item{data}{matrix, or data frame containing numeric columns}

\item{tour_path}{tour path generator, defaults to 2d grand tour}
}
\description{
Animate a nD tour path with star glyphs.
}
\details{
Currently, scaling doesn't seem to be computed absolutely correctly, as
centres move around as well as outside points.
}
\examples{
animate_stars(flea[1:10, 1:6])
animate_stars(flea[1:10, 1:6], grand_tour(5))
animate_stars(flea[, 1:6], grand_tour(5))
animate_stars(flea[1:10, 1:6], grand_tour(5),
  col.stars = rep("grey50", 10), radius = FALSE
)
}
\seealso{
\code{\link{animate}} for options that apply to all animations
}
\keyword{hplot}
