\name{run_spatial_app}
\alias{run_spatial_app}
\title{Launch the SpatialOmicsMSI Shiny Application}
\description{
Launches the Shiny application bundled with the installed package.
}
\usage{
run_spatial_app(launch.browser = TRUE, port = NULL, host = "127.0.0.1")
}
\arguments{
  \item{launch.browser}{Whether to open a browser when the application starts.}
  \item{port}{Optional integer TCP port from 1 to 65535. When \code{NULL},
  Shiny selects an available port.}
  \item{host}{Address on which to listen. Use \code{"127.0.0.1"} with an SSH
  tunnel, or \code{"0.0.0.0"} to accept connections on all server interfaces.}
}
\value{
The value returned by \code{shiny::runApp()}; the function normally blocks
until the application stops.
}
