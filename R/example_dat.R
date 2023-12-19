#' A List of Simulated Longitudinal Data
#'
#' The \code{example_dat} is a list of two simulated longitudinal data.
#'
#' @details The first object is the \code{biomk} data which mimics
#' the scaled cerebrospinal fluid biomarker level distribution,
#' containing the baseline age, gender, biomarker levels, and visit times
#' for 100 participants.
#' The second list object is the \code{primary} data which mimics
#' the distribution of the cognitive scores,
#' containing the baseline age, baseline education, gender, cognitive scores,
#' time for each available score, and a categorical variable indicating if an
#' individual has more than a copy of a certain gene for the same 100
#' participants.
#'
#'
#' @format The data frames contain the following columns:
#'
#' \describe{
#' \item{biomk: id}{participant identifer}
#' \item{biomk: time}{date of the biomarker measurement}
#' \item{biomk: gender}{an indicator of gender (1 = male, 0 = female)}
#' \item{biomk: baseline_age}{baseline age}
#' \item{biomk: biomk}{biomarker level}
#'
#' \item{primary: id}{participant identifer}
#' \item{primary: time}{date of cognitive test}
#' \item{primary: gender}{an indicator of gender (1 = male, 0 = female)}
#' \item{primary: baseline_age}{baseline age}
#' \item{primary: gene}{an indicator of gene copy number (1 = at least one copy,
#'  0 = 0 copy)}
#' \item{primary: baseline_edu}{number of educational years completed at baseline}
#' \item{primary: score}{cognitive score}
#' }
#'
#' @keywords datasets
#' @examples
#' \dontrun{
#' head(example_dat$biomk)
#' head(example_dat$primary)
#' }
#'
"example_dat"
