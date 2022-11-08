#' Standardized macroeconomic quarterly data for 8 countries.
#'
#' @source CEIC
#' @format The list of multivariate time series.
#' @importFrom 'Rdpack' reprompt
#' @references
#'     \insertRef{xu2021testing}{eigTest}
#'     \insertRef{ceic}{eigTest}
#'     \insertRef{oecd}{eigTest}
#' @examples
#' \dontrun{
#'  countryMacro
#' }
"countryMacro"

#' Coefficient matrices of 8 countries' VAR(1) models.
#'
#' @format The array of coefficient matrices for country-wise macroeconomic VAR(1) models.
#' @importFrom 'Rdpack' reprompt
#' @references
#' \insertRef{xu2021testing}{eigTest}
#' @examples
#' \dontrun{
#'  countryCoeff
#' }
"countryCoeff"

#' Asymptotic covariance matrices of the 8 coefficient matrices' estimates.
#'
#' @format The array of covariance matrices of the estimated VAR(1) coefficient matrices.
#' @importFrom 'Rdpack' reprompt
#' @references
#' \insertRef{xu2021testing}{eigTest}
#' @examples
#' \dontrun{
#'  countryCovar
#' }
"countryCovar"

#' Hudson river daily discharge table.
#'
#' @format The table of the Hudson river daily discharge with daily quantile and categories.
#' @importFrom 'Rdpack' reprompt
#' @references
#'     \insertRef{xu2021testing}{eigTest}
#'     \insertRef{usgs}{eigTest}
#' @examples
#' \dontrun{
#'  hudsonDaily
#' }
"hudsonDaily"

#' Hudson river weekly discharge table.
#'
#' @format The table of the Hudson river weekly discharge with weekly quantile and categories.
#' @importFrom 'Rdpack' reprompt
#' @references
#'     \insertRef{xu2021testing}{eigTest}
#'     \insertRef{usgs}{eigTest}
#' @examples
#' \dontrun{
#'  hudsonWeekly
#' }
"hudsonWeekly"
