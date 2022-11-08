#' Macroeconomic Data
#'
#' Standardized Macroeconomic Quarterly Data for 8 Countries.
#' Available from \insertCite{ceic,oecd;textual}{eigTest}.
#' See \insertCite{xu2021testing;textual}{eigTest}.
#'
#' @source CEIC
#' @format The list of multivariate time series.
#' @importFrom 'Rdpack' reprompt
#' @references
#' \insertAllCited{}
#' @examples
#' \dontrun{
#'  countryMacro
#' }
"countryMacro"

#' Macroeconomic VAR Model Data
#'
#' Estimated coefficient matrices of 8 countries' VAR(1) models.
#' See \insertCite{xu2021testing;textual}{eigTest}.
#'
#' @format The array of coefficient matrices for country-wise macroeconomic VAR(1) models.
#' @importFrom 'Rdpack' reprompt
#' @references
#' \insertAllCited{}
#' @examples
#' \dontrun{
#'  countryCoeff
#' }
"countryCoeff"

#' Macroeconomic VAR Model Data
#'
#' Asymptotic covariance matrices of the 8 coefficient matrices' estimates.
#' See \insertCite{xu2021testing;textual}{eigTest}.
#'
#' @format The array of covariance matrices of the estimated VAR(1) coefficient matrices.
#' @importFrom 'Rdpack' reprompt
#' @references
#' \insertAllCited{}
#' @examples
#' \dontrun{
#'  countryCovar
#' }
"countryCovar"

#' Streamflow Discharge Data
#'
#' Hudson river daily discharge table available from \insertCite{usgs}{eigTest}.
#' See \insertCite{xu2021testing;textual}{eigTest}.
#'
#' @format The table of the Hudson river daily discharge with daily quantile and categories.
#' @importFrom 'Rdpack' reprompt
#' @references
#' \insertAllCited{}
#' @examples
#' \dontrun{
#'  hudsonDaily
#' }
"hudsonDaily"

#' Streamflow Discharge Data
#'
#' Hudson river weekly discharge table available from \insertCite{usgs}{eigTest}.
#' See \insertCite{xu2021testing;textual}{eigTest}.
#'
#' @format The table of the Hudson river weekly discharge with weekly quantile and categories.
#' @importFrom 'Rdpack' reprompt
#' @references
#' \insertAllCited{}
#' @examples
#' \dontrun{
#'  hudsonWeekly
#' }
"hudsonWeekly"
