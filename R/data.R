#' Macroeconomic Data
#'
#' @description
#' \itemize{
#' \item \code{countryMacro} - Standardized Macroeconomic Quarterly Data for 8 Countries.
#' \item \code{countryCoeff} - Estimated coefficient matrices of 8 countries' VAR(1) models.
#' \item \code{countryCovar} - Estimated covariance matrices of the 8 coefficient matrices' estimates.
#' }
#' Data are accessed from \insertCite{ceic,oecd;textual}{eigTest}.
#' See \insertCite{xu2021testing;textual}{eigTest} for processing details.
#'
#' @source CEIC, OECD
#' @format `countryMacro` - The length-8 \code{list} of multivariate time series.
#'
#' @importFrom 'Rdpack' reprompt
#' @rdname MacroeconomicData
#' @references
#' \insertAllCited{}
#' @examples countryMacro; countryCoeff; countryCovar
"countryMacro"

#' Macroeconomic VAR Model Data
#'
#' @format `countryCoeff` - The 8 x 3 x 3 \code{array} of coefficient matrices for country-wise macroeconomic VAR(1) models.
#' @rdname MacroeconomicData
"countryCoeff"

#' Macroeconomic VAR Model Data
#'
#' @format `countryCovar` - The 8 x 9 x 9 \code{array} of covariance matrices of the estimated VAR(1) coefficient matrices.
#' @rdname MacroeconomicData
"countryCovar"

#' Streamflow Discharge Data
#'
#' @description Hudson river discharges with quantiles and categories, consisting 1827 rows and 6 columns.
#' \itemize{
#' \item \code{hudsonDaily} - The discharge table of daily data.
#' \item \code{hudsonWeekly} - The discharge table of weekly data.
#' }
#' Data are accessed from \insertCite{usgs;textual}{eigTest}.
#' See \insertCite{xu2021testing;textual}{eigTest} for processing details.
#'
#' @source USGS
#' @format NULL
#'
#' @importFrom 'Rdpack' reprompt
#' @rdname StreamflowDischarge
#' @references
#' \insertAllCited{}
#' @examples hudsonDaily; hudsonWeekly
"hudsonDaily"

#' Streamflow Discharge Data
#'
#' @format NULL
#' @rdname StreamflowDischarge
"hudsonWeekly"
