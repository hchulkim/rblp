# Dataset documentation and loader functions for the bundled empirical datasets
# (Nevo 2000 cereal data and BLP 1995 automobile data). These datasets serve
# as standard benchmarks for validating the BLP estimator against published
# results and pyblp output.

#' Nevo (2000) Cereal Product Data
#'
#' Product-level data from the Nevo (2000) study of the US ready-to-eat cereal
#' market. Contains 2256 product-market observations across 94 markets
#' (47 cities x 2 quarters).
#'
#' @format A data frame with 2256 rows and the following columns:
#' \describe{
#'   \item{market_ids}{Market identifier (city-quarter combination)}
#'   \item{city_ids}{City identifier}
#'   \item{quarter}{Quarter (1 or 2)}
#'   \item{product_ids}{Product identifier within market}
#'   \item{firm_ids}{Firm identifier}
#'   \item{brand_ids}{Brand identifier}
#'   \item{shares}{Market share}
#'   \item{prices}{Product price (dollars per serving)}
#'   \item{sugar}{Sugar content (grams per serving)}
#'   \item{mushy}{Mushiness indicator (0/1)}
#'   \item{demand_instruments0}{Excluded demand-side instrument 0}
#'   \item{demand_instruments1}{Excluded demand-side instrument 1}
#'   \item{demand_instruments2}{Excluded demand-side instrument 2}
#'   \item{demand_instruments3}{Excluded demand-side instrument 3}
#'   \item{demand_instruments4}{Excluded demand-side instrument 4}
#'   \item{demand_instruments5}{Excluded demand-side instrument 5}
#'   \item{demand_instruments6}{Excluded demand-side instrument 6}
#'   \item{demand_instruments7}{Excluded demand-side instrument 7}
#'   \item{demand_instruments8}{Excluded demand-side instrument 8}
#'   \item{demand_instruments9}{Excluded demand-side instrument 9}
#'   \item{demand_instruments10}{Excluded demand-side instrument 10}
#'   \item{demand_instruments11}{Excluded demand-side instrument 11}
#'   \item{demand_instruments12}{Excluded demand-side instrument 12}
#'   \item{demand_instruments13}{Excluded demand-side instrument 13}
#'   \item{demand_instruments14}{Excluded demand-side instrument 14}
#'   \item{demand_instruments15}{Excluded demand-side instrument 15}
#'   \item{demand_instruments16}{Excluded demand-side instrument 16}
#'   \item{demand_instruments17}{Excluded demand-side instrument 17}
#'   \item{demand_instruments18}{Excluded demand-side instrument 18}
#'   \item{demand_instruments19}{Excluded demand-side instrument 19}
#' }
#'
#' @references
#' Nevo, A. (2000). A Practitioner's Guide to Estimation of
#' Random-Coefficients Logit Models of Demand. \emph{Journal of Economics &
#' Management Strategy}, 9(4), 513-548.
#'
#' @source \url{https://pyblp.readthedocs.io/en/stable/}
#' @examples
#' nevo_products <- load_nevo_products()
#' head(nevo_products)
"nevo_products"

#' Nevo (2000) Cereal Agent Data
#'
#' Agent-level data (simulated consumers) for the Nevo (2000) cereal study.
#' Contains integration nodes, weights, and demographic variables for
#' 20 agents per market.
#'
#' @format A data frame with 1880 rows and the following columns:
#' \describe{
#'   \item{market_ids}{Market identifier (must match nevo_products)}
#'   \item{city_ids}{City identifier}
#'   \item{quarter}{Quarter}
#'   \item{weights}{Integration weights (sum to 1 within each market)}
#'   \item{nodes0}{Integration node dimension 0}
#'   \item{nodes1}{Integration node dimension 1}
#'   \item{nodes2}{Integration node dimension 2}
#'   \item{nodes3}{Integration node dimension 3}
#'   \item{income}{Log income}
#'   \item{income_squared}{Squared log income}
#'   \item{age}{Age (transformed)}
#'   \item{child}{Child presence indicator (transformed)}
#' }
#'
#' @references
#' Nevo, A. (2000). A Practitioner's Guide to Estimation of
#' Random-Coefficients Logit Models of Demand. \emph{Journal of Economics &
#' Management Strategy}, 9(4), 513-548.
#'
#' @source \url{https://pyblp.readthedocs.io/en/stable/}
"nevo_agents"

#' BLP (1995) Automobile Product Data
#'
#' Product-level data from the Berry, Levinsohn, and Pakes (1995) study of
#' the US automobile market. Contains 2217 product-market observations across
#' 20 annual markets (1971-1990).
#'
#' @format A data frame with 2217 rows and the following columns:
#' \describe{
#'   \item{market_ids}{Market year}
#'   \item{clustering_ids}{Clustering identifier for standard errors}
#'   \item{car_ids}{Vehicle identifier}
#'   \item{firm_ids}{Manufacturer identifier}
#'   \item{region}{Region of origin (US, Japan, Europe)}
#'   \item{shares}{Market share}
#'   \item{prices}{Price (10,000s of 1983 dollars)}
#'   \item{hpwt}{Horsepower/weight ratio}
#'   \item{air}{Air conditioning indicator}
#'   \item{mpd}{Miles per dollar}
#'   \item{mpg}{Miles per gallon}
#'   \item{space}{Interior space (length x width)}
#'   \item{trend}{Time trend}
#'   \item{demand_instruments0}{Excluded demand-side instrument 0}
#'   \item{demand_instruments1}{Excluded demand-side instrument 1}
#'   \item{demand_instruments2}{Excluded demand-side instrument 2}
#'   \item{demand_instruments3}{Excluded demand-side instrument 3}
#'   \item{demand_instruments4}{Excluded demand-side instrument 4}
#'   \item{demand_instruments5}{Excluded demand-side instrument 5}
#'   \item{demand_instruments6}{Excluded demand-side instrument 6}
#'   \item{demand_instruments7}{Excluded demand-side instrument 7}
#'   \item{supply_instruments0}{Excluded supply-side instrument 0}
#'   \item{supply_instruments1}{Excluded supply-side instrument 1}
#'   \item{supply_instruments2}{Excluded supply-side instrument 2}
#'   \item{supply_instruments3}{Excluded supply-side instrument 3}
#'   \item{supply_instruments4}{Excluded supply-side instrument 4}
#'   \item{supply_instruments5}{Excluded supply-side instrument 5}
#'   \item{supply_instruments6}{Excluded supply-side instrument 6}
#'   \item{supply_instruments7}{Excluded supply-side instrument 7}
#'   \item{supply_instruments8}{Excluded supply-side instrument 8}
#'   \item{supply_instruments9}{Excluded supply-side instrument 9}
#'   \item{supply_instruments10}{Excluded supply-side instrument 10}
#'   \item{supply_instruments11}{Excluded supply-side instrument 11}
#' }
#'
#' @references
#' Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile Prices in Market
#' Equilibrium. \emph{Econometrica}, 63(4), 841-890.
#'
#' @source \url{https://pyblp.readthedocs.io/en/stable/}
"blp_products"

#' BLP (1995) Automobile Agent Data
#'
#' Agent-level data for the BLP (1995) automobile study. Contains integration
#' nodes, weights, and income demographics.
#'
#' @format A data frame with 4000 rows and the following columns:
#' \describe{
#'   \item{market_ids}{Market year}
#'   \item{weights}{Integration weights}
#'   \item{nodes0}{Integration node dimension 0}
#'   \item{nodes1}{Integration node dimension 1}
#'   \item{nodes2}{Integration node dimension 2}
#'   \item{nodes3}{Integration node dimension 3}
#'   \item{nodes4}{Integration node dimension 4}
#'   \item{income}{Simulated income}
#' }
#'
#' @references
#' Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile Prices in Market
#' Equilibrium. \emph{Econometrica}, 63(4), 841-890.
#'
#' @source \url{https://pyblp.readthedocs.io/en/stable/}
"blp_agents"

#' Load Nevo (2000) cereal product data
#'
#' @return Data frame with cereal product data
#' @export
load_nevo_products <- function() {
  find_extdata_("nevo_products.csv")
}

#' Load Nevo (2000) cereal agent data
#'
#' @return Data frame with cereal agent data
#' @export
load_nevo_agents <- function() {
  find_extdata_("nevo_agents.csv")
}

#' Load BLP (1995) automobile product data
#'
#' @return Data frame with automobile product data
#' @export
load_blp_products <- function() {
  find_extdata_("blp_products.csv")
}

#' Load BLP (1995) automobile agent data
#'
#' @return Data frame with automobile agent data
#' @export
load_blp_agents <- function() {
  find_extdata_("blp_agents.csv")
}

#' Find and load extdata file (works in both installed and development mode)
#' @param filename CSV file name
#' @return Data frame
#' @keywords internal
find_extdata_ <- function(filename) {
  # Try system.file first (installed package)
  path <- system.file("extdata", filename, package = "rblp")
  if (nzchar(path)) return(utils::read.csv(path, stringsAsFactors = FALSE))

  # Development mode: look relative to working directory and common test paths
  candidates <- c(
    file.path("inst", "extdata", filename),
    file.path("..", "inst", "extdata", filename),
    file.path("..", "..", "inst", "extdata", filename),
    file.path("extdata", filename)
  )
  for (p in candidates) {
    if (file.exists(p)) return(utils::read.csv(p, stringsAsFactors = FALSE))
  }

  stop(sprintf("Cannot find '%s'. Is the rblp package installed or are you in the package root?",
               filename))
}
