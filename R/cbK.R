#' Cost-based Ripley's K and L functions
#'
#' @param loc SpatialPoins* object with locations
#' @param cs raster object representing the cost surface in
#'   variable \code{band1}.
#' @param directions Number of directions in which cells are connected: 4
#'   (rook's case), 8 (queen's case), 16 (knight or king's case) or other (see
#'   \code{\link[raster]{adjacent}})
#'
#' @return
#'   A list with elements:
#'   \itemize{
#'     \item \code{dm} cost-based distance matrix
#'     \item \code{mds} the MDS fitted object
#'     \item \code{loc} the MDS locations as a point pattern (ppp) object
#'     \item \code{K} Ripley's K function
#'     \item \code{env_K} a Monte Carlo envelope for the K estimate
#'     \item \code{L} the L function
#'     \item \code{env_L} a Monte Carlo envelope for the L function
#'     \item \code{eucl} a list with the same objects for the Euclidean distance
#'   }
#' @import geoRcb
#' @import spatstat
#' @export
#'
#' @examples
#'   phase1 <- cbK(
#'     loc = roman_settlements
#'     cs  = tortosa_cs,
#'     directions = 4
#'   )
cbK <- function(loc, cs, directions = 8) {
  ## Compute the conductivity surface in raster format
  cond <- 1/cs

  ## Compute cost-based distance matrix among locations
  dm <- distmatGen(loc, cond, ret = "obs", directions = directions)

  ## MDS
  fit <- cmdscale(dm, k =2, eig = TRUE)

  ## Point-pattern object setup with rectangular encompassing area
  cb_pp <- ppp(
    x = fit$points[, 1],
    y = fit$points[, 2],
    xrange = range(fit$points[, 1]),
    yrange = range(fit$points[, 2])
  )

  K <- Kest(cb_pp)
  EK <- envelope(cb_pp, Kest, verbose = FALSE)

  L <- Lest(cb_pp)
  EL <- envelope(cb_pp, Lest, verbose = FALSE)

  ## Euclidean reference
  eu_dm <- dist(coordinates(loc))

  eu_pp <- ppp(
    x = coordinates(loc)[, 1],
    y = coordinates(loc)[, 2],
    xrange = range(coordinates(loc)[, 1]),
    yrange = range(coordinates(loc)[, 2])
  )

  K_eu <- Kest(eu_pp)
  EK_eu <- envelope(eu_pp, Kest, verbose = FALSE)

  L_eu <- Lest(eu_pp)
  EL_eu <- envelope(eu_pp, Lest, verbose = FALSE)

  objects <-
    list(
      dm = dm,
      mds = fit,
      loc = cb_pp,
      K = K,
      env_K = EK,
      L = L,
      env_L = EL,
      eucl = list(
        dm = eu_dm,
        loc = eu_pp,
        K = K_eu,
        env_K = EK_eu,
        L = L_eu,
        env_L = EL_eu
      )
    )
}