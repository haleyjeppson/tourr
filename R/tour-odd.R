#' An outlier tour path.
#'
#' This method generates a sample of target bases by randomly sampling
#' on the space of all d-dimensional planes in p-space and selecting
#' the basis most likely to be an outlier.
#'
#' Usually, you will not call this function directly, but will pass it to
#' a method that works with tour paths like \code{\link{animate}},
#' \code{\link{save_history}} or \code{\link{render}}.
#'
#' @param d target dimensionality
#' @param ... arguments sent to the generator
#' @export
#' @examples
#' # All animation methods use the grand tour path by default
#' animate_xy(flea[, 1:6], odd_tour(n_points = 5), col = flea$species)
odd_tour <- function(d = 2, n_points = 10, n_hist = 3, ...) {
  generator <- function(current, data, ...) {
    if (is.null(current)) {
      current <- basis_init_odd(ncol(data), d, data, n_points, n_hist)
      return(current)
    }

    # target <- basis_random(ncol(data), d)
    target <- basis_random_odd(ncol(data), d, data, n_points, n_hist)
    list(target = target)
  }

  new_geodesic_path("grand", generator)
}

#' Generate initial basis.
#'
#' First two variables are projected on first two axes.
#'
#' @keywords internal
#' @param n dimensionality of data
#' @param d dimensionality of target projection
#' @export
basis_init_odd <- function(n = ncol(data), d = 2, data, n_points = n_points, n_hist = n_hist, ...) {
  # browser()
  n <- ncol(data)
  d <- d
  rand_bases <- tibble(id = 1:(n_points*n_hist)) %>%
    mutate(bases = map(id, ~basis_random(n = ncol(data), d = d)),
           projections = map(bases, calc_proj, data),
           scags = map(projections, get_scags)) %>%
    unnest(scags) %>%
    spread(scags_type, scags)

  pc <- oddstream::get_pc_space(rand_bases[,4:12], robust = FALSE)
  # check for outliers in training set
  proj_outliers_stray <- stray::find_HDoutliers(pc$pcnorm, alpha = 0.1)
  hscv <- ks::Hscv(x = pc$pcnorm)

  current <- rand_bases %>%
    filter(id == which.max(proj_outliers_stray$out_score)) %>%
    select(bases) %>% unnest(cols = c(bases)) %>% pull(bases)

  rcd_env <- parent.frame(n = 4)
  rcd_env[["record_odd"]] <- dplyr::add_row(
    rcd_env[["record_odd"]],
    basis = list(current),
    scags = list(rand_bases[,4:12]),
    info = "new_basis",
    pc = list(pc),
    hscv = list(hscv)
  )
  return(current)
}


#' Generate a sample of random bases
#'
#' @keywords internal
#' @param n dimensionality of data
#' @param d dimensionality of target projection
#' @param n_points number of projections to generate
#' @export
basis_random_odd <- function(n = ncol(data), d = 2, data, n_points = n_points, n_hist = n_hist) {
  # browser()
  rcd_env <- parent.frame(n = 4)
  train <- tail(rcd_env[["record_odd"]], 1)
  rand_bases <- tibble(id = 1:n_points) %>%
    mutate(bases = map(id, ~basis_random(n = n,  d = d)),
           projections = map(bases, calc_proj, data),
           scags = map(projections, get_scags)) %>%
    unnest(scags) %>%
    spread(scags_type, scags)


  pctest <- scale(rand_bases[,4:12], train$pc[[1]]$center, train$pc[[1]]$scale) %*%
    train$pc[[1]]$rotation
  values <- ks::kde(x = train$pc[[1]]$pcnorm, H = train$hscv[[1]], compute.cont = TRUE,
                    eval.points = pctest[,1:2])$estimate

  current <- rand_bases %>%
    filter(id == which.min(values)) %>%
    select(bases) %>% unnest(cols = c(bases)) %>% pull(bases)

  scags <- train$scags[[1]][(n_points+1):(n_points*n_hist), ]
  scags <- rbind(scags, rand_bases[,4:12])

  pc <- oddstream::get_pc_space(scags, robust = FALSE)
  hscv <- ks::Hscv(x = pc$pcnorm)

  rcd_env[["record_odd"]] <- dplyr::add_row(
    rcd_env[["record_odd"]],
    basis = list(current),
    scags = list(scags),
    info = "new_basis",
    pc = list(pc),
    hscv = list(hscv)
  )
  return(current)
}

#' Calculate projection
#'
#' First
#'
#' @keywords internal
#' @param base projection basis
#' @param data data matrix
#' @export
calc_proj <- function(base, data){
  data %*% base
}


#' Calculate scagnostics
#'
#' First
#'
#' @keywords internal
#' @param data the projection data frame
#' @export
get_scags <- function(data){
  scags = scagnostics::scagnostics(data)
  scags_type = names(scags)
  tibble(scags = scags, scags_type = scags_type)
}
