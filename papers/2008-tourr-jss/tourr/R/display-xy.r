#' Display tour path with a scattploer
#'
#' Animate a 2D tour path with a scatterplot.
#'
#' @param data matrix, or data frame containing numeric columns
#' @param axes position of the axes: center, bottomleft or off
#' @param center if TRUE, centers projected data to (0,0).  This pins the 
#'  center of data cloud and make it easier to focus on the changing shape
#'  rather than position.
#' @param limit limits of the projected data.  Defaults to 3 * square root
#'  of the largest eigenvalue.
#' @param col color to be plotted.  Defaults to "black"
#' @param pch size of the point to be plotted.  Defaults to 20.
#' @param ...  other arguments passed on to \code{\link{animate}} and
#'   \code{\link{display_xy}}
#' @param tour_path tour path, defaults to the grand tour
#' @aliases display_xy animate_xy
#' @usage display_xy(data, center = TRUE, axes = "center", limit = NULL, col = "black", pch  = 20, ...)
#'        animate_xy(data, tour_path = grand_tour(), ...)
#' @examples
#' animate_xy(flea[, 1:6])
#' animate(flea[, 1:6], grand_tour(), display_xy())
#' animate(flea[, 1:6], grand_tour(), display_xy(axes = "bottomleft"))
#' animate(flea[, 1:6], grand_tour(), display_xy(limits = c(-3, 3)))
#' animate_xy(flea[, 1:6], little_tour())
#' animate_xy(flea[, 1:3], guided_tour(holes), sphere = TRUE)
#' animate_xy(flea[, 1:6], center = FALSE)
#'
#' # The default axes are centered, like a biplot, but there are other options
#' animate_xy(flea[, 1:6], axes = "bottomleft")
#' animate_xy(flea[, 1:6], axes = "off")
#' animate_xy(flea[, 1:6], dependence_tour(c(1, 2, 1, 2, 1, 2)),
#'   axes = "bottomleft")
display_xy <- function(data, center = TRUE, axes = "center", limit = NULL, col = "black", pch  = 20, ...) {
  
  labels <- rng <- limit <- NULL
  init <- function(data) {
    if (is.null(limit)) {
      first_eigen <- sqrt(eigen(var(data))$values[1])
      limit <<- 3 * first_eigen
    }
    rng <<- c(-limit, limit)    
    labels <<- abbreviate(colnames(data), 3)
  }
  
  render_frame <- function() {
    par(pty = "s", mar = rep(1,4))
    blank_plot(xlim = rng, ylim = rng)
  }
  render_transition <- function() {
    rect(-limit, -limit, limit, limit, col="#FFFFFFE6", border=NA)
  }
  render_data <- function(data, proj, geodesic) {
    draw_tour_axes(proj, labels, limit, axes)

    # Render projected points
    x <- data %*% proj
    if (center) x <- scale(x, center = TRUE, scale = FALSE)    
    points(x, col = col, pch = pch)
  }
  
  list(
    init = init,
    render_frame = render_frame,
    render_transition = render_transition,
    render_data = render_data,
    render_target = nul
  )
}

# not being documented.  already aliased somewhere else
animate_xy <- function(data, tour_path = grand_tour(), ...) {
  animate(data, tour_path, display_xy(data, ...), ...)
}



#' Draw tour axes with base graphics
#' @keywords internal
draw_tour_axes <- function(proj, labels, limits, position) {
  position <- match.arg(position, c("center", "bottomleft", "off"))
  if (position == "off") return()
  
  if (position == "center") {
    axis_scale <- 2 * limits / 3
    axis_pos <- 0
  } else if (position == "bottomleft") {
    axis_scale <- limits / 6
    axis_pos <- -2/3 * limits
  }
 
  adj <- function(x) axis_pos + x * axis_scale

  segments(adj(0), adj(0), adj(proj[, 1]), adj(proj[, 2]), col="grey50")
  theta <- seq(0, 2 * pi, length = 50)
  lines(adj(cos(theta)), adj(sin(theta)), col = "grey50")
  text(adj(proj[, 1]), adj(proj[, 2]), label = labels, col = "grey50")       
}
