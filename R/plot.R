utils::globalVariables(c("x", "y", "xend", "yend", "label", "dimension", "r", "lab",
                         "space1", "space2", "n_axes", "mAD", "value", "axis", "group",
                         "level", "xr", "yr", "residual", "col", "species", "n"))

.dim_palette <- function(dims) {
  base <- c("#2E7D32", "#1565C0", "#C62828", "#6A1B9A", "#EF6C00", "#00838F",
            "#AD1457", "#4E342E", "#558B2F", "#283593")
  stats::setNames(rep(base, length.out = length(dims)), dims)
}

#' Plot a functional space
#'
#' Species in the plane of two PCoA axes, with an optional 2-D kernel density
#' of species (contours at the 50 % and 99 % quantiles, in the spirit of the
#' TPD-based figures of Toussaint et al. 2026) and trait correlation arrows
#' coloured by functional dimension.
#'
#' @param space An `fs_space`.
#' @param axes Two axes to plot.
#' @param density Draw the species density.
#' @param arrows Draw trait correlation arrows.
#' @param min_cor Arrows are drawn for traits with |r| >= `min_cor` on at
#'   least one of the plotted axes.
#' @param highlight Optional character vector of species to highlight.
#' @param color_by Optional named vector (names = species) used to colour
#'   points (e.g. IUCN category or a metric).
#' @param point_size,point_alpha Point aesthetics.
#' @param discrete Spaces where many species share identical positions (e.g.
#'   few binary or multichoice traits) are drawn as bubbles whose area is
#'   proportional to the number of species, without kernel density. `"auto"`
#'   (default) does so when the number of distinct positions is below 10 % of
#'   the species; `TRUE`/`FALSE` force the choice.
#' @param title Plot title (default: space name).
#' @return A `ggplot` object.
#' @export
plot_space <- function(space, axes = c(1, 2), density = TRUE, arrows = TRUE,
                       min_cor = 0.3, highlight = NULL, color_by = NULL,
                       point_size = 1, point_alpha = 0.5, title = NULL,
                       discrete = "auto") {
  stopifnot(inherits(space, "fs_space"))
  if (max(axes) > ncol(space$coords)) stop("Only ", ncol(space$coords), " axes available.")
  X <- space$coords[, axes, drop = FALSE]
  df <- data.frame(species = rownames(X), x = X[, 1], y = X[, 2], stringsAsFactors = FALSE)
  e <- space$eigen
  xl <- sprintf("PCoA %d (%.1f%%)", axes[1], 100 * e$relative[axes[1]])
  yl <- sprintf("PCoA %d (%.1f%%)", axes[2], 100 * e$relative[axes[2]])
  if (is.null(title)) title <- space$name

  n_pos <- nrow(unique(round(df[, c("x", "y")], 8)))
  is_discrete <- if (identical(discrete, "auto")) n_pos < 0.1 * nrow(df) else isTRUE(discrete)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y))
  if (is_discrete && is.null(color_by)) {
    cnt <- stats::aggregate(list(n = rep(1L, nrow(df))),
                            by = list(x = round(df$x, 8), y = round(df$y, 8)), FUN = sum)
    p <- p + ggplot2::geom_point(data = cnt, ggplot2::aes(x = x, y = y, size = n),
                                 shape = 21, fill = "#E8A33D", colour = "#7F0000", alpha = 0.6,
                                 inherit.aes = FALSE) +
      ggplot2::scale_size_area(max_size = 14, name = "Species")
  } else if (density && nrow(df) >= 10) {
    p <- p + ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(level)),
                                      geom = "polygon", alpha = 0.35, bins = 12,
                                      show.legend = FALSE) +
      ggplot2::scale_fill_gradient(low = "#FFF3C4", high = "#B22222")
    cl <- .density_contours(df$x, df$y, probs = c(0.5, 0.99))
    if (!is.null(cl)) p <- p + ggplot2::geom_path(data = cl, ggplot2::aes(x = x, y = y, group = group),
                                                  colour = "#7F0000", linewidth = 0.4, inherit.aes = FALSE)
  }
  if (is_discrete && is.null(color_by)) {
    # bubbles already drawn
  } else if (is.null(color_by)) {
    p <- p + ggplot2::geom_point(size = point_size, alpha = point_alpha, colour = "grey25")
  } else {
    df$col <- color_by[df$species]
    p <- p + ggplot2::geom_point(data = df, ggplot2::aes(colour = col),
                                 size = point_size, alpha = point_alpha) +
      ggplot2::labs(colour = NULL)
  }
  if (!is.null(highlight)) {
    h <- df[df$species %in% highlight, ]
    p <- p + ggplot2::geom_point(data = h, shape = 21, size = point_size * 2.5,
                                 colour = "black", fill = "gold", stroke = 0.6) +
      ggplot2::geom_text(data = h, ggplot2::aes(label = species), size = 2.6,
                         vjust = -1, fontface = "italic")
  }
  if (arrows) {
    cc <- space$cor
    cn <- colnames(space$coords)[axes]
    if (all(cn %in% names(cc))) {
      a <- data.frame(label = cc$trait, dimension = cc$dimension,
                      xend = cc[[cn[1]]], yend = cc[[cn[2]]], stringsAsFactors = FALSE)
      a <- a[!is.na(a$xend) & pmax(abs(a$xend), abs(a$yend)) >= min_cor, ]
      if (nrow(a)) {
        sc <- 0.9 * min(max(abs(df$x)), max(abs(df$y))) / max(sqrt(a$xend^2 + a$yend^2))
        a$xend <- a$xend * sc; a$yend <- a$yend * sc
        pal <- .dim_palette(unique(space$cor$dimension))
        p <- p + ggplot2::geom_segment(data = a, ggplot2::aes(x = 0, y = 0, xend = xend, yend = yend,
                                                              colour = dimension),
                                       arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm")),
                                       inherit.aes = FALSE, linewidth = 0.5) +
          ggplot2::geom_text(data = a, ggplot2::aes(x = xend * 1.08, y = yend * 1.08,
                                                    label = label, colour = dimension),
                             size = 3, inherit.aes = FALSE, show.legend = FALSE)
        if (is.null(color_by)) p <- p + ggplot2::scale_colour_manual(values = pal, name = "Dimension")
      }
    }
  }
  p + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.18)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.08)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = xl, y = yl, title = title) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

# contour lines enclosing `probs` of the kernel density mass
.density_contours <- function(x, y, probs = c(0.5, 0.99), n = 100) {
  if (length(x) < 10 || stats::sd(x) == 0 || stats::sd(y) == 0) return(NULL)
  kd <- MASS::kde2d(x, y, n = n,
                    lims = c(range(x) + c(-1, 1) * diff(range(x)) * 0.15,
                             range(y) + c(-1, 1) * diff(range(y)) * 0.15))
  z <- kd$z / sum(kd$z)
  zs <- sort(as.vector(z), decreasing = TRUE)
  cs <- cumsum(zs)
  lv <- vapply(probs, function(p) zs[which(cs >= p)[1]], numeric(1))
  cl <- grDevices::contourLines(kd$x, kd$y, z, levels = lv)
  if (!length(cl)) return(NULL)
  do.call(rbind, lapply(seq_along(cl), function(i)
    data.frame(x = cl[[i]]$x, y = cl[[i]]$y, group = i)))
}

#' Heat map of pairwise Procrustes correlations
#' @param x An `fs_comparison` from [compare_spaces()].
#' @return A `ggplot` object.
#' @export
plot_comparison <- function(x) {
  stopifnot(inherits(x, "fs_comparison"))
  tab <- x$table
  nm <- rownames(x$matrix)
  df <- rbind(
    data.frame(space1 = tab$space1, space2 = tab$space2, r = tab$r,
               lab = sprintf("%.2f%s", tab$r, ifelse(tab$signif %in% c("", "ns"), "", tab$signif))),
    data.frame(space1 = nm, space2 = nm, r = NA, lab = "")
  )
  df$space1 <- factor(df$space1, levels = nm)
  df$space2 <- factor(df$space2, levels = rev(nm))
  ggplot2::ggplot(df, ggplot2::aes(x = space1, y = space2, fill = r)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = lab), size = 4) +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 1),
                                 na.value = "grey95", name = "Procrustes r") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Congruence between functional spaces",
                  subtitle = sprintf("%s Procrustes test, %d permutations; axes: %s",
                                     x$settings$method, x$settings$nrepet,
                                     paste(x$settings$n_axes, collapse = ","))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

#' Procrustes superimposition of two spaces
#'
#' Arrows link each species' position in the first space to its rotated
#' position in the second (first two axes of the superimposition); the
#' `n_label` species with the largest residuals are labelled.
#' @param residuals Output of [procrustes_residuals()].
#' @param n_label Number of species labelled.
#' @return A `ggplot` object.
#' @export
plot_procrustes <- function(residuals, n_label = 8) {
  fit <- attr(residuals, "procrustes")
  sp <- attr(residuals, "spaces")
  X <- fit$X[, 1:2, drop = FALSE]; Y <- fit$Yrot[, 1:2, drop = FALSE]
  df <- data.frame(species = rownames(fit$X), x = X[, 1], y = X[, 2],
                   xend = Y[, 1], yend = Y[, 2], stringsAsFactors = FALSE)
  df$residual <- residuals$residual[match(df$species, residuals$species)]
  top <- df[df$species %in% utils::head(residuals$species, n_label), ]
  ggplot2::ggplot(df) +
    ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = residual),
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.08, "cm")), linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y), size = 0.6) +
    ggplot2::geom_text(data = top, ggplot2::aes(x = xend, y = yend, label = species),
                       size = 2.6, fontface = "italic", vjust = -0.6, check_overlap = TRUE) +
    ggplot2::scale_colour_viridis_c(option = "magma", end = 0.9, direction = -1, name = "Residual") +
    ggplot2::labs(x = "Dimension 1", y = "Dimension 2",
                  title = sprintf("Procrustes: %s -> %s", sp[1], sp[2]),
                  subtitle = sprintf("Procrustes r = %.3f", sqrt(1 - fit$ss))) +
    ggplot2::theme_classic(base_size = 11)
}

#' Plot space quality against the number of axes
#' @param quality Output of [space_quality()], or a named list of them.
#' @return A `ggplot` object.
#' @export
plot_quality <- function(quality) {
  if (is.data.frame(quality)) quality <- list(space = quality)
  quality <- quality[!vapply(quality, is.null, logical(1))]
  if (!length(quality)) stop("No quality table to plot.")
  df <- do.call(rbind, lapply(names(quality), function(n) cbind(group = n, quality[[n]])))
  ggplot2::ggplot(df, ggplot2::aes(x = n_axes, y = mAD, colour = group)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::scale_x_continuous(breaks = unique(df$n_axes)) +
    ggplot2::labs(x = "Number of PCoA axes", y = "mAD (lower = better)", colour = "Space",
                  title = "Quality of the functional spaces") +
    ggplot2::theme_classic(base_size = 11)
}

#' Scree plot of the inertia captured by PCoA axes
#' @param spaces An `fs_spaces` object or a single `fs_space`.
#' @param n_axes Number of axes shown.
#' @return A `ggplot` object.
#' @export
plot_eigen <- function(spaces, n_axes = 10) {
  if (inherits(spaces, "fs_space")) spaces <- list(spaces)
  df <- do.call(rbind, lapply(spaces, function(s) {
    e <- utils::head(s$eigen, n_axes)
    data.frame(group = s$name, axis = seq_len(nrow(e)), value = 100 * e$relative)
  }))
  ggplot2::ggplot(df, ggplot2::aes(x = axis, y = value, colour = group)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::scale_x_continuous(breaks = seq_len(n_axes)) +
    ggplot2::labs(x = "PCoA axis", y = "Inertia (%)", colour = "Space") +
    ggplot2::theme_classic(base_size = 11)
}
