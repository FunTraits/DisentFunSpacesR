# Interactive figures in the Shiny application --------------------------------
# The plot_*() functions return ggplots (used as such in scripts, reports and
# the PDF downloads). In the app they are converted to 'plotly' widgets: hover
# shows the species (and its group, value...), zoom / pan, legend entries can
# be hidden by a click, and the camera button exports an SVG. 'plotly' is
# suggested, not imported: without it the app falls back to static plots.

.has_plotly <- function() requireNamespace("plotly", quietly = TRUE)

# A point layer with a hover text. ggplot2 ignores the 'text' aesthetic (and
# warns at construction); plotly::ggplotly() shows it in the tooltip.
.geom_point_hover <- function(...) suppressWarnings(ggplot2::geom_point(...))

# Output and render functions of the app: plotly when available, else static.
.fig_output <- function(id, height = "400px") {
  if (.has_plotly()) plotly::plotlyOutput(id, height = height)
  else shiny::plotOutput(id, height = height)
}

.render_fig <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) expr <- substitute(expr)
  if (.has_plotly()) {
    plotly::renderPlotly(bquote(.(.to_plotly)(.(expr))), env = env, quoted = TRUE)
  } else shiny::renderPlot(expr, env = env, quoted = TRUE, res = 96)
}

#' Convert a ggplot of the package into an interactive plotly figure
#'
#' Used by the Shiny application; it can also be called on the output of any
#' `plot_*()` function of the package. Points show the species name (and its
#' group or value) on hover; arrows (trait correlations) keep their heads; the
#' subtitle is kept under the title; the camera button exports an SVG.
#'
#' @param p A ggplot.
#' @param filename Name of the exported image (without extension).
#' @return A `plotly` htmlwidget.
#' @export
#' @examplesIf requireNamespace("plotly", quietly = TRUE)
#' sp <- build_spaces(prepare_traits(birds_traits[1:200, ], birds_trait_info),
#'                    spaces = list(Diet = "Diet"))
#' as_plotly(plot_space(sp$Diet))
as_plotly <- function(p, filename = "figure") .to_plotly(p, filename = filename)

.to_plotly <- function(p, filename = "figure") {
  if (!.has_plotly()) stop("Install the 'plotly' package for interactive figures.", call. = FALSE)
  if (inherits(p, "plotly")) return(p)
  stopifnot(inherits(p, "ggplot"))
  # Procrustes arrows: ggplotly() makes one trace per species (continuous colour
  # on segments), which blocks the app for minutes with thousands of species
  if (!is.null(attr(p, "fs_residuals")))
    return(.plotly_procrustes(attr(p, "fs_residuals"), attr(p, "fs_n_label"), filename))
  if (!is.null(attr(p, "fs_comparison")))
    return(.plotly_comparison(attr(p, "fs_comparison"), filename))
  has_text <- any(vapply(p$layers, function(l) "text" %in% names(l$mapping), logical(1)))
  # arrow heads are lost by ggplotly(): redrawn below as annotations
  arr <- .arrow_layers(p)
  g <- suppressWarnings(plotly::ggplotly(p, tooltip = if (has_text) "text" else "all"))
  g <- suppressWarnings(plotly::plotly_build(g))
  for (i in seq_along(g$x$data)) {
    tr <- g$x$data[[i]]
    # layers without a hover text (densities, hulls, contours): no empty tooltip
    if (has_text && (is.null(tr$text) || all(!nzchar(as.character(tr$text)))))
      g$x$data[[i]]$hoverinfo <- "skip"
    # "(Aquarium,outside reference space)" -> "Aquarium, outside reference space"
    if (!is.null(tr$name) && grepl("^\\(.*\\)$", tr$name)) {
      nm <- sub("^\\((.*)\\)$", "\\1", tr$name)
      nm <- sub(",1$", "", gsub(",", ", ", nm, fixed = TRUE))
      g$x$data[[i]]$name <- nm
    }
  }
  # labels above their point (ggplotly() ignores vjust)
  for (i in seq_along(g$x$data))
    if (identical(g$x$data[[i]]$mode, "text")) {
      g$x$data[[i]]$textposition <- "top center"
      g$x$data[[i]]$cliponaxis <- FALSE
    }
  # coord_equal(): same scale on both axes
  if (inherits(p$coordinates, "CoordFixed") && !is.null(g$x$layout$yaxis)) {
    g$x$layout$yaxis$scaleanchor <- "x"
    g$x$layout$yaxis$scaleratio <- if (is.null(p$coordinates$ratio)) 1 else p$coordinates$ratio
  }
  if (length(arr)) {
    ann <- c(if (!is.null(g$x$layout$annotations)) g$x$layout$annotations, arr)
    g$x$layout$annotations <- ann
  }
  .finish_plotly(g, p$labels$title, p$labels$subtitle, filename)
}

# Segments drawn with an arrow head (at most 200 per plot) -> plotly annotations
.arrow_layers <- function(p) {
  idx <- which(vapply(p$layers, function(l)
    inherits(l$geom, "GeomSegment") && !is.null(l$geom_params$arrow), logical(1)))
  if (!length(idx)) return(list())
  out <- list()
  for (i in idx) {
    d <- tryCatch(ggplot2::layer_data(p, i), error = function(e) NULL)
    if (is.null(d) || !nrow(d) || nrow(d) > 200) next
    d <- d[is.finite(d$x) & is.finite(d$xend), , drop = FALSE]
    col <- if (is.null(d$colour)) rep("black", nrow(d)) else d$colour
    rgb <- grDevices::col2rgb(col)
    col <- sprintf("rgb(%d,%d,%d)", rgb[1, ], rgb[2, ], rgb[3, ])
    out <- c(out, lapply(seq_len(nrow(d)), function(j) list(
      x = d$xend[j], y = d$yend[j], ax = d$x[j], ay = d$y[j],
      xref = "x", yref = "y", axref = "x", ayref = "y", text = "",
      showarrow = TRUE, arrowhead = 2, arrowsize = 1, arrowwidth = 1.2,
      arrowcolor = col[j], standoff = 0)))
  }
  out
}

# Title + subtitle, fonts, toolbar and SVG export of every figure
.finish_plotly <- function(g, ttl = NULL, sub = NULL, filename = "figure") {
  if (!is.null(sub)) {
    sub <- sprintf("<span style='font-size:12px;color:#566a70'>%s</span>", sub)
    txt <- if (is.null(ttl)) sub else sprintf("<b>%s</b><br>%s", ttl, sub)
    g$x$layout$title <- list(text = txt, x = 0, xref = "paper", xanchor = "left", font = list(size = 14))
    g$x$layout$margin$t <- max(c(g$x$layout$margin$t, 60), na.rm = TRUE)
  } else if (!is.null(ttl) && is.null(g$x$layout$title)) {
    g$x$layout$title <- list(text = sprintf("<b>%s</b>", ttl), x = 0, xref = "paper", xanchor = "left", font = list(size = 14))
  }
  g$x$layout$font <- utils::modifyList(
    if (is.null(g$x$layout$font)) list() else g$x$layout$font,
    list(family = "'IBM Plex Sans', system-ui, -apple-system, 'Segoe UI', Arial, sans-serif"))
  g$x$layout$hoverlabel <- list(font = list(size = 12))
  plotly::config(g, displaylogo = FALSE,
                 modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d",
                                            "hoverClosestCartesian", "hoverCompareCartesian"),
                 toImageButtonOptions = list(format = "svg", filename = filename))
}

# Procrustes superimposition built directly with plotly: segments grouped in
# classes of residuals (one trace per class), species points with hover text.
.plotly_procrustes <- function(residuals, n_label = 8, filename = "procrustes", n_class = 6) {
  if (is.null(n_label)) n_label <- 8
  fit <- attr(residuals, "procrustes"); sp <- attr(residuals, "spaces")
  X <- fit$X[, 1:2, drop = FALSE]; Y <- fit$Yrot[, 1:2, drop = FALSE]
  sps <- rownames(fit$X)
  res <- residuals$residual[match(sps, residuals$species)]
  br <- unique(stats::quantile(res, seq(0, 1, length.out = n_class + 1), na.rm = TRUE, names = FALSE))
  cl <- if (length(br) > 1) cut(res, br, include.lowest = TRUE) else factor(rep("all", length(res)))
  pal <- rev(grDevices::hcl.colors(nlevels(cl) + 1, "Inferno"))[-1]   # light (small) to dark (large)
  lab <- if (length(br) > 1) sprintf("%.3g - %.3g", utils::head(br, -1), br[-1]) else "all"
  g <- plotly::plot_ly()
  for (k in seq_len(nlevels(cl))) {
    i <- which(cl == levels(cl)[k])
    if (!length(i)) next
    xs <- as.vector(rbind(X[i, 1], Y[i, 1], NA)); ys <- as.vector(rbind(X[i, 2], Y[i, 2], NA))
    g <- plotly::add_trace(g, type = "scatter", mode = "lines", x = xs, y = ys,
                           line = list(color = pal[k], width = 1), name = lab[k],
                           legendgroup = lab[k], hoverinfo = "skip")
  }
  g <- plotly::add_trace(g, type = "scatter", mode = "markers", x = X[, 1], y = X[, 2],
                         marker = list(size = 3.5, color = "black"), name = "species", showlegend = FALSE,
                         text = sprintf("%s<br>residual = %.3g", sps, res), hoverinfo = "text")
  top <- match(utils::head(residuals$species, n_label), sps)
  top <- top[!is.na(top)]
  if (length(top)) g <- plotly::add_trace(g, type = "scatter", mode = "text", x = Y[top, 1], y = Y[top, 2],
                                          text = sps[top], textposition = "top center", cliponaxis = FALSE,
                                          textfont = list(size = 10), hoverinfo = "skip", showlegend = FALSE)
  ax <- function(t) list(title = list(text = t), zeroline = FALSE, showgrid = FALSE, showline = TRUE,
                         ticks = "outside", linecolor = "black")
  g <- plotly::layout(g, xaxis = ax("Dimension 1"), yaxis = ax("Dimension 2"),
                      legend = list(title = list(text = "Residual")), plot_bgcolor = "white")
  g <- plotly::plotly_build(g)
  .finish_plotly(g, sprintf("Procrustes: %s -> %s", sp[1], sp[2]),
                 sprintf("Procrustes r = %.3f", sqrt(1 - fit$ss)), filename)
}

# Matrix of Procrustes correlations built directly with plotly (same layout as
# plot_comparison(): lower triangle, fixed 0-1 colour scale, empty diagonal).
.plotly_comparison <- function(x, filename = "procrustes") {
  tab <- x$table; nm <- rownames(x$matrix); k <- length(nm)
  z <- matrix(NA_real_, k, k, dimnames = list(rev(nm), nm))   # rows = y (top to bottom: rev)
  lab <- hov <- matrix("", k, k)
  for (i in seq_len(nrow(tab))) {
    a <- match(tab$space1[i], nm); b <- match(tab$space2[i], rev(nm))
    z[b, a] <- tab$r[i]
    lab[b, a] <- sprintf("%.2f%s", tab$r[i], if (tab$signif[i] %in% c("", "ns")) "" else tab$signif[i])
    hov[b, a] <- sprintf("%s vs %s<br>r = %.3f<br>p = %.3g<br>%s species", tab$space1[i], tab$space2[i],
                         tab$r[i], tab$p_value[i], tab$n_species[i])
  }
  g <- plotly::plot_ly(x = nm, y = rev(nm), z = z, type = "heatmap", zmin = 0, zmax = 1,
                       colorscale = list(list(0, "#F7FBFF"), list(1, "#08519C")),
                       xgap = 1, ygap = 1, text = hov, hoverinfo = "text", hoverongaps = FALSE,
                       colorbar = list(title = list(text = "Procrustes r")))
  xx <- rep(nm, each = k); yy <- rep(rev(nm), times = k); ll <- as.vector(lab)
  keep <- nzchar(ll)
  if (any(keep)) g <- plotly::add_trace(g, type = "scatter", mode = "text", x = xx[keep], y = yy[keep],
                                        text = ll[keep], textfont = list(size = 13), hoverinfo = "skip",
                                        showlegend = FALSE, inherit = FALSE)
  ax <- function(lev) list(title = "", type = "category", categoryorder = "array", categoryarray = lev,
                           showgrid = FALSE, zeroline = FALSE)
  g <- plotly::layout(g, xaxis = c(ax(nm), list(tickangle = -30)), yaxis = ax(rev(nm)),
                      plot_bgcolor = "white")
  g <- plotly::plotly_build(g)
  .finish_plotly(g, "Congruence between functional spaces",
                 sprintf("%s Procrustes test, %d permutations; axes: %s", x$settings$method,
                         x$settings$nrepet, paste(x$settings$n_axes, collapse = ",")), filename)
}
