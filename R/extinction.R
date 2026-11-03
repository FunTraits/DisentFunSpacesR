## ---------------------------------------------------------------------------
## Extinction scenarios: what would each functional dimension lose if
## threatened species disappeared?
## ---------------------------------------------------------------------------

## Internal: extinction probabilities from IUCN categories or numbers
.extinction_prob <- function(p, species) {
  if (is.null(names(p))) stop("'p' must be named by species.", call. = FALSE)
  v <- p[match(species, names(p))]
  if (is.numeric(v)) return(list(P = stats::setNames(as.numeric(v), species), cat = NULL))
  cat <- toupper(trimws(as.character(v)))
  P <- unname(iucn_probabilities()[cat])
  list(P = stats::setNames(as.numeric(P), species), cat = stats::setNames(cat, species))
}

## Internal: metrics of the remaining species, relative to the full set
.remaining_metrics <- function(X, keep, ha, base, D = NULL, metrics) {
  Xk <- X[keep, , drop = FALSE]; n <- nrow(Xk)
  out <- stats::setNames(rep(NA_real_, length(metrics)), metrics)
  if (n < 2) return(out)
  if ("FRic" %in% metrics) out["FRic"] <- .hull_vol(Xk[, seq_len(ha), drop = FALSE]) / base[["FRic"]]
  if ("FDis" %in% metrics) out["FDis"] <- mean(sqrt(rowSums(sweep(Xk, 2, colMeans(Xk))^2))) / base[["FDis"]]
  if ("FOri" %in% metrics) out["FOri"] <- mean(FNN::get.knn(Xk, k = 1)$nn.dist[, 1]) / base[["FOri"]]
  if ("MPD" %in% metrics) out["MPD"] <- if (!is.null(D)) {
    i <- which(keep); sum(D[i, i]) / (n * (n - 1)) / base[["MPD"]]
  } else .mean_pairwise_dist(Xk) / base[["MPD"]]
  out
}

#' Extinction scenarios in functional spaces
#'
#' Removes species in decreasing order of extinction risk and follows the
#' functional structure of the remaining species in every space, compared
#' with the removal of the same numbers of species at random. This answers
#' "which functional dimensions would lose most if threatened species went
#' extinct, and more than expected from species richness alone?".
#'
#' Scenarios:
#' * `risk` (ordered): species are removed by decreasing `priority` (default:
#'   extinction probability; ties broken at random in each replicate);
#' * `random`: the same numbers of species removed at random among the
#'   species with a probability (null model);
#' * stochastic (`stochastic = TRUE`): in each replicate, every species goes
#'   extinct with its probability `P` (independent Bernoulli draws; with IUCN
#'   categories the EDGE2 probabilities of [iucn_probabilities()], i.e. over
#'   about 50 years), compared with the loss of the same number of species
#'   drawn at random.
#'
#' Metrics are relative to the full set of species (1 = no loss): `FRic`
#' (convex hull volume on the first `hull_axes` axes), `FDis` (mean
#' distance to the centroid), `FOri` (mean nearest-neighbour distance) and
#' `MPD` (mean pairwise distance), computed on the first `n_axes` axes.
#'
#' @param spaces An `fs_spaces` object, a named list of `fs_space`, or one
#'   `fs_space`.
#' @param p Named vector (names = species): IUCN categories (`LC`, `NT`,
#'   `VU`, `EN`, `CR`, `DD`; other values = no probability) or extinction
#'   probabilities. Species without a probability are never removed.
#' @param priority Optional named numeric vector ordering the removals
#'   (largest first), e.g. FUSE or a probability from another source.
#'   Default: `p`.
#' @param max_n Largest number of species removed (default: the threatened
#'   species, VU + EN + CR, with IUCN categories; otherwise half of the
#'   species with a probability).
#' @param steps Number of removal steps from 0 to `max_n`.
#' @inheritParams group_overlap
#' @param metrics Metrics followed.
#' @param nrep Replicates (ties of the ordered scenario, random scenario,
#'   stochastic draws).
#' @param stochastic Also run the stochastic scenario?
#' @return An object of class `fs_extinction`: list with `curve` (space,
#'   scenario, n_removed, prop_removed, metric, mean, lower, upper: 2.5 and
#'   97.5% quantiles), `summary` (per space and metric at `max_n`: loss
#'   under the ordered and random scenarios, `ses` and one-sided `p` that
#'   the ordered loss is larger), `stochastic` (per space and metric: mean
#'   number of species lost, loss, loss expected for as many random
#'   species, `ses`, `p`) and `settings`.
#' @export
#' @examples
#' \donttest{
#' data(birds_traits); data(birds_trait_info)
#' tr <- drop_incomplete(prepare_traits(birds_traits, birds_trait_info, "species"))
#' sp <- build_spaces(tr, spaces = list(Diet = "Diet", Locomotion = "Locomotion"))
#' ext <- extinction_scenarios(sp, setNames(birds_traits$iucn, birds_traits$species),
#'                             hull_axes = 3, nrep = 19, steps = 5)
#' ext$summary
#' plot_extinction(ext)
#' }
extinction_scenarios <- function(spaces, p, priority = NULL, max_n = NULL, steps = 10,
                                 n_axes = NULL, hull_axes = 4,
                                 metrics = c("FRic", "FDis", "FOri", "MPD"),
                                 nrep = 99, stochastic = TRUE, seed = NULL) {
  if (inherits(spaces, "fs_space")) spaces <- list(space = spaces)
  if (!is.null(seed)) set.seed(seed)
  metrics <- match.arg(metrics, c("FRic", "FDis", "FOri", "MPD"), several.ok = TRUE)
  q <- function(v, pr) if (all(is.na(v))) NA_real_ else stats::quantile(v, pr, na.rm = TRUE, names = FALSE)
  res <- lapply(names(spaces), function(nm) {
    X <- .space_coords(spaces[[nm]], n_axes)
    sp <- rownames(X); n <- nrow(X)
    ep <- .extinction_prob(p, sp)
    P <- ep$P
    elig <- which(is.finite(P))
    if (length(elig) < 2) stop("Fewer than two species have an extinction probability.", call. = FALSE)
    pri <- if (is.null(priority)) P else stats::setNames(as.numeric(priority[sp]), sp)
    mx <- if (!is.null(max_n)) max_n else if (!is.null(ep$cat)) sum(ep$cat %in% c("VU", "EN", "CR")) else
      floor(length(elig) / 2)
    mx <- min(mx, length(elig), n - 3)
    ns <- unique(round(seq(0, mx, length.out = steps + 1)))
    ha <- min(hull_axes, ncol(X))
    D <- if ("MPD" %in% metrics && n <= getOption("DisentFunSpaces.max_n_dist_matrix", 10000))
      as.matrix(stats::dist(X)) else NULL
    base <- c(FRic = .hull_vol(X[, seq_len(ha), drop = FALSE]),
              FDis = mean(sqrt(rowSums(sweep(X, 2, colMeans(X))^2))),
              FOri = mean(FNN::get.knn(X, k = 1)$nn.dist[, 1]),
              MPD = if (!is.null(D)) sum(D) / (n * (n - 1)) else .mean_pairwise_dist(X))
    path <- function(ord) {                   # metrics along one removal order
      vapply(ns, function(k) {
        keep <- rep(TRUE, n); if (k > 0) keep[ord[seq_len(k)]] <- FALSE
        .remaining_metrics(X, keep, ha, base, D, metrics)
      }, numeric(length(metrics)))
    }
    pe <- pri[elig]
    run <- function(type) {
      sims <- replicate(nrep, {
        ord <- if (type == "risk") elig[order(-pe, stats::runif(length(elig)))] else sample(elig)
        path(ord)
      }, simplify = "array")                  # metrics x steps x nrep
      if (length(metrics) == 1) sims <- array(sims, c(1, length(ns), nrep))
      sims
    }
    sr <- run("risk"); sn <- run("random")
    curve <- do.call(rbind, lapply(list(risk = sr, random = sn), function(a) {
      do.call(rbind, lapply(seq_along(metrics), function(m) data.frame(
        n_removed = ns, prop_removed = ns / n, metric = metrics[m],
        mean = apply(a[m, , , drop = FALSE], 2, mean, na.rm = TRUE),
        lower = apply(a[m, , , drop = FALSE], 2, q, 0.025),
        upper = apply(a[m, , , drop = FALSE], 2, q, 0.975), stringsAsFactors = FALSE)))
    }))
    curve <- cbind(space = nm, scenario = rep(c("risk", "random"), each = length(ns) * length(metrics)),
                   curve, stringsAsFactors = FALSE)
    last <- length(ns)
    summ <- do.call(rbind, lapply(seq_along(metrics), function(m) {
      lr <- 1 - sr[m, last, ]; ln <- 1 - sn[m, last, ]
      data.frame(space = nm, metric = metrics[m], n_removed = ns[last], prop_removed = ns[last] / n,
                 loss_risk = mean(lr, na.rm = TRUE), loss_random = mean(ln, na.rm = TRUE),
                 ses = if (stats::sd(ln, na.rm = TRUE) > 0) (mean(lr, na.rm = TRUE) - mean(ln, na.rm = TRUE)) /
                   stats::sd(ln, na.rm = TRUE) else NA,
                 p = (sum(ln >= mean(lr, na.rm = TRUE), na.rm = TRUE) + 1) / (sum(is.finite(ln)) + 1),
                 stringsAsFactors = FALSE)
    }))
    stoch <- NULL
    if (stochastic) {
      st <- replicate(nrep, {
        lost <- elig[stats::runif(length(elig)) < P[elig]]
        lost <- utils::head(lost, n - 3)
        rnd <- sample(elig, length(lost))
        k1 <- rep(TRUE, n); k1[lost] <- FALSE
        k2 <- rep(TRUE, n); k2[rnd] <- FALSE
        c(length(lost), 1 - .remaining_metrics(X, k1, ha, base, D, metrics),
          1 - .remaining_metrics(X, k2, ha, base, D, metrics))
      })
      nl <- st[1, ]; L <- st[1 + seq_along(metrics), , drop = FALSE]
      R <- st[1 + length(metrics) + seq_along(metrics), , drop = FALSE]
      stoch <- do.call(rbind, lapply(seq_along(metrics), function(m) {
        dlt <- L[m, ] - R[m, ]
        data.frame(space = nm, metric = metrics[m], n_lost_mean = mean(nl),
                   loss = mean(L[m, ], na.rm = TRUE), loss_lower = q(L[m, ], 0.025),
                   loss_upper = q(L[m, ], 0.975), loss_random = mean(R[m, ], na.rm = TRUE),
                   ses = if (stats::sd(R[m, ], na.rm = TRUE) > 0)
                     (mean(L[m, ], na.rm = TRUE) - mean(R[m, ], na.rm = TRUE)) / stats::sd(R[m, ], na.rm = TRUE) else NA,
                   p = (sum(dlt <= 0, na.rm = TRUE) + 1) / (sum(is.finite(dlt)) + 1),
                   stringsAsFactors = FALSE)
      }))
    }
    list(curve = curve, summary = summ, stochastic = stoch)
  })
  bind <- function(w) { d <- do.call(rbind, lapply(res, `[[`, w)); if (!is.null(d)) rownames(d) <- NULL; d }
  structure(list(curve = bind("curve"), summary = bind("summary"), stochastic = bind("stochastic"),
                 settings = list(nrep = nrep, steps = steps, hull_axes = hull_axes, n_axes = n_axes,
                                 metrics = metrics, priority = if (is.null(priority)) "extinction probability" else "user")),
            class = "fs_extinction")
}

#' @export
print.fs_extinction <- function(x, ...) {
  cat(sprintf("Extinction scenarios (%d replicates): loss after removing the most threatened species\n",
              x$settings$nrep))
  print(x$summary, row.names = FALSE, digits = 3)
  invisible(x)
}

#' Plot extinction scenarios
#'
#' Relative value of a metric of [extinction_scenarios()] as species are
#' removed by decreasing extinction risk (red) or at random (grey), with
#' 95% envelopes, one panel per space.
#'
#' @param x Output of [extinction_scenarios()].
#' @param metric Metric to plot (`"FRic"`, `"FDis"`, `"FOri"`, `"MPD"`).
#' @return A ggplot.
#' @export
plot_extinction <- function(x, metric = "FRic") {
  d <- x$curve[x$curve$metric == metric, ]
  if (!nrow(d)) stop("Metric '", metric, "' not available.", call. = FALSE)
  d$space <- factor(d$space, levels = unique(d$space))
  d$scenario <- factor(d$scenario, levels = c("risk", "random"),
                       labels = c("by extinction risk", "at random"))
  ggplot2::ggplot(d, ggplot2::aes(n_removed, mean, colour = scenario, fill = scenario)) +
    ggplot2::geom_hline(yintercept = 1, colour = "grey70") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~space) +
    ggplot2::scale_colour_manual(values = c("by extinction risk" = "#B2182B", "at random" = "grey40"), name = NULL) +
    ggplot2::scale_fill_manual(values = c("by extinction risk" = "#B2182B", "at random" = "grey40"), name = NULL) +
    ggplot2::labs(x = "Species removed", y = paste(metric, "(relative to all species)")) +
    ggplot2::theme_classic(base_size = 11) + ggplot2::theme(legend.position = "bottom")
}
