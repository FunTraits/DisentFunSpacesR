# Look of the Shiny application ----------------------------------------------
# One bslib theme plus a stylesheet and a script shipped in inst/app/www (the
# IBM Plex fonts are there too, so that the app looks the same offline). The
# family look of the FunTraits applications (intraitR, Rfishmorph): white
# cards on a pale page, hairlines, condensed headings.

.APP_PALETTE <- c(
  primary = "#1f4e79",   # the accent (buttons, links, current step)
  ink     = "#1a2a2f",   # text
  muted   = "#566a70",   # secondary text
  line    = "#d3dbe0",   # hairlines
  wash    = "#eef2f5",   # the page behind the cards
  gold    = "#c8901a"    # the rare second accent (active tab, focal species)
)

.app_theme <- function() {
  pal <- .APP_PALETTE
  sans <- "'IBM Plex Sans', system-ui, -apple-system, 'Segoe UI', 'Helvetica Neue', Arial, sans-serif"
  bslib::bs_theme(
    version = 5,
    bg = pal[["wash"]], fg = pal[["ink"]],
    primary = pal[["primary"]], secondary = pal[["muted"]],
    success = "#2f7a4f", info = pal[["primary"]], warning = pal[["gold"]], danger = "#a83a2c",
    base_font = sans,
    heading_font = paste0("'IBM Plex Sans Condensed', ", sans),
    "font-size-base" = "0.9375rem",
    "headings-font-weight" = 600,
    "border-radius" = "0.375rem",
    "border-radius-lg" = "0.5rem",
    "border-color" = pal[["line"]],
    "card-bg" = "#ffffff",
    "card-border-color" = pal[["line"]],
    "card-cap-bg" = "#ffffff",
    "card-box-shadow" = "none",
    "input-bg" = "#ffffff",
    "input-border-color" = "#b9c6ce",
    "link-decoration" = "none",
    "link-hover-decoration" = "underline",
    "table-border-color" = pal[["line"]])
}

.app_deps <- function() {
  htmltools::htmlDependency(
    "DisentFunSpaces-app", as.character(utils::packageVersion("DisentFunSpaces")),
    src = c(file = system.file("app", "www", package = "DisentFunSpaces")),
    stylesheet = "dfs.css", script = "dfs.js", all_files = TRUE)
}

# A navbar entry: step number + label (the number turns green when done).
.step_title <- function(n, label)
  shiny::tagList(shiny::tags$span(class = "step-n", `data-step` = n, n), label)

# The brand mark: three dots, one per functional dimension.
.brand_mark <- function() {
  shiny::HTML(paste0(
    '<svg class="brand-mark" viewBox="0 0 26 26" aria-hidden="true">',
    '<circle cx="8" cy="9" r="4.2" fill="#2E7D32"/><circle cx="18" cy="9" r="4.2" fill="#1565C0"/>',
    '<circle cx="13" cy="18" r="4.2" fill="#C62828"/>',
    '<circle cx="13" cy="12.2" r="2.1" fill="#c8901a" stroke="#fff" stroke-width="1"/></svg>'))
}

# Home-page illustration: three disaggregated spaces (one per dimension) and
# the combined space; the same focal species, linked across them. Positions
# come from a fixed pseudo-random sequence (no use of the R random stream).
.home_art_svg <- function() {
  u <- function(i) { v <- sin(i * 12.9898 + 4.1414) * 43758.5453; v - floor(v) }
  cloud <- function(n, seed, sx, sy, rho) {
    a <- u(seq_len(n) + seed); b <- u(seq_len(n) + seed + 500)
    z1 <- sqrt(-2 * log(pmax(a, 1e-3))) * cos(2 * pi * b)
    z2 <- sqrt(-2 * log(pmax(a, 1e-3))) * sin(2 * pi * b)
    cbind(sx * z1, sy * (rho * z1 + sqrt(1 - rho^2) * z2))
  }
  panels <- list(
    list(lab = "Locomotion", col = "#2E7D32", x = 8,   y = 22, w = 128, h = 92, seed = 11, rho = 0.6, focal = c(0.78, 0.22)),
    list(lab = "Diet",       col = "#1565C0", x = 156, y = 22, w = 128, h = 92, seed = 37, rho = -0.3, focal = c(0.22, 0.30)),
    list(lab = "Reproduction", col = "#C62828", x = 304, y = 22, w = 128, h = 92, seed = 73, rho = 0.1, focal = c(0.70, 0.75)))
  comb <- list(lab = "Combined", col = "#1f4e79", x = 96, y = 158, w = 248, h = 96, seed = 97, rho = 0.35,
               focal = c(0.86, 0.26))
  delay <- 0.2
  draw <- function(P, n) {
    xy <- cloud(n, P$seed, 0.16, 0.16, P$rho)
    cx <- P$x + P$w / 2; cy <- P$y + P$h / 2
    px <- pmin(pmax(cx + xy[, 1] * P$w, P$x + 6), P$x + P$w - 6)
    py <- pmin(pmax(cy + xy[, 2] * P$h, P$y + 6), P$y + P$h - 6)
    ang <- atan(P$rho) * 180 / pi * -1
    pts <- vapply(seq_len(n), function(i) {
      delay <<- delay + 0.012
      sprintf('<circle class="sp" cx="%.1f" cy="%.1f" r="2.6" fill="%s" style="animation-delay:%.2fs"/>',
              px[i], py[i], P$col, delay)
    }, "")
    fx <- P$x + P$focal[1] * P$w; fy <- P$y + P$focal[2] * P$h
    c(sprintf('<rect class="panel" x="%d" y="%d" width="%d" height="%d" rx="6"/>', P$x, P$y, P$w, P$h),
      sprintf('<text class="panel-lab" x="%d" y="%d">%s</text>', P$x + 2, P$y - 6, P$lab),
      sprintf('<ellipse class="contour" cx="%.1f" cy="%.1f" rx="%.1f" ry="%.1f" stroke="%s" transform="rotate(%.1f %.1f %.1f)"/>',
              cx, cy, P$w * 0.33, P$h * 0.27, P$col, ang, cx, cy),
      sprintf('<ellipse class="contour" cx="%.1f" cy="%.1f" rx="%.1f" ry="%.1f" stroke="%s" stroke-dasharray="2 3" transform="rotate(%.1f %.1f %.1f)"/>',
              cx, cy, P$w * 0.45, P$h * 0.42, P$col, ang, cx, cy),
      pts,
      sprintf('<circle class="focal" cx="%.1f" cy="%.1f" r="4.4"/>', fx, fy))
  }
  focal_xy <- function(P) c(P$x + P$focal[1] * P$w, P$y + P$focal[2] * P$h)
  fc <- focal_xy(comb)
  links <- vapply(seq_along(panels), function(i) {
    f <- focal_xy(panels[[i]])
    sprintf('<path class="link" style="animation-delay:%.2fs" d="M%.1f,%.1f C%.1f,%.1f %.1f,%.1f %.1f,%.1f"/>',
            1.4 + 0.15 * i, f[1], f[2], f[1], f[2] + 60, fc[1], fc[2] - 70, fc[1], fc[2])
  }, "")
  focal <- vapply(c(panels, list(comb)), function(P) {
    f <- focal_xy(P); sprintf('<circle class="focal" cx="%.1f" cy="%.1f" r="4.4"/>', f[1], f[2])
  }, "")
  body <- c(unlist(lapply(panels, draw, n = 26)), draw(comb, 40), links, focal)
  shiny::HTML(paste0(
    '<svg class="home-art" viewBox="0 0 440 262" role="img" ',
    'aria-label="Three functional spaces, one per dimension, and the combined space; ',
    'the same species highlighted in each">',
    paste(body, collapse = ""), "</svg>"))
}
