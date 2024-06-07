library(targets)
library(rlang)

tar_option_set(
  packages = c(
    "dplyr",
    "tidyr",
    "tibble",
    "purrr",
    "ggplot2",
    "stringr",
    "readr",
    "jsonlite",
    "rstan",
    "brms",
    "bridgesampling",
    "bayestestR"
  ),
  format = "qs"
)

options(mc.cores = parallel::detectCores())

optional <- function(cond, ...) if (cond) list2(...) else list()

tar_stanModel <- function(name) {
  nameChr <- as_string(ensym(name))
  fileChr <- paste0(nameChr, "StanFile")
  fileSym <- sym(fileChr)

  list2(
    tar_target_raw(fileChr, paste0("stan/", nameChr, ".stan"), format = "file"),
    tar_target_raw(paste0(nameChr, "Model"), expr(stan_model(!!fileSym))),
  )
}

tar_stanFit <- function(name, data, suffix1 = "", suffix2 = "", summary = TRUE, bridge = FALSE, ...) {
  nameChr <- as_string(ensym(name))
  fullNameChr <- paste0(nameChr, suffix1, suffix2)
  dataChr <- paste0(fullNameChr, "StanData")
  dataSym <- sym(dataChr)
  fitChr <- paste0(fullNameChr, "Fit")
  fitSym <- sym(fitChr)
  parsSym <- sym(paste0(nameChr, suffix1, "Pars"))
  dataExpr <- enexpr(data)
  args <- enexprs(...)

  list2(
    tar_target_raw(dataChr, dataExpr),
    tar_target_raw(fitChr, expr(
      sampling(
        !!sym(paste0(nameChr, "Model")),
        !!dataSym,
        !!!args
      )
    )),
    optional(
      summary,
      tar_target_raw(paste0(fullNameChr, "Summary"), expr(
        getFitSummary(!!fitSym, !!parsSym)
      ))
      # tar_target_raw(paste0(fullNameChr, "HDI"), expr(
      #   getFitHDI(!!fitSym, !!parsSym)
      # ))
    ),
    optional(
      bridge,
      tar_target_raw(paste0(fullNameChr, "Bridge"), expr(
        bridge_sampler(!!fitSym, sampling(
          !!sym(paste0(nameChr, "Model")),
          !!dataSym,
          chains = 0,
          !!!args
        ))
      ))
    )
  )
}

assign_if <- function(x, b, where, value) if (b) assign_in(x, where, value) else x

discard_if <- function(x, b, at) if (b) discard_at(x, at) else x

sameName <- function(x, y) trimws(tolower(x)) == trimws(tolower(y))

getCensored <- function(x, lb, ub) {
  case_when(
    x >= ub ~ "right",
    x <= lb ~ "left",
    .default = "none"
  )
}

getCensoredIndices <- function(columns, lb, ub) {
  flattened <- list_c(list_transpose(columns))
  iu <- which(flattened >= ub)
  il <- which(flattened <= lb)
  list(
    Nu = length(iu),
    Nl = length(il),
    iu = iu,
    il = il
  )
}

consolidateData <- function(dataDir) {
  demog <- read_csv(file.path(dataDir, "demographic.csv"))
  result <- list.files(dataDir, "prolific.*\\.json") |>
    map(\(fileName) {
      demog_ <- demog |>
        filter(`Participant id` == str_sub(fileName, 10, -6)) |>
        arrange(`Total approvals`) |>
        head(1)
      processPttData(demog_, read_json(file.path(dataDir, fileName)))
    })
  write_json(result, dataFilePath, pretty = 2, auto_unbox = T, digits = I(6))
  dataFilePath
}

getData <- function(dataFile) {
  read_json(dataFile) |>
    keep(\(p) p$passAttChecks)
}

getDemog <- function(data) {
  map(data, \(p) as_tibble(p$demog)) |>
    list_rbind()
}

getFitSummary <- function(fit, pars) {
  summary(
    fit,
    pars = pars,
    probs = c(.025, .5, .975)
  )$summary
}

getFitHDI <- function(fit, pars) {
  extract(fit, pars) |>
    as.data.frame() |>
    hdi()
}

formatPoints <- function(...) {
  vectors <- list2(...)
  maxMags <- map_dbl(vectors, \(v) max(abs(v)))
  digits <- map_int(maxMags, \(m) max(0, 3 - round(log10(m))))
  vectors |>
    transpose() |>
    map_chr(
      \(p) p |>
        imap_chr(\(n, i) ifelse(
          is.integer(n),
          formatC(n, format = "d"),
          formatC(n, digits = digits[i], format = "f")
        )) |>
        paste(collapse = " ")
    ) |>
    paste(collapse = "\n")
}

augment <- function(r, n) {
  1 - 1 / (1 + n * r / (1 - r))
}

getOutputFile <- function(...) {
  if (!dir.exists("output")) {
    dir.create("output")
  }
  paste0("output/", ..., ".txt")
}
