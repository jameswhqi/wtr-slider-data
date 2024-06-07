rawDataPath <- "../expt3/data/1"

dataFilePath <- "expt2.json"

processPttData <- function(demog, json) {
  trials <- json$trials |>
    imap(
      \(trial, i) trial |>
        assign_in("trialNumber", i) |>
        assign_if(
          trial$memory,
          "memoryCorrect",
          sameName(trial$name, trial$recalledName)
        ) |>
        assign_in("targetRank", match(trial$name, json$orderedNames)) |>
        discard_if(trial$kind != "catch", c("name", "recalledName"))
    )

  memory <- trials |>
    keep(\(trial) trial$memory) |>
    map_lgl(\(trial) trial$memoryCorrect) |>
    sum()
  left <- trials |>
    keep(\(trial) trial$kind == "catch" && trial$name == "Left") |>
    map_lgl(\(trial) trial$lambda < -1.9) |>
    sum()
  right <- trials |>
    keep(\(trial) trial$kind == "catch" && trial$name == "Right") |>
    map_lgl(\(trial) trial$lambda > 1.9) |>
    sum()

  list(
    demog = list(
      age = demog$Age,
      sex = demog$Sex
    ),
    targetRanks = match(json$names, json$orderedNames),
    trials = trials,
    times = json$times,
    passAttChecks = memory + left + right >= 8
  )
}

getChiData <- function(data) {
  data |>
    imap(
      \(ptt, i) ptt$trials |>
        keep(\(trial) trial$kind != "catch") |>
        map(as_tibble) |>
        list_rbind() |>
        rename(slider = kind) |>
        select(targetRank, slider, lambda) |>
        mutate(measurement = 1:2, .by = c(targetRank, slider)) |>
        mutate(slider = case_match(slider, "base" ~ "b", "pos" ~ "p", "neg" ~ "n")) |>
        mutate(chi = (lambda + case_match(slider, "b" ~ 2, "p" ~ 1.25, "n" ~ 2.75)) / 4) |>
        mutate(censored = getCensored(chi, 0, 1)) |>
        mutate(pttId = i)
    ) |>
    list_rbind()
}

getChiStanData <- function(chiData, muOffset = 0) {
  data <- chiData |>
    pivot_wider(
      id_cols = c(pttId, targetRank),
      names_from = c(slider, measurement),
      names_sep = "",
      values_from = chi
    )

  columns <- with(data, list(b1, b2, p1, p2, n1, n2))
  rows <- list_transpose(columns)
  flattened <- list_c(rows)
  lb <- 0
  ub <- 1
  il <- which(flattened <= lb)
  iu <- which(flattened >= ub)
  list(
    muOffset = muOffset,
    lb = lb,
    ub = ub,
    N = length(rows),
    x = rows,
    Nu = length(iu),
    Nl = length(il),
    iu = iu,
    il = il,
    muPriorMean = 0.5,
    muPriorSd = 0.5,
    sigmaPriorLogMean = log(0.25)
  )
}

chiPars <- c("mu", "sigma", "corrMat[1,2]", "corrMat[1,3]", "corrMat[1,5]", "corrMat[3,4]", "corrMat[3,5]", "corrMat[5,6]")

getChiPlots <- function(chiData, chiSummary, lambdaSummary) {
  summ <- tibble(
    chi = chiSummary[, "mean"],
    lambda = lambdaSummary[, "mean"]
  ) |>
    t() |>
    as_tibble(rownames = "hypo", .name_repair = "minimal") |>
    set_names(c("hypo", rownames(chiSummary))) |>
    rename(bp = `corrMat[1,3]`, bn = `corrMat[1,5]`, pn = `corrMat[3,5]`) |>
    select(hypo, mu, sigma, bp, bn, pn) |>
    pivot_longer(
      bp:pn,
      names_to = "comp",
      values_to = "rho"
    ) |>
    mutate(
      s1 = str_sub(comp, 1, 1),
      s2 = str_sub(comp, 2, 2),
      offset1 = case_when(
        hypo == "lambda" & s1 == "p" ~ -.1875,
        .default = 0
      ),
      offset2 = case_when(
        hypo == "lambda" & s2 == "p" ~ -.1875,
        hypo == "lambda" & s2 == "n" ~ .1875,
        .default = 0
      )
    )

  map(
    c("bp", "bn", "pn"),
    \(comp_) {
      s1 <- str_sub(comp_, 1, 1)
      s2 <- str_sub(comp_, 2, 2)
      data <- chiData |>
        pivot_wider(
          id_cols = c(pttId, targetRank, measurement),
          names_from = slider,
          values_from = chi
        )

      layers <- map2(
        c("chi", "lambda"),
        c("red", "blue"),
        \(hypo_, color) {
          summ_ <- summ |>
            filter(hypo == hypo_, comp == comp_)

          line <- tibble(
            x1 = c(0, 1) + summ_$offset1,
            x2 = c(0, 1) + summ_$offset2
          )

          ellipse <- expand.grid(
            theta = seq(0, 2, .02),
            k = 1:2
          ) |>
            mutate(
              x1_ = k * cospi(theta) * sqrt(1 + summ_$rho),
              x2_ = k * sinpi(theta) * sqrt(1 - summ_$rho),
              x1 = summ_$mu + summ_$offset1 + (x1_ - x2_) / sqrt(2) * summ_$sigma,
              x2 = summ_$mu + summ_$offset2 + (x1_ + x2_) / sqrt(2) * summ_$sigma
            )

          list(
            geom_path(aes(x1, x2, group = k), ellipse, color = color),
            geom_line(aes(x1, x2), line, color = color)
          )
        }
      ) |>
        list_flatten()

      ggplot() +
        geom_point(aes(!!sym(s1), !!sym(s2)), data, alpha = .3) +
        layers +
        scale_x_continuous(breaks = scales::breaks_width(.5)) +
        scale_y_continuous(breaks = scales::breaks_width(.5)) +
        theme_bw()
    }
  )
}

getChiOutput <- function(chiData) {
  file <- getOutputFile("chi")
  data <- chiData |>
    pivot_wider(
      id_cols = c(pttId, targetRank, measurement),
      names_from = slider,
      values_from = chi
    )
  write_file(formatPoints(data$b, data$p, data$n), file)
  file
}

getDistFit <- function(chiData) {
  # get_prior(lambda | cens(censored) ~ mo(targetRank) + slider + (mo(targetRank) | pttId), chiData)
  brm(
    lambda | cens(censored) ~ mo(targetRank) + slider + (mo(targetRank) | pttId),
    chiData,
    prior = c(
      prior(normal(0, 0.5), class = b, coef = motargetRank),
      prior(normal(0, 1), class = b, coef = sliderp),
      prior(normal(0, 1), class = b, coef = slidern),
      prior(normal(1, 2), class = Intercept),
      prior(lognormal(0, 1), class = sd, group = pttId, coef = Intercept),
      prior(lognormal(-1, 1), class = sd, group = pttId, coef = motargetRank),
      prior(lognormal(0, 1), class = sigma)
    )
  )
}
