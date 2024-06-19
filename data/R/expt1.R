rawDataPath <- "../expt1/data/2"

dataFilePath <- "expt1.json"

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

getLambdaData <- function(data) {
  data |>
    imap(
      \(ptt, i) ptt$trials |>
        keep(\(trial) trial$kind == "lambda") |>
        map(as_tibble) |>
        list_rbind() |>
        select(targetRank, lambda) |>
        mutate(measurement = 1:2, .by = targetRank) |>
        mutate(censored = getCensored(lambda, -2, 2)) |>
        mutate(pttId = i)
    ) |>
    list_rbind()
}

getLambdaTrData <- function(lambdaData) {
  lambdaData |>
    pivot_wider(
      id_cols = c(pttId, targetRank),
      names_from = measurement,
      names_prefix = "x",
      values_from = lambda
    )
}

svoLb <- -16.26
svoUb <- 61.389

getSvoData <- function(data) {
  data |>
    imap(
      \(ptt, i) ptt$trials |>
        keep(\(trial) trial$kind == "svo") |>
        map(as_tibble) |>
        list_rbind() |>
        select(targetRank, sliderId, selfPay, oppPay, sliderX) |>
        mutate(measurement = 1:2, .by = c(targetRank, sliderId)) |>
        summarize(
          svo = atan((mean(oppPay) - 50) / (mean(selfPay) - 50)) * 180 / pi,
          .by = c(targetRank, measurement)
        ) |>
        mutate(censored = getCensored(svo, svoLb, svoUb)) |>
        mutate(pttId = i)
    ) |>
    list_rbind()
}

getSvoTrData <- function(svoData) {
  svoData |>
    pivot_wider(
      id_cols = c(pttId, targetRank),
      names_from = measurement,
      names_prefix = "x",
      values_from = svo
    )
}

getTrStanData <- function(data, lb, ub, muPriorMean, muPriorSd, sigmaPriorLogMean) {
  columns <- with(data, list(x1, x2))
  rows <- list_transpose(columns)
  flattened <- list_c(rows)
  iu <- which(flattened >= ub)
  il <- which(flattened <= lb)
  list(
    ub = ub,
    lb = lb,
    N = length(rows),
    x = rows,
    Nu = length(iu),
    Nl = length(il),
    iu = iu,
    il = il,
    muPriorMean = muPriorMean,
    muPriorSd = muPriorSd,
    sigmaPriorLogMean = sigmaPriorLogMean
  )
}

trPars <- c("mu", "sigma", "a")

getTrPlot <- function(data, fitSummary) {
  means <- fitSummary[, "mean"]
  mu <- means["mu"]
  sigma <- means["sigma"]
  rho <- means["a"]

  ellipse <- expand.grid(
    theta = seq(0, 2, .02),
    k = 1:2
  ) |>
    mutate(
      x1_ = k * cospi(theta) * sqrt(1 + rho),
      x2_ = k * sinpi(theta) * sqrt(1 - rho),
      x1 = mu + (x1_ - x2_) / sqrt(2) * sigma,
      x2 = mu + (x1_ + x2_) / sqrt(2) * sigma
    )

  ggplot() +
    geom_point(aes(x1, x2), data, alpha = .5) +
    geom_path(aes(x1, x2, group = k), ellipse) +
    theme_bw()
}

getLambdaDistFit <- function(lambdaData) {
  # get_prior(lambda | cens(censored) ~ mo(targetRank) + (mo(targetRank) | pttId), lambdaData)

  brm(
    lambda | cens(censored) ~ mo(targetRank) + (mo(targetRank) | pttId),
    lambdaData |> mutate(targetRank = ordered(targetRank)),
    prior = c(
      prior(normal(0, 0.5), class = b, coef = motargetRank),
      prior(normal(1, 2), class = Intercept),
      prior(lognormal(0, 1), class = sd, group = pttId, coef = Intercept),
      prior(lognormal(-1, 1), class = sd, group = pttId, coef = motargetRank),
      prior(lognormal(0, 1), class = sigma)
    )
  )
}

getSvoDistFit <- function(svoData) {
  # get_prior(svo | cens(censored) ~ mo(targetRank) + (mo(targetRank) | pttId), svoData)

  brm(
    svo | cens(censored) ~ mo(targetRank) + (mo(targetRank) | pttId),
    svoData |> mutate(targetRank = ordered(targetRank)),
    prior = c(
      prior(normal(0, 30), class = b, coef = motargetRank),
      prior(normal(40, 40), class = Intercept),
      prior(lognormal(3, 1), class = sd, group = pttId, coef = Intercept),
      prior(lognormal(2, 1), class = sd, group = pttId, coef = motargetRank),
      prior(lognormal(3, 1), class = sigma)
    ),
    init = \() list2(
      Intercept = 40,
      sigma = 20,
      `sd_1[1]` = 20,
      `sd_1[2]` = 7,
    )
  )
}

getDistPlot <- function(data, fit, response) {
  respSym <- ensym(response)

  p <- ggplot(data, aes(targetRank, !!respSym)) +
    geom_jitter(width = .3, height = 0, alpha = .3) +
    stat_summary() +
    theme_bw()
  eff <- conditional_effects(fit)$targetRank |>
    mutate(targetRank = as.integer(as.character(targetRank))) |>
    select(!(!!respSym)) |>
    rename(!!respSym := estimate__)
  p +
    geom_smooth(aes(ymin = lower__, ymax = upper__), eff, stat = "identity")
}

getSvoDistPlot <- function(data, fit) {
  p <- ggplot(data, aes(targetRank, svo)) +
    geom_jitter(width = 1, height = 0, alpha = .3) +
    stat_summary() +
    theme_bw()
  eff <- conditional_effects(fit)$targetRank |>
    mutate(targetRank = as.integer(as.character(targetRank))) |>
    select(!svo) |>
    rename(svo = estimate__)
  p +
    geom_smooth(aes(ymin = lower__, ymax = upper__), eff, stat = "identity")
}

getConvergentStanData <- function(lambdaData, svoData) {
  data <- inner_join(
    lambdaData |>
      pivot_wider(
        id_cols = c(pttId, targetRank),
        names_from = measurement,
        names_prefix = "l",
        values_from = lambda
      ),
    svoData |>
      pivot_wider(
        id_cols = c(pttId, targetRank),
        names_from = measurement,
        names_prefix = "s",
        values_from = svo
      ),
    by = c("pttId", "targetRank")
  )

  columns <- with(data, list(l1, l2, s1, s2))
  rows <- list_transpose(columns)

  iL <- getCensoredIndices(with(data, list(l1, l2)), -2, 2)
  iS <- getCensoredIndices(with(data, list(s1, s2)), svoLb, svoUb)

  list(
    lb = c(-2, svoLb),
    ub = c(2, svoUb),
    N = length(rows),
    x = rows,
    Nu = c(iL$Nu, iS$Nu),
    Nl = c(iL$Nl, iS$Nl),
    iu1 = iL$iu,
    il1 = iL$il,
    iu2 = iS$iu,
    il2 = iS$il,
    muPriorMean = c(0, 20),
    muPriorSd = c(2, 40),
    sigmaPriorLogMean = c(log(1), log(20))
  )
}

convergentPars <- c("mu", "sigma", "corrMat[1,2]", "corrMat[1,3]", "corrMat[3,4]")

getConvergentPlot <- function(lambdaData, svoData, fitSummary) {
  means <- fitSummary[, "mean"]
  mu1 <- means["mu[1]"]
  mu2 <- means["mu[2]"]
  sigma1 <- means["sigma[1]"]
  sigma2 <- means["sigma[2]"]
  rho <- means["corrMat[1,3]"]

  ellipse <- expand.grid(
    theta = seq(0, 2, .02),
    k = 1:2
  ) |>
    mutate(
      x1_ = k * cospi(theta) * sqrt(1 + rho),
      x2_ = k * sinpi(theta) * sqrt(1 - rho),
      x1 = mu1 + (x1_ - x2_) / sqrt(2) * sigma1,
      x2 = mu2 + (x1_ + x2_) / sqrt(2) * sigma2
    )

  data <- inner_join(
    lambdaData,
    svoData,
    by = c("pttId", "targetRank", "measurement")
  )

  ggplot() +
    geom_point(aes(lambda, svo), data, alpha = .5) +
    geom_path(aes(x1, x2, group = k), ellipse) +
    theme_bw()
}

getTrOutput <- function(data, name) {
  file <- getOutputFile(name)
  write_file(formatPoints(data$x1, data$x2), file)
  file
}

getDistOutput <- function(data, fit, response, name) {
  respSym <- ensym(response)
  respStr <- as_string(respSym)

  rawFile <- getOutputFile(name, "-raw")
  seFile <- getOutputFile(name, "-se")
  effFile <- getOutputFile(name, "-eff")

  write_file(formatPoints(data$targetRank, data[[respStr]]), rawFile)

  seData <- data |>
    summarize(mean_se(!!respSym), .by = targetRank) |>
    mutate(se = y - ymin)
  write_file(formatPoints(seData$targetRank, seData$y, seData$se), seFile)

  effData <- conditional_effects(fit)$targetRank |>
    mutate(targetRank = as.integer(as.character(targetRank)))
  write_file(with(effData, formatPoints(targetRank, estimate__, lower__, upper__)), effFile)

  c(rawFile, seFile, effFile)
}

getConvergentOutput <- function(lambdaData, svoData) {
  data <- inner_join(
    lambdaData,
    svoData,
    by = c("pttId", "targetRank", "measurement")
  )
  file <- getOutputFile("convergent")
  write_file(formatPoints(data$lambda, data$svo), file)
  file
}

getTrAugment <- function(n, fitSummary) {
  postMedian <- fitSummary["a", "50%"]
  postLow <- fitSummary["a", "2.5%"]
  postHigh <- fitSummary["a", "97.5%"]
  augment(c(postMedian, postLow, postHigh), n)
}
