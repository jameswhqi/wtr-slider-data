tar_donate <- function(slider) {
  list2(
    tar_stanFit(
      donate,
      suffix2 = slider,
      getDonateStanData(filter(donateData, slider == !!slider)),
      bridge = T,
      iter = 5000
    ),
    tar_stanFit(
      donate,
      suffix1 = "Null",
      suffix2 = slider,
      getDonateStanData(filter(donateData, slider == !!slider), null = T),
      bridge = T,
      iter = 5000
    ),
    tar_target_raw(
      paste0("donate", slider, "BF"),
      expr(bayes_factor(
        !!sym(paste0("donate", slider, "Bridge")),
        !!sym(paste0("donateNull", slider, "Bridge"))
      ))
    ),
  )
}

rawDataPath <- "../expt4/data/raw"

dataFilePath <- "expt3.json"

processPttData <- function(demog, json) {
  if (any(str_detect(json$names, regex("arica", ignore_case = T)))) {
    stop("Arica detected in ", json$client$workerId)
  }

  allNames <- c(json$names, "arica")
  trials <- list(
    json$trialConfigs,
    json$trials,
    seq_along(json$trials)
  ) |>
    pmap(
      \(config, trial, i) c(config, trial) |>
        assign_in("trialNumber", i) |>
        assign_if(
          config$kind == "Lambda" && config$memory,
          "memoryCorrect",
          sameName(allNames[[config$target + 1]], trial$recalledNames[[1]])
        ) |>
        set_names(\(name) replace(name, name == "lambdaKind", "slider")) |>
        discard_at("recalledNames")
    )

  memory <- trials |>
    keep(\(trial) trial$kind == "Lambda" && trial$memory) |>
    map_lgl(\(trial) trial$memoryCorrect) |>
    sum()
  left <- trials |>
    keep(\(trial) trial$kind == "CatchLeft") |>
    map_lgl(\(trial) trial$sliderX < .05) |>
    sum()
  right <- trials |>
    keep(\(trial) trial$kind == "CatchRight") |>
    map_lgl(\(trial) trial$sliderX > .95) |>
    sum()

  list(
    demog = list(
      age = demog$Age,
      sex = demog$Sex
    ),
    trials = trials,
    bonusSliderX = json$bonusSliderX,
    debrief = json$debrief,
    passAttChecks = memory + left + right >= 4
  )
}

getDonations <- function(dataDir) {
  fileName <- "output/donations.csv"
  list.files(dataDir, "prolific.*\\.json") |>
    map(\(fileName) {
      json <- read_json(file.path(dataDir, fileName))
      tibble(
        id = json$client$workerId |> str_sub(-5),
        donation = round(json$bonusSliderX * 2, 2)
      )
    }) |>
    list_rbind() |>
    write_csv(fileName)
  fileName
}

getVarcovData <- function(data) {
  data |>
    imap(
      \(ptt, i) ptt$trials |>
        keep(\(trial) trial$kind == "Lambda") |>
        map(as_tibble) |>
        list_rbind() |>
        mutate(measurement = 1:2, .by = c(target, slider)) |>
        mutate(slider = case_match(slider, "Balanced" ~ "b", "SelfMore" ~ "s", "TargetMore" ~ "t")) |>
        mutate(lambda = sliderX * 4 - 2) |>
        pivot_wider(
          id_cols = target,
          names_from = c(slider, measurement),
          names_sep = "",
          values_from = lambda
        ) |>
        mutate(pttId = i)
    ) |>
    list_rbind()
}

getDistData <- function(data) {
  data |>
    imap(
      \(ptt, i) ptt$trials |>
        keep(\(trial) trial$kind == "Lambda" && trial$target != 5) |>
        map(as_tibble) |>
        list_rbind() |>
        mutate(lambda = sliderX * 4 - 2) |>
        mutate(censored = getCensored(lambda, -2, 2)) |>
        mutate(pttId = i) |>
        select(target, slider, pttId, lambda, censored)
    ) |>
    list_rbind()
}

getDonateData <- function(data) {
  data |>
    imap(
      \(ptt, i) ptt$trials |>
        keep(\(trial) trial$kind == "Lambda" && trial$target == 5) |>
        map(as_tibble) |>
        list_rbind() |>
        mutate(measurement = 1:2, .by = slider) |>
        mutate(lambda = sliderX * 4 - 2) |>
        pivot_wider(
          id_cols = slider,
          names_from = measurement,
          names_prefix = "lambda",
          values_from = lambda
        ) |>
        mutate(pttId = i, donation = ptt$bonusSliderX * 2)
    ) |>
    list_rbind()
}

getTrPlots <- function(varcovData) {
  map(set_names(c("b", "s", "t")), \(slider) {
    ggplot(varcovData, aes(!!sym(paste0(slider, "1")), !!sym(paste0(slider, "2")))) +
      geom_jitter(alpha = .5, width = .1, height = .1)
  })
}

getVarcovStanData <- function(varcovData, oneRho = FALSE) {
  columns <- with(varcovData, list(b1, b2, s1, s2, t1, t2))
  rows <- list_transpose(columns)
  flattened <- list_c(rows)
  ub <- 2
  lb <- -2
  iu <- which(flattened >= ub)
  il <- which(flattened <= lb)
  list(
    oneRho = oneRho,
    ub = ub,
    lb = lb,
    N = length(rows),
    x = rows,
    Nu = length(iu),
    Nl = length(il),
    iu = iu,
    il = il,
    muPriorMean = 0,
    muPriorSd = 2,
    sigmaPriorLogMean = log(1)
  )
}

getDonateStanData <- function(donateData, null = FALSE) {
  columns <- with(donateData, list(lambda1, lambda2, donation))
  rows <- list_transpose(columns)

  iL <- getCensoredIndices(with(donateData, list(lambda1, lambda2)), -2, 2)
  iD <- getCensoredIndices(with(donateData, list(donation)), 0, 2)

  list(
    null = null,
    lb = c(-2, 0),
    ub = c(2, 2),
    N = length(rows),
    x = rows,
    Nu = c(iL$Nu, iD$Nu),
    Nl = c(iL$Nl, iD$Nl),
    iu1 = iL$iu,
    il1 = iL$il,
    iu2 = iD$iu,
    il2 = iD$il,
    muPriorMean = c(0, 1),
    muPriorSd = c(2, 1),
    sigmaPriorLogMean = c(log(1), log(1))
  )
}

varcovPars <- c("mu", "sigma", "corrMat[1,2]", "corrMat[1,3]", "corrMat[1,5]", "corrMat[3,4]", "corrMat[3,5]", "corrMat[5,6]")
varcovOneRhoPars <- c("mu", "sigma", "corrMat[1,2]", "corrMat[1,3]", "corrMat[1,5]", "corrMat[3,5]")

onecorPars <- c("mu", "sigma", "a")

donatePars <- c("mu", "sigma", "corrMat[1,2]", "corrMat[1,3]")
donateNullPars <- c("mu", "sigma", "corrMat[1,2]")

getDistFit <- function(distData) {
  brm(
    lambda | cens(censored) ~ mo(target) + slider + (mo(target) | pttId),
    distData,
    prior = c(
      prior(normal(0, 0.5), class = b, coef = motarget),
      prior(normal(0, 1), class = b, coef = sliderSelfMore),
      prior(normal(0, 1), class = b, coef = sliderTargetMore),
      prior(normal(1, 2), class = Intercept),
      prior(lognormal(0, 1), class = sd, group = pttId, coef = Intercept),
      prior(lognormal(-1, 1), class = sd, group = pttId, coef = motarget),
      prior(lognormal(0, 1), class = sigma)
    ),
    iter = 5000
  )
}

getIneqData <- function(data) {
  data |>
    imap(
      \(ptt, i) ptt$trials |>
        keep(\(trial) trial$kind == "Lambda") |>
        map(as_tibble) |>
        list_rbind() |>
        # filter(slider != "Balanced") |>
        mutate(measurement = 1:2, .by = c(target, slider)) |>
        mutate(lambda = sliderX * 4 - 2) |>
        mutate(pttId = i) |>
        select(target, slider, pttId, lambda)
    ) |>
    list_rbind()
}

getIneqStanData <- function(ineqData) {
  df <- ineqData |>
    mutate(
      pttId = as.integer(as.factor(pttId)),
      target = as.integer(as.factor(target))
    )
  list(
    Np = max(df$pttId),
    Nt = max(df$target),
    N = nrow(df),
    x = df$lambda,
    p = df$pttId,
    t = df$target,
    s = case_match(df$slider, "Balanced" ~ 1, "SelfMore" ~ 2, "TargetMore" ~ 3)
  )
}

ineqPars <- c("muL", "lambda", "kappa", "beta_", "llhd")

getIneqPlot <- function(ineqSummary) {
  df <- as_tibble(ineqSummary, rownames = "par")
  df |>
    filter(str_starts(par, "kappa")) |>
    arrange(mean) |>
    mutate(i = 1:n()) |>
    ggplot() +
    geom_segment(aes(x = `2.5%`, xend = `97.5%`, y = i, yend = i)) +
    geom_point(aes(x = mean, y = i))
}

getIneqOutput <- function(ineqData, ineqSummary) {
  kappaFile <- getOutputFile("ineq-kappa")

  summ <- as_tibble(ineqSummary, rownames = "par")
  kappa <- summ |>
    filter(str_starts(par, "kappa")) |>
    mutate(id = str_extract(par, "[0-9]+")) |>
    arrange(`50%`)
  # llhd <- summ |>
  #   filter(str_starts(par, "llhd")) |>
  #   mutate(id = str_extract(par, "[0-9]+"))
  # combined <- inner_join(kappa, llhd, by = "id", suffix = c(".k", ".l")) |>
  #   select(id, mean.k, `50%.k`, `50%.l`, `2.5%.k`, `97.5%.k`)
  lambda <- summ |>
    filter(str_starts(par, "lambda")) |>
    mutate(
      id = str_extract(par, "[0-9]+"),
      target = str_extract(par, ",([0-9]+)", group = 1)
    )

  write_file(formatPoints(kappa$`50%`, kappa$`2.5%`, kappa$`97.5%`), kappaFile)

  n <- c(3, 52, 73)
  rawFiles <- getOutputFile("ineq-raw-", n)
  predFiles <- getOutputFile("ineq-pred-", n)

  for (i in 1:3) {
    id <- kappa$id[n[i]]
    df <- ineqData |>
      filter(pttId == id) |>
      mutate(slider = case_match(slider, "Balanced" ~ 1L, "SelfMore" ~ 2L, "TargetMore" ~ 3L)) |>
      arrange(target, slider) |>
      group_by(target) |>
      mutate(x = as.integer(c(2, 3, 5, 6, 8, 9)))
    write_file(formatPoints(df$target + 1L, df$x, df$slider, df$lambda), rawFiles[i])

    lambdas <- lambda |>
      filter(id == .env$id) |>
      arrange(target) |>
      pull(`50%`)
    k <- kappa$`50%`[n[i]]
    df <- expand_grid(
      target = 1:6,
      slider = 1:3
    ) |>
      mutate(pred = map2_dbl(slider, target, \(s, t) getIneqPred(s, lambdas[t], k)))
    write_file(formatPoints(df$target, df$slider, df$pred), predFiles[i])
  }

  c(kappaFile, rawFiles, predFiles)
}

getIneqPred <- function(slider, l, k) {
  mu1 <- (l + k) / (1 - k)
  mu2 <- (l - k) / (1 + k)
  xc <- sqrt(3) - 1
  case_match(
    slider,
    1 ~ if (mu1 < -2) {
      -2
    } else if (mu1 < xc) {
      mu1
    } else if (mu2 < xc) {
      xc
    } else if (mu2 < 2) {
      mu2
    } else {
      2
    },
    2 ~ if (mu1 < -2) {
      -2
    } else if (mu1 < 2) {
      mu1
    } else {
      2
    },
    3 ~ if (mu2 < -2) {
      -2
    } else if (mu2 < 2) {
      mu2
    } else {
      2
    }
  )
}

getIneqDonateStanData <- function(data, ineqSummary) {
  donate <- data |>
    imap(
      \(ptt, i) tibble(id = as.character(i), donation = ptt$bonusSliderX * 2)
    ) |>
    list_rbind()

  summ <- as_tibble(ineqSummary, rownames = "par")
  kappa <- summ |>
    filter(str_starts(par, "kappa")) |>
    mutate(id = str_extract(par, "[0-9]+")) |>
    select(id, `50%`) |>
    rename(kappa = `50%`)

  df <- inner_join(donate, kappa)

  list(
    N = nrow(df),
    x = list_transpose(list(abs(df$donation - 1), df$kappa)),
    muPriorMean = c(.5, .5),
    muPriorSd = c(.5, .5),
    sigmaPriorLogMean = c(log(.5), log(.5))
  )
}

getIneqDonatePred <- function(l, k) {
  if (l + 2 * k - 1 < 0) {
    0
  } else if (l - 2 * k - 1 < 0) {
    .5
  } else {
    1
  }
}

corIneqDonatePars <- c("mu", "sigma", "a")
