source("R/common.R")
source("R/expt3.R")

list2(
  if (dir.exists(rawDataPath)) {
    list2(
      tar_target(rawDataDir, rawDataPath, format = "file"),
      tar_target(donations, getDonations(rawDataDir), format = "file"),
      tar_target(dataFile, consolidateData(rawDataDir), format = "file"),
    )
  } else {
    tar_target(dataFile, dataFilePath)
  },
  #
  tar_target(data, getData(dataFile)),
  # robustness of Lambda Slider (varcov = variance-covariance)
  tar_target(varcovData, getVarcovData(data)),
  tar_target(varcovDataSmall, slice_sample(varcovData, prop = .15)),
  tar_target(trPlots, getTrPlots(varcovData)),
  tar_stanModel(varcov),
  #
  tar_stanFit(varcov, getVarcovStanData(varcovData), iter = 3000),
  tar_stanFit(
    varcov,
    suffix2 = "Small",
    getVarcovStanData(varcovDataSmall),
    summary = F,
    bridge = T,
    iter = 10000
  ),
  #
  tar_stanFit(
    varcov,
    suffix1 = "OneRho",
    getVarcovStanData(varcovData, oneRho = T),
    iter = 3000
  ),
  tar_stanFit(
    varcov,
    suffix1 = "OneRho",
    suffix2 = "Small",
    getVarcovStanData(varcovDataSmall, oneRho = T),
    summary = F,
    bridge = T,
    iter = 10000
  ),
  #
  tar_target(varcovOneRhoBF, bayes_factor(varcovSmallBridge, varcovOneRhoSmallBridge)),
  # lambda vs. social distance
  tar_target(distData, getDistData(data)),
  tar_target(distFit, getDistFit(distData)),
  tar_target(distHDI, hdi(distFit)),
  # lambda vs. donation
  tar_target(donateData, getDonateData(data)),
  tar_stanModel(donate),
  #
  purrr::map(c("Balanced", "SelfMore", "TargetMore"), tar_donate),
  # inequity aversion
  tar_target(ineqData, getIneqData(data)),
  tar_stanModel(ineq),
  tar_stanFit(ineq, getIneqStanData(ineqData)),
  tar_target(ineqPlot, getIneqPlot(ineqSummary)),
  tar_target(ineqOutput, getIneqOutput(ineqData, ineqSummary), format = "file"),
)
