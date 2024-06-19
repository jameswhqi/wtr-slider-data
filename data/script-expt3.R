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
  tar_target(demog, getDemog(data)),
  # robustness of Lambda Slider (varcov = variance-covariance)
  tar_target(varcovData, getVarcovData(data)),
  tar_target(trPlots, getTrPlots(varcovData)),
  tar_stanModel(varcov),
  #
  tar_stanFit(varcov, getVarcovStanData(varcovData), iter = 3000),
  tar_target(varcovOneRhoBF, getVarcovOneRhoBF(varcovModel, varcovData)),
  # lambda vs. social distance
  tar_target(distData, getDistData(data)),
  tar_target(distFit, getDistFit(distData)),
  tar_target(distHDI, hdi(distFit)),
  # lambda vs. social distance + sex
  tar_target(distSexData, getDistSexData(distData)),
  tar_target(distSexFit, getDistSexFit(distSexData)),
  tar_target(distSexNull1Fit, getDistSexNull1Fit(distSexData)),
  tar_target(distSexNull2Fit, getDistSexNull2Fit(distSexData)),
  tar_target(distSexBF1, bayes_factor(distSexFit, distSexNull1Fit)),
  tar_target(distSexBF2, bayes_factor(distSexNull1Fit, distSexNull2Fit)),
  # lambda vs. donation
  tar_target(donateData, getDonateData(data)),
  tar_stanModel(donate),
  #
  purrr::map(c("Balanced", "SelfMore", "TargetMore"), tar_donate),
  # donation vs. sex
  tar_target(donateSexData, getDonateSexData(donateData)),
  # inequity aversion
  tar_target(ineqData, getIneqData(data)),
  tar_stanModel(ineq),
  tar_stanFit(ineq, getIneqStanData(ineqData)),
  tar_target(ineqPlot, getIneqPlot(ineqSummary)),
  tar_target(ineqOutput, getIneqOutput(ineqData, ineqSummary), format = "file"),
  tar_stanModel(cor),
  tar_stanFit(cor, getIneqDonateStanData(data, ineqSummary), suffix1 = "IneqDonate"),
)
