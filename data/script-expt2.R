source("R/common.R")
source("R/expt2.R")

list2(
  if (dir.exists(rawDataPath)) {
    list2(
      tar_target(rawDataDir, rawDataPath, format = "file"),
      tar_target(dataFile, consolidateData(rawDataDir), format = "file"),
    )
  } else {
    tar_target(dataFile, dataFilePath)
  },
  #
  tar_target(data, getData(dataFile)),
  tar_target(demog, getDemog(data)),
  # H_\lambda vs H_\chi
  tar_target(chiData, getChiData(data)),
  tar_stanModel(chi),
  tar_stanFit(
    chi,
    suffix2 = "Chi",
    getChiStanData(chiData),
    bridge = T,
    iter = 5000
  ),
  tar_stanFit(
    chi,
    suffix2 = "Lambda",
    getChiStanData(chiData, muOffset = .1875),
    bridge = T,
    iter = 5000
  ),
  tar_target(chiBF, bayes_factor(chiLambdaBridge, chiChiBridge, log = T)),
  tar_target(chiPlots, getChiPlots(chiData, chiChiSummary, chiLambdaSummary)),
  tar_target(chiOutput, getChiOutput(chiData), format = "file"),
  # lambda vs. social distance
  tar_target(distFit, getDistFit(chiData)),
  tar_target(distHDI, hdi(distFit)),
)
