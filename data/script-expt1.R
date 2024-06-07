source("R/common.R")
source("R/expt1.R")

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
  #
  tar_target(lambdaData, getLambdaData(data)),
  tar_target(lambdaTrData, getLambdaTrData(lambdaData)),
  tar_target(svoData, getSvoData(data)),
  tar_target(svoTrData, getSvoTrData(svoData)),
  # test-retest reliability
  tar_stanModel(tr),
  tar_stanFit(
    tr,
    suffix2 = "Lambda",
    getTrStanData(
      lambdaTrData,
      lb = -2,
      ub = 2,
      muPriorMean = 0,
      muPriorSd = 2,
      sigmaPriorLogMean = log(1)
    )
  ),
  tar_stanFit(
    tr,
    suffix2 = "Svo",
    getTrStanData(
      svoTrData,
      lb = svoLb,
      ub = svoUb,
      muPriorMean = 20,
      muPriorSd = 40,
      sigmaPriorLogMean = log(20)
    )
  ),
  #
  tar_target(lambdaTrAugment, getTrAugment(3, trLambdaSummary)),
  tar_target(lambdaTrPlot, getTrPlot(lambdaTrData, trLambdaSummary)),
  tar_target(lambdaTrOutput, getTrOutput(lambdaTrData, "tr-lambda"), format = "file"),
  tar_target(svoTrPlot, getTrPlot(svoTrData, trSvoSummary)),
  tar_target(svoTrOutput, getTrOutput(svoTrData, "tr-svo"), format = "file"),
  # lambda/svo vs. social distance
  tar_target(lambdaDistFit, getLambdaDistFit(lambdaData)),
  tar_target(svoDistFit, getSvoDistFit(svoData)),
  #
  tar_target(lambdaDistPlot, getDistPlot(lambdaData, lambdaDistFit, lambda)),
  tar_target(lambdaDistOutput, getDistOutput(lambdaData, lambdaDistFit, lambda, "dist-lambda"), format = "file"),
  tar_target(svoDistPlot, getSvoDistPlot(svoData, svoDistFit)),
  tar_target(svoDistOutput, getDistOutput(svoData, svoDistFit, svo, "dist-svo"), format = "file"),
  # lambda vs svo
  tar_stanModel(convergent),
  tar_stanFit(convergent, getConvergentStanData(lambdaData, svoData)),
  tar_target(convergentPlot, getConvergentPlot(lambdaData, svoData, convergentSummary)),
  tar_target(convergentOutput, getConvergentOutput(lambdaData, svoData), format = "file"),
)
