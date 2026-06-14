library(targets)
library(crew)
# gcloud storage rsync -r Results/ gs://dfm-regions-results/Results/
# gcloud storage rsync -r gs://dfm-regions-results/Results/ Results/

tar_option_set(
  packages = c(
    "tidyverse",
    "readxl",
    "writexl",
    "tsibble",
    "foreach",
    "seastests",
    "deseats",
    "seasonal",
    "tempdisagg",
    "zoo",
    "xts",
    "glue",
    "flextable",
    "huxtable",
    "imputeTS",
    "writexl",
    "Rcpp",
    "RcppArmadillo",
    "scoringRules",
    "Matrix",
    "coda"
  ),
  controller = crew_controller_local(
    workers = 16,
    seconds_launch = 120,
    seconds_interval = 1,
    crashes_max = 1000L,
    retry_tasks = FALSE,
    options_local = crew::crew_options_local(log_directory = "logs/crew")
  ),
  seed = 12345,
  error = "null"
)

options(scipen = 9999)

source("R/functions.R")

list(
  # ---- Configuration ---------------------------------------------------------
  tar_target(runTag, "09_04"),
  tar_target(
    cfg,
    list(
      sample_start = as.Date("2008-01-01"),
      sample_end = as.Date("2025-11-01"),
      estim_start = as.Date("2009-01-01")
    )
  ),
  tar_target(
    estSettings,
    list(
      hor = 12,
      M0 = 5000,
      M = 5000,
      L = 2,
      Lu = 2,
      T00 = 1,
      R00 = 0.04,
      T0 = 0.75,
      R0 = 0.04,
      T0p = 0,
      R0p = 0.04,
      V0 = 10,
      D0 = 0.1 * (10 - 1),
      estim_start = as.Date("2009-01-01"),
      QQ = c(0.16, 0.84),
      fix_lam1 = TRUE,
      fix_psi_grp = FALSE,
      fix_sig2f = FALSE,
      V0f = 20,
      D0f = 3 * (20 - 1),
      progressEvery = 50,
      smoother = "precision",
      sig2_tau = 0.0001
    )
  ),
  tar_target(estVargroups, c(1, 2)),
  tar_target(
    regionOverrides,
    # list(
    #   AR = list(V0f = 20, D0f = 0.1 (20 - 1)),
    #   TO = list(V0f = 20, D0f = 0.1 (20 - 1))
    # )
    NULL
  ),
  tar_target(estFolder, {
    d <- paste0("Results/EST_", runTag, "_NOGDP")
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    invisible(file.remove(list.files(
      d,
      pattern = "^\\.status_",
      full.names = TRUE,
      all.files = TRUE
    )))
    d
  }),
  tar_target(estFolderGdp, {
    d <- paste0("Results/EST_", runTag, "_GDP")
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    invisible(file.remove(list.files(
      d,
      pattern = "^\\.status_",
      full.names = TRUE,
      all.files = TRUE
    )))
    d
  }),
  tar_target(
    evalSettings,
    list(
      horizons = c(-6, -5, -4, -3, -2, -1, 0, 3, 6, 12),
      minGrpObs = 10,
      grpReleaseLag = 15,
      gdpReleaseLag = 3,
      folder = paste0("Results/EVAL_", runTag, "_NOGDP")
    )
  ),
  tar_target(evalSettingsGdp, {
    s <- evalSettings
    s$folder <- paste0("Results/EVAL_", runTag, "_GDP")
    s
  }),
  tar_target(
    evalMcmcSettings,
    list(M0 = 5000, M = 5000)
  ),
  tar_target(
    cutoffGrid,
    defineCutoffGrid(
      evalSettings$horizons,
      evalSettings$minGrpObs,
      evalSettings$grpReleaseLag,
      cfg$estim_start,
      cfg$sample_end
    )
  ),
  tar_target(evalCutoffIds, cutoffGrid$cutoffs$cutoff_id),

  # ---- File Tracking ---------------------------------------------------------
  tar_target(inputFile, "Data/РегДанные.xlsx"),

  # ---- Metadata --------------------------------------------------------------
  tar_target(codes, readCodes(inputFile)),
  tar_target(codesTable, readCodesTable(inputFile)),
  tar_target(metadata, readMetadata(inputFile)),

  # ---- Monthly: Inflation & Investment Prices (need SA) ----------------------
  tar_target(inflation, readInflation(inputFile)),
  tar_target(
    inflationSA,
    seasAdj(inflation) %>%
      filter(date >= cfg$sample_start & date <= cfg$sample_end)
  ),
  tar_target(invprice, readInvprice(inputFile)),
  tar_target(
    invpriceSA,
    seasAdj(invprice) %>%
      filter(date >= cfg$sample_start & date <= cfg$sample_end)
  ),

  # ---- Monthly: Simple YoY / raw indicators ----------------------------------
  tar_target(ivbom, readIvbom(inputFile, codes, cfg)),
  tar_target(industry, readIndustry(inputFile, cfg)),
  tar_target(wholesale, readWholesale(inputFile, cfg)),
  tar_target(retail, readRetail(inputFile, cfg)),
  tar_target(construction, readConstruction(inputFile, cfg)),
  tar_target(services, readServices(inputFile, cfg)),
  tar_target(rest, readRest(inputFile, cfg)),

  # ---- Monthly: Wages (deflated, YoY) ----------------------------------------
  tar_target(rwage, {
    wage <- readWage(inputFile, codes)
    toReal(wage, inflation) %>%
      group_by(region, OKATO_id) %>%
      mutate(value = value / lag(value, 12) * 100) %>%
      filter(date >= cfg$sample_start & date <= cfg$sample_end)
  }),

  # ---- Monthly: Labour Market ------------------------------------------------
  tar_target(labourRaw, readLabour(inputFile, codes, cfg)),
  tar_target(
    unempSA,
    seasAdj(labourRaw$unemp) %>%
      filter(date >= cfg$sample_start & date <= cfg$sample_end)
  ),
  tar_target(
    vacancylfSA,
    seasAdj(labourRaw$vacancylf) %>%
      filter(date >= cfg$sample_start & date <= cfg$sample_end)
  ),
  tar_target(
    vacancyuSA,
    seasAdj(labourRaw$vacancyu) %>%
      filter(date >= cfg$sample_start & date <= cfg$sample_end)
  ),
  tar_target(vacancyyoy, labourRaw$vacancyyoy),

  # ---- Monthly: Corporate Profits (quarterly freq, SA, deflated) -------------
  tar_target(profits, readProfits(inputFile, codes)),
  tar_target(profitsSA, seasAdj(profits, freq = 4, prelog = FALSE)),
  tar_target(
    rprofitsSA,
    toReal(profitsSA, inflationSA) %>%
      filter(date >= cfg$sample_start & date <= cfg$sample_end)
  ),

  # ---- Quarterly Variables ---------------------------------------------------
  tar_target(ivboq, readIvboq(inputFile, codes, cfg)),
  tar_target(rincome, readRincome(inputFile, codes, cfg)),
  tar_target(rinvestment, {
    inv <- readInvestment(inputFile, codes)
    toReal(inv, invpriceSA, freq = 4) %>%
      group_by(region, OKATO_id) %>%
      mutate(value = value / lag(value, 4) * 100) %>%
      filter(date >= cfg$sample_start & date <= cfg$sample_end)
  }),
  tar_target(transport, readTransport(inputFile, codes, cfg)),
  tar_target(transportation, readTransportation(inputFile, codes, cfg)),
  tar_target(agri, readAgri(inputFile, codes, cfg)),

  # ---- Sentiment Indicators (quarterly, SA) ----------------------------------
  tar_target(sentconstrSA, {
    raw <- readSentConstr(inputFile, codes)
    filled <- fillNA(raw, freq = 4)
    seasAdj(filled, freq = 4, prelog = FALSE) %>%
      filter(date >= cfg$sample_start & date <= cfg$sample_end)
  }),
  tar_target(sentwholesaleSA, {
    raw <- readSentWholesale(inputFile, codes)
    filled <- fillNA(raw, freq = 4)
    seasAdj(filled, freq = 4, prelog = FALSE) %>%
      filter(date >= cfg$sample_start & date <= cfg$sample_end)
  }),
  tar_target(sentretailSA, {
    raw <- readSentRetail(inputFile, codes)
    filled <- fillNA(raw, freq = 4)
    seasAdj(filled, freq = 4, prelog = FALSE) %>%
      filter(date >= cfg$sample_start & date <= cfg$sample_end)
  }),
  tar_target(sentservicesSA, {
    raw <- readSentServices(inputFile, codes)
    seasAdj(raw, freq = 4, prelog = FALSE) %>%
      filter(date >= cfg$sample_start & date <= cfg$sample_end)
  }),

  # ---- Annual ----------------------------------------------------------------
  tar_target(grp, readGrp(inputFile, codes)),
  tar_target(ivboYtd, readIvboYtd(inputFile, codes)),
  tar_target(gdpNat, readGdpNat(inputFile)),

  # ---- Assembly: join all variables ------------------------------------------
  tar_target(
    econ,
    joinAllVariables(
      list(
        inflation = inflationSA,
        invprice = invpriceSA,
        ivbom = ivbom,
        industry = industry,
        wholesale = wholesale,
        retail = retail,
        construction = construction,
        services = services,
        rest = rest,
        rwage = rwage,
        unemp = unempSA,
        vacancylf = vacancylfSA,
        vacancyu = vacancyuSA,
        vacancyyoy = vacancyyoy,
        rprofits = rprofitsSA,
        ivboq = ivboq,
        rincome = rincome,
        rinvestment = rinvestment,
        transport = transport,
        transportation = transportation,
        agri = agri,
        sentconstr = sentconstrSA,
        sentwholesale = sentwholesaleSA,
        sentretail = sentretailSA,
        sentservices = sentservicesSA,
        grp = grp
      ),
      metadata,
      codes
    )
  ),

  # ---- Summary Statistics ----------------------------------------------------
  tar_target(
    summaryStats,
    buildSummaryStats(econ, metadata, codes, codesTable)
  ),
  tar_target(
    summaryStatsFile,
    {
      saveRDS(summaryStats, "Data/summaryStats.RDS")
      "Data/summaryStats.RDS"
    },
    format = "file"
  ),

  # ---- State-Space Panel (monthly.xlsx) --------------------------------------
  tar_target(
    econSS,
    buildEconSS(econ, metadata, cfg$estim_start, cfg$sample_end)
  ),
  tar_target(
    monthlyXlsx,
    {
      writexl::write_xlsx(
        list(econSS = econSS, transformations = metadata$meta),
        "Data/monthly.xlsx"
      )
      "Data/monthly.xlsx"
    },
    format = "file"
  ),

  # ---- Weight Matrices -------------------------------------------------------
  tar_target(
    weightGrp,
    {
      w <- buildWeightMatrix(inputFile, codes, cfg$estim_start, cfg$sample_end)
      saveRDS(w, "Data/weight_grp.RDS")
      "Data/weight_grp.RDS"
    },
    format = "file"
  ),
  tar_target(
    weightEqual,
    {
      n_regs <- nrow(codes)
      n_months <- length(seq.Date(
        cfg$estim_start,
        cfg$sample_end,
        by = "month"
      ))
      w <- matrix(1 / n_regs, nrow = n_regs, ncol = n_months)
      saveRDS(w, "Data/weight_equal.RDS")
      "Data/weight_equal.RDS"
    },
    format = "file"
  ),

  # ---- Source Data (read the panel produced by data pipeline) -----------------
  tar_target(sourceDataRaw, read_excel(monthlyXlsx)),
  tar_target(stdParams, computeStdParams(sourceDataRaw, noStdVars = "^grp\\.")),
  tar_target(sourceData, applyStandardization(sourceDataRaw, stdParams)),
  tar_target(sourceDataRawGdp, addGdpToSourceData(sourceDataRaw, gdpNat)),
  tar_target(
    stdParamsGdp,
    computeStdParams(sourceDataRawGdp, noStdVars = "^grp\\.")
  ),
  tar_target(
    sourceDataGdp,
    applyStandardization(sourceDataRawGdp, stdParamsGdp)
  ),

  # ---- Pre-compile Rcpp (once, before parallel workers) ----------------------
  tar_target(cppCompiled, {
    dir.create(".rcpp_cache", showWarnings = FALSE)
    Rcpp::sourceCpp("R/functions.cpp", cacheDir = ".rcpp_cache")
    Rcpp::sourceCpp("R/precision_sampler.cpp", cacheDir = ".rcpp_cache")
    TRUE
  }),

  # ---- Region List -----------------------------------------------------------
  tar_target(regionList, {
    vn <- colnames(sourceData)
    regs <- unique(sapply(strsplit(vn[4:length(vn)], "\\."), function(x) x[2]))
    # regs <- c('MW', 'MO')
  }),

  # ---- Per-Region Estimation (dynamic branching) -----------------------------
  tar_target(
    estResultNoGdp,
    {
      force(cppCompiled)
      outpath <- estimRegionTarget(
        regionList,
        sourceData,
        estSettings,
        estVargroups,
        estFolder,
        grp,
        stdParams,
        regionOverrides
      )
      read_excel(outpath, sheet = "factors") %>%
        select(date, grpBase, grpMoM, grpYoY) %>%
        mutate(region = regionList, .before = 1)
    },
    pattern = map(regionList)
  ),
  tar_target(
    estResultGdp,
    {
      force(cppCompiled)
      outpath <- estimRegionTarget(
        regionList,
        sourceDataGdp,
        estSettings,
        estVargroups,
        estFolderGdp,
        grp,
        stdParamsGdp,
        regionOverrides
      )
      read_excel(outpath, sheet = "factors") %>%
        select(date, grpBase, grpMoM, grpYoY) %>%
        mutate(region = regionList, .before = 1)
    },
    pattern = map(regionList)
  ),

  # ---- MCMC Diagnostics -------------------------------------------------------
  tar_target(
    diagNoGdp,
    {
      estResultNoGdp
      computeAllDiagnostics(estFolder)
    }
  ),
  tar_target(
    diagGdp,
    {
      estResultGdp
      computeAllDiagnostics(estFolderGdp)
    }
  ),

  # ---- Evaluation: Per-Region-Cutoff Estimation ------------------------------
  tar_target(
    evalResultNoGdp,
    {
      force(cppCompiled)
      ct <- cutoffGrid$cutoffs %>% filter(cutoff_id == evalCutoffIds)
      truncRaw <- truncateSourceData(
        sourceDataRaw,
        ct$cutoff_date,
        ct$grp_available_through
      )
      truncStdParams <- computeStdParams(truncRaw, noStdVars = "^grp\\.")
      truncData <- applyStandardization(truncRaw, truncStdParams)
      folder <- file.path(
        evalSettings$folder,
        sprintf("CUT_%s", format(ct$cutoff_date, "%Y_%m"))
      )
      dir.create(folder, showWarnings = FALSE, recursive = TRUE)
      evalEstSettings <- estSettings
      evalEstSettings$M0 <- evalMcmcSettings$M0
      evalEstSettings$M <- evalMcmcSettings$M
      outpath <- estimRegionTarget(
        regionList,
        truncData,
        evalEstSettings,
        estVargroups,
        folder,
        grp,
        truncStdParams
      )
      fitted <- read_excel(outpath, sheet = "fitted")
      grpCol <- names(fitted)[grepl("^grp\\.", names(fitted))][1]
      if (is.na(grpCol)) {
        return(tibble(
          region = character(),
          target_year = integer(),
          horizon = integer(),
          cutoff_id = integer(),
          model_forecast = double(),
          draws = list()
        ))
      }
      fittedDates <- as.Date(fitted$date)
      grpDraws <- read_excel(outpath, sheet = "grpmm")
      evalTargets <- cutoffGrid$evalMap %>% filter(cutoff_id == evalCutoffIds)
      purrr::map_dfr(seq_len(nrow(evalTargets)), function(j) {
        targetDec <- as.Date(sprintf("%d-12-01", evalTargets$target_year[j]))
        idx <- match(targetDec, fittedDates)
        if (is.na(idx)) {
          return(NULL)
        }
        tibble(
          region = regionList,
          target_year = evalTargets$target_year[j],
          horizon = evalTargets$horizon[j],
          cutoff_id = evalCutoffIds,
          model_forecast = as.numeric(fitted[idx, grpCol]),
          draws = list(as.numeric(grpDraws[[as.character(targetDec)]]))
        )
      })
    },
    pattern = cross(evalCutoffIds, regionList),
    garbage_collection = TRUE
  ),
  tar_target(
    evalResultGdp,
    {
      force(cppCompiled)
      ct <- cutoffGrid$cutoffs %>% filter(cutoff_id == evalCutoffIds)
      truncRaw <- truncateSourceData(
        sourceDataRawGdp,
        ct$cutoff_date,
        ct$grp_available_through
      )
      truncRaw <- truncateGdpInSourceData(
        truncRaw,
        ct$cutoff_date,
        evalSettingsGdp$gdpReleaseLag
      )
      truncStdParams <- computeStdParams(truncRaw, noStdVars = "^grp\\.")
      truncData <- applyStandardization(truncRaw, truncStdParams)
      folder <- file.path(
        evalSettingsGdp$folder,
        sprintf("CUT_%s", format(ct$cutoff_date, "%Y_%m"))
      )
      dir.create(folder, showWarnings = FALSE, recursive = TRUE)
      evalEstSettings <- estSettings
      evalEstSettings$M0 <- evalMcmcSettings$M0
      evalEstSettings$M <- evalMcmcSettings$M
      outpath <- estimRegionTarget(
        regionList,
        truncData,
        evalEstSettings,
        estVargroups,
        folder,
        grp,
        truncStdParams
      )
      fitted <- read_excel(outpath, sheet = "fitted")
      grpCol <- names(fitted)[grepl("^grp\\.", names(fitted))][1]
      if (is.na(grpCol)) {
        return(tibble(
          region = character(),
          target_year = integer(),
          horizon = integer(),
          cutoff_id = integer(),
          model_forecast = double(),
          draws = list()
        ))
      }
      fittedDates <- as.Date(fitted$date)
      grpDraws <- read_excel(outpath, sheet = "grpmm")
      evalTargets <- cutoffGrid$evalMap %>% filter(cutoff_id == evalCutoffIds)
      purrr::map_dfr(seq_len(nrow(evalTargets)), function(j) {
        targetDec <- as.Date(sprintf("%d-12-01", evalTargets$target_year[j]))
        idx <- match(targetDec, fittedDates)
        if (is.na(idx)) {
          return(NULL)
        }
        tibble(
          region = regionList,
          target_year = evalTargets$target_year[j],
          horizon = evalTargets$horizon[j],
          cutoff_id = evalCutoffIds,
          model_forecast = as.numeric(fitted[idx, grpCol]),
          draws = list(as.numeric(grpDraws[[as.character(targetDec)]]))
        )
      })
    },
    pattern = cross(evalCutoffIds, regionList),
    garbage_collection = TRUE
  ),

  # ---- Evaluation: Score Forecasts -------------------------------------------
  tar_target(
    evalResultsNoGdp,
    scoreFromForecasts(
      bind_rows(evalResultNoGdp),
      cutoffGrid,
      grp,
      ivboYtd,
      weightGrp,
      codes,
      cfg
    )
  ),
  tar_target(
    evalResultsGdp,
    scoreFromForecasts(
      bind_rows(evalResultGdp),
      cutoffGrid,
      grp,
      ivboYtd,
      weightGrp,
      codes,
      cfg
    )
  ),

  # ---- Consolidated Factor Output ---------------------------------------------
  tar_target(
    factorsDfNoGdp,
    bind_rows(estResultNoGdp)
  ),
  tar_target(
    factorsFileNoGdp,
    {
      saveRDS(factorsDfNoGdp, file.path(estFolder, "factors.RDS"))
      file.path(estFolder, "factors.RDS")
    },
    format = "file"
  ),
  tar_target(
    factorsDfGdp,
    bind_rows(estResultGdp)
  ),
  tar_target(
    factorsFileGdp,
    {
      saveRDS(factorsDfGdp, file.path(estFolderGdp, "factors.RDS"))
      file.path(estFolderGdp, "factors.RDS")
    },
    format = "file"
  )
)
