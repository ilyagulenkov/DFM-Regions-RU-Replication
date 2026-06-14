# Helper Functions ----
hpdInterval <- function(draws, prob) {
  draws <- sort(draws)
  n <- length(draws)
  window <- floor(prob * n)
  widths <- draws[(window + 1):n] - draws[1:(n - window)]
  i <- which.min(widths)
  list(l = draws[i], u = draws[i + window])
}

safeChol <- function(V) {
  res <- tryCatch(chol(V), error = function(e) NULL)
  if (!is.null(res)) return(res)
  R <- suppressWarnings(chol(V, pivot = TRUE))
  piv <- attr(R, "pivot")
  R[, order(piv)]
}

lag0 <- function(x, p) {
  R <- nrow(x)
  C <- ncol(x)
  x1 <- x[1:(R - p), , drop = FALSE]
  out <- rbind(matrix(0, p, C), x1)
  return(out)
}
delif <- function(y, cond) {
  x <- y[cond == 0]
  return(x)
}

# Data Preparation ----
seasAdj <- function(df, freq = 12, prelog = TRUE, method = 'BV4.1') {
  uniqs <- df %>%
    distinct(region, OKATO_id)
  dfSA <- foreach(
    i = 1:nrow(uniqs),
    .combine = rbind,
    .packages = c('tidyverse', 'readxl', 'seasonal', 'deseats', 'seastests')
  ) %do%
    {
      uniqrow <- uniqs[i, ]
      sub <- df %>%
        filter(OKATO_id == uniqrow$OKATO_id) %>%
        mutate(err = F)
      TT <- nrow(sub)
      sub2 <- sub %>%
        filter(!is.na(value))
      TT2 <- nrow(sub2)
      periodfun <- ifelse(freq == 12, month, quarter)
      if (prelog) {
        unsa <- ts(
          log(sub2$value),
          frequency = freq,
          start = c(year(first(sub2$date)), periodfun(first(sub2$date)))
        )
      } else {
        unsa <- ts(
          sub2$value,
          frequency = freq,
          start = c(year(first(sub2$date)), periodfun(first(sub2$date)))
        )
      }
      if (method == 'BV4.1') {
        sa <- try(BV4.1(unsa) %>% deseasonalize())
      } else {
        sa2 <- try(RJDemetra::x13(unsa))
        sa2 <- sa2$final$series[, 'sa']
        sa <- c(rep(NA, TT - TT2), sa2)
      }
      if (class(sa) == 'try-error') {
        sub$err <- T
        sub$value <- NA
      } else {
        sub$seasdumm <- seasdum(sa, freq = freq)$Pval
        sub$kw <- kw(sa, freq = freq)$Pval
        if (prelog) {
          sub$value <- c(rep(NA, TT - TT2), exp(as.numeric(sa)))
        } else {
          sub$value <- c(rep(NA, TT - TT2), as.numeric(sa))
        }
      }
      return(sub)
    }
  return(dfSA)
}
toReal <- function(df, cpi, baseyear = 2023, freq = 12) {
  cpi <- cpi %>%
    select(region, date, cpi = value) %>%
    group_by(region) %>%
    mutate(
      cpi = cpi /
        mean(cpi[
          date >= as.Date(glue('{baseyear}-01-01')) &
            date < as.Date(glue('{baseyear+1}-01-01'))
        ])
    )
  if (freq == 4) {
    cpi <- cpi %>%
      group_by(region, date = as.Date(yearquarter(date)) + months(2)) %>%
      summarize(cpi = mean(cpi))
  }
  df %>%
    left_join(cpi) %>%
    mutate(value = value / cpi) %>%
    select(-cpi)
}
fillNA <- function(df, freq = 12) {
  uniqs <- df %>%
    distinct(region, OKATO_id)
  dfImp <- foreach(
    i = 1:nrow(uniqs),
    .combine = rbind,
    .packages = c('tidyverse', 'readxl', 'imputeTS')
  ) %do%
    {
      uniqrow <- uniqs[i, ]
      sub <- df %>%
        filter(OKATO_id == uniqrow$OKATO_id) %>%
        mutate(err = F)
      periodfun <- ifelse(freq == 12, month, quarter)
      series <- ts(
        sub$value,
        frequency = freq,
        start = c(year(first(sub$date)), periodfun(first(sub$date)))
      )
      series_imputed <- na_kalman(series)
      sub$value <- as.numeric(series_imputed)
      return(sub)
    }
  return(dfImp)
}
transformData <- function(
  df,
  extra_vars,
  levels = NULL,
  logLevels = NULL,
  logLevels_by100 = NULL,
  diffM = NULL,
  diffQ = NULL,
  diffY = NULL,
  lDiffM = NULL,
  lDiffQ = NULL,
  lDiffY = NULL
) {
  df %>%
    pivot_longer(-c(date, all_of(extra_vars))) %>%
    group_by(across(all_of(c(extra_vars, 'name')))) %>%
    mutate(
      value = suppressWarnings(case_when(
        name %in% levels ~ value,
        name %in% logLevels ~ log(value) * 100,
        name %in% logLevels_by100 ~ log(value / 100) * 100,
        name %in% diffM ~ difference(value, 1),
        name %in% diffQ ~ difference(value, 3),
        name %in% diffY ~ difference(value, 12),
        name %in% lDiffM ~ difference(log(value), 1) * 100,
        name %in% lDiffQ ~ difference(log(value), 3) * 100,
        name %in% lDiffY ~ difference(log(value), 12) * 100
      ))
    )
}
tabForm <- function(x, cap = NULL, dig = 2, na_str = "\u2013") {
  x %>%
    # autofit() %>%
    set_table_properties(layout = "autofit", width = 1) %>%
    flextable::font(fontname = 'Times New Roman', part = 'all') %>%
    padding(padding = 0, part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::align(align = 'justify', part = 'footer') %>%
    colformat_double(big.mark = "", digits = dig, na_str = na_str) %>%
    flextable::set_caption(
      cap,
      align_with_table = F,
      word_stylename = "Table Caption",
      fp_p = officer::fp_par(text.align = "justify")
    )
}

# Metadata Readers ----
readCodes <- function(path) {
  read_excel(path, sheet = 'КОДЫ') %>%
    select(region = name_official, Abbrev = abbrev, OKATO_id) %>%
    arrange(Abbrev) %>%
    filter(!is.na(Abbrev))
}
readCodesTable <- function(path) {
  read_excel(path, sheet = 'КОДЫ') %>%
    select(OKATO_id, region = name_rus)
}
readMetadata <- function(path) {
  meta <- read_xlsx(path, sheet = 'КодыПеременные')
  renameVector <- deframe(meta[, 1:2])
  renameVectorEng <- deframe(meta[, c(4, 2)])
  renameVectorEngShort <- deframe(meta[, c(5, 2)])
  renameVectorRu <- deframe(meta[, c(3, 2)])
  meta <- meta %>%
    column_to_rownames('name')
  idx_nsa <- which(meta$nsa == 1)
  list(
    meta = meta,
    renameVector = renameVector,
    renameVectorEng = renameVectorEng,
    renameVectorEngShort = renameVectorEngShort,
    renameVectorRu = renameVectorRu,
    idx_nsa = idx_nsa
  )
}

# Monthly Data Readers ----
readInflation <- function(path) {
  read_excel(
    path,
    sheet = 'ИПЦ (2008-2025)',
    range = cell_limits(c(110, 2), c(NA, NA))
  ) %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    filter(!is.na(OKATO_id)) %>%
    select(region, OKATO_id, date, value) %>%
    mutate(date = as.Date(as.numeric(date), origin = '1899-12-30'))
}
readInvprice <- function(path) {
  read_excel(
    path,
    sheet = 'ИнвестЦены (2005-2025)',
    range = cell_limits(c(103, 4), c(NA, NA))
  ) %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    filter(!is.na(OKATO_id)) %>%
    select(region, OKATO_id, date, value) %>%
    mutate(date = as.Date(as.numeric(date), origin = '1899-12-30'))
}
readIvbom <- function(path, codes, cfg) {
  read_excel(
    path,
    sheet = 'РегИВБО (2022-2025)',
    range = cell_limits(c(1, 2), c(NA, NA))
  ) %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(date = as.Date(as.numeric(date), origin = '1899-12-30')) %>%
    left_join(codes, by = c("region", "OKATO_id")) %>%
    filter(!is.na(OKATO_id)) %>%
    filter(date >= cfg$sample_start & date <= cfg$sample_end)
}
readIndustry <- function(path, cfg) {
  read_excel(path, sheet = 'Пром г-г (2009-2024)') %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(date = as.Date(as.numeric(date), origin = '1899-12-30')) %>%
    filter(!is.na(OKATO_id)) %>%
    filter(date >= cfg$sample_start & date <= cfg$sample_end)
}
readWholesale <- function(path, cfg) {
  read_excel(
    path,
    sheet = 'Опт г-г (2009-2025)',
    range = cell_limits(c(1, 2), c(NA, NA))
  ) %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(date = as.Date(as.numeric(date), origin = '1899-12-30')) %>%
    filter(!is.na(OKATO_id)) %>%
    filter(date >= cfg$sample_start & date <= cfg$sample_end)
}
readRetail <- function(path, cfg) {
  read_excel(path, sheet = 'Розница г-г (2009-2024)') %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(date = as.Date(as.numeric(date), origin = '1899-12-30')) %>%
    filter(!is.na(OKATO_id)) %>%
    filter(date >= cfg$sample_start & date <= cfg$sample_end)
}
readConstruction <- function(path, cfg) {
  read_excel(path, sheet = 'Стройка г-г (2009-2024)') %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(date = as.Date(as.numeric(date), origin = '1899-12-30')) %>%
    filter(!is.na(OKATO_id)) %>%
    filter(date >= cfg$sample_start & date <= cfg$sample_end) %>%
    filter(!region == 'Республика Ингушетия')
}
readServices <- function(path, cfg) {
  read_excel(path, sheet = 'Услуги г-г (2009-2024)') %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(date = as.Date(as.numeric(date), origin = '1899-12-30')) %>%
    filter(!is.na(OKATO_id)) %>%
    filter(date >= cfg$sample_start & date <= cfg$sample_end)
}
readRest <- function(path, cfg) {
  read_excel(
    path,
    sheet = 'Общепит г-г (2016-2025)',
    range = cell_limits(c(1, 2), c(NA, NA))
  ) %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(date = as.Date(as.numeric(date), origin = '1899-12-30')) %>%
    filter(!is.na(OKATO_id)) %>%
    filter(date >= cfg$sample_start & date <= cfg$sample_end)
}
readWage <- function(path, codes) {
  read_excel(path, sheet = 'Зарплата (2006-2025)') %>%
    pivot_longer(-c(region), names_to = 'date') %>%
    mutate(date = as.Date(as.numeric(date), origin = '1899-12-30')) %>%
    left_join(codes, by = "region") %>%
    filter(!is.na(OKATO_id))
}
readLabour <- function(path, codes, cfg) {
  ldemand <- read_excel(path, sheet = 'СпросТруд (2005-2025)') %>%
    pivot_longer(-c(region), names_to = 'date') %>%
    mutate(date = as.Date(as.numeric(date), origin = '1899-12-30')) %>%
    rename(c('ldemand' = 'value'))
  labour <- read_excel(path, sheet = 'РынокТруда (2009-2025)') %>%
    pivot_longer(-c(region, name), names_to = 'date') %>%
    mutate(date = as.Date(as.numeric(date), origin = '1899-12-30')) %>%
    pivot_wider(names_from = 'name', values_from = 'value') %>%
    full_join(ldemand, by = c("region", "date")) %>%
    arrange(region, date) %>%
    left_join(codes, by = "region") %>%
    group_by(region, OKATO_id) %>%
    mutate(
      u = U / LF * 100,
      ld = ldemand / lag(ldemand, 12) * 100,
      v = ldemand / LF / 10^3 * 100,
      vu = ldemand / 10^3 / U
    ) %>%
    filter(!is.na(OKATO_id))
  unemp <- labour %>%
    select(region, OKATO_id, date, value = u)
  vacancylf <- labour %>%
    select(region, OKATO_id, date, value = v)
  vacancyu <- labour %>%
    select(region, OKATO_id, date, value = vu)
  vacancyyoy <- labour %>%
    select(region, OKATO_id, date, value = ld) %>%
    filter(date >= cfg$sample_start & date <= cfg$sample_end)
  list(
    unemp = unemp,
    vacancylf = vacancylf,
    vacancyu = vacancyu,
    vacancyyoy = vacancyyoy
  )
}
readProfits <- function(path, codes) {
  read_excel(path, sheet = 'Прибыли (2009-2025)') %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(date = as.Date(as.numeric(date), origin = '1899-12-30')) %>%
    left_join(codes, by = c("region", "OKATO_id")) %>%
    filter(!is.na(OKATO_id))
}

# Quarterly Data Readers ----
readIvboq <- function(path, codes, cfg) {
  read_excel(
    path,
    sheet = 'РегИВБО (2018-2025)',
    range = cell_limits(c(1, 2), c(NA, NA))
  ) %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(
      date = as.Date(as.numeric(date), origin = '1899-12-30') + months(2)
    ) %>%
    left_join(codes, by = c("region", "OKATO_id")) %>%
    filter(!is.na(OKATO_id)) %>%
    filter(date >= cfg$sample_start & date <= cfg$sample_end)
}
readRincome <- function(path, codes, cfg) {
  read_excel(
    path,
    sheet = 'ДенДох г-г (2009-2025)',
    range = cell_limits(c(1, 2), c(NA, NA))
  ) %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(
      date = as.Date(as.numeric(date), origin = '1899-12-30') + months(2)
    ) %>%
    left_join(codes, by = c("region", "OKATO_id")) %>%
    filter(!is.na(OKATO_id)) %>%
    filter(date >= cfg$sample_start & date <= cfg$sample_end)
}
readInvestment <- function(path, codes) {
  read_excel(path, sheet = 'Инвестиции (2009-2025)') %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(
      date = as.Date(as.numeric(date), origin = '1899-12-30') + months(2)
    ) %>%
    left_join(codes, by = c("region", "OKATO_id")) %>%
    filter(!is.na(OKATO_id))
}
readTransport <- function(path, codes, cfg) {
  read_excel(path, sheet = 'ПеревозкаГрузов (2008-2025)') %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(
      date = as.Date(as.numeric(date), origin = '1899-12-30') + months(2)
    ) %>%
    left_join(codes, by = c("region", "OKATO_id")) %>%
    filter(!is.na(OKATO_id)) %>%
    group_by(region, OKATO_id) %>%
    mutate(value = value / lag(value, 4) * 100) %>%
    filter(date >= cfg$sample_start & date <= cfg$sample_end)
}
readTransportation <- function(path, codes, cfg) {
  read_excel(path, sheet = 'Грузооборот (2008-2025)') %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(
      date = as.Date(as.numeric(date), origin = '1899-12-30') + months(2)
    ) %>%
    left_join(codes, by = c("region", "OKATO_id")) %>%
    filter(!is.na(OKATO_id)) %>%
    group_by(region, OKATO_id) %>%
    mutate(value = value / lag(value, 4) * 100) %>%
    filter(date >= cfg$sample_start & date <= cfg$sample_end)
}
readAgri <- function(path, codes, cfg) {
  read_excel(
    path,
    sheet = 'СельХоз г-г (2011-2025)',
    range = cell_limits(c(1, 2), c(NA, NA))
  ) %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(
      date = as.Date(as.numeric(date), origin = '1899-12-30') + months(2)
    ) %>%
    left_join(codes, by = c("region", "OKATO_id")) %>%
    filter(!is.na(OKATO_id)) %>%
    filter(!OKATO_id == '40000000000') %>%
    filter(!OKATO_id == '45000000000') %>%
    filter(date >= cfg$sample_start & date <= cfg$sample_end)
}

# Sentiment Indicator Readers ----
readSentConstr <- function(path, codes) {
  read_excel(
    path,
    sheet = 'ИПУСтроит (2005-2025)',
    range = cell_limits(c(1, 2), c(NA, NA))
  ) %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(
      date = as.Date(as.numeric(date), origin = '1899-12-30') + months(2)
    ) %>%
    left_join(codes, by = c("region", "OKATO_id")) %>%
    filter(!is.na(OKATO_id)) %>%
    filter(date >= as.Date('2010-04-01'))
}
readSentWholesale <- function(path, codes) {
  raw <- read_excel(
    path,
    sheet = 'ИПУОпт (2013-2025)',
    range = cell_limits(c(312, 3), c(NA, NA))
  ) %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(
      date = as.Date(as.numeric(date), origin = '1899-12-30') + months(2)
    ) %>%
    left_join(codes, by = c("region", "OKATO_id")) %>%
    filter(!is.na(OKATO_id))
  excl <- raw %>%
    group_by(region, OKATO_id) %>%
    summarize(n = sum(is.na(value)), .groups = "drop") %>%
    filter(n > 5) %>%
    pull(OKATO_id)
  raw %>% filter(!OKATO_id %in% excl)
}
readSentRetail <- function(path, codes) {
  read_excel(
    path,
    sheet = 'ИПУРозн (2010-2025)',
    range = cell_limits(c(1, 2), c(NA, NA))
  ) %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(
      date = as.Date(as.numeric(date), origin = '1899-12-30') + months(2)
    ) %>%
    left_join(codes, by = c("region", "OKATO_id")) %>%
    filter(!is.na(OKATO_id)) %>%
    filter(date >= as.Date('2011-01-01'))
}
readSentServices <- function(path, codes) {
  read_excel(
    path,
    sheet = 'ИПУУслуги (2012-2024)',
    range = cell_limits(c(1, 2), c(NA, NA))
  ) %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(
      date = as.Date(as.numeric(date), origin = '1899-12-30') + months(2)
    ) %>%
    left_join(codes, by = c("region", "OKATO_id")) %>%
    filter(!is.na(OKATO_id))
}

# Annual Data Readers ----
readGrp <- function(path, codes) {
  read_excel(path, sheet = 'ВРППостЦены (2006-2023)') %>%
    pivot_longer(-c(region, OKATO_id), names_to = 'date') %>%
    mutate(
      date = as.Date(as.numeric(date), origin = '1899-12-30') + months(11)
    ) %>%
    filter(date >= as.Date('2008-01-01')) %>%
    left_join(codes, by = c("region", "OKATO_id"))
}

readIvboYtd <- function(path, codes) {
  read_excel(path, sheet = "РегИВБО YTD (2018-2025)") %>%
    select(-1) %>%
    pivot_longer(-c(region, OKATO_id), names_to = "date") %>%
    filter(!is.na(value), value > 0) %>%
    mutate(
      date = as.Date(as.numeric(date), origin = "1899-12-30"),
      year = year(date),
      ivboYtdGrowth = log(value / 100) * 100
    ) %>%
    left_join(codes, by = c("region", "OKATO_id")) %>%
    filter(!is.na(Abbrev)) %>%
    select(Abbrev, OKATO_id, year, ivboYtdGrowth)
}

readGdpNat <- function(inputFile) {
  read_excel(inputFile, sheet = "ФЕДЕРАЛЬНЫЕ") %>%
    select(date, gdp) %>%
    mutate(date = as.Date(date))
}

# Assembly Functions ----
joinAllVariables <- function(varList, metadata, codes) {
  meta <- metadata$meta
  renameVector <- metadata$renameVector
  idx_nsa <- metadata$idx_nsa
  extra_vars <- c('region', 'OKATO_id', 'date')
  join_by <- c('region', 'OKATO_id', 'date')
  varNames <- renameVector
  df <- NULL
  for (i in seq_along(varNames)) {
    displayName <- names(varNames)[i]
    internalName <- varNames[i]
    if (is.null(varList[[internalName]])) {
      next
    }
    piece <- varList[[internalName]] %>%
      ungroup() %>%
      select(all_of(extra_vars), !!internalName := value) %>%
      distinct()
    if (is.null(df)) {
      df <- piece
    } else {
      df <- df %>% full_join(piece, by = join_by)
    }
  }
  df <- df %>%
    arrange(across(all_of(join_by))) %>%
    left_join(codes, by = c("region", "OKATO_id"))
  df
}
buildSummaryStats <- function(econ, metadata, codes, codesTable) {
  meta <- metadata$meta
  renameVector <- metadata$renameVector
  renameVectorRu <- metadata$renameVectorRu
  renameVectorEng <- metadata$renameVectorEng
  renameVectorEngShort <- metadata$renameVectorEngShort
  summaryStats1 <- transformData(
    econ,
    extra_vars = c('Abbrev', 'region', 'OKATO_id'),
    levels = meta$variable[meta$transf_summary == 0],
    logLevels = meta$variable[meta$transf_summary == 2],
    logLevels_by100 = meta$variable[meta$transf_summary == 3],
    diffM = meta$variable[meta$transf_summary == 1],
    diffQ = meta$variable[meta$transf_summary == 4],
    diffY = meta$variable[meta$transf_summary == 12],
    lDiffM = meta$variable[meta$transf_summary == 10],
    lDiffQ = meta$variable[meta$transf_summary == 40],
    lDiffY = meta$variable[meta$transf_summary == 120]
  ) %>%
    summarize(value = mean(value, na.rm = T)) %>%
    pivot_wider(names_from = 'name', values_from = 'value') %>%
    ungroup() %>%
    select(-region) %>%
    left_join(codesTable, by = "OKATO_id") %>%
    select(-c(OKATO_id, Abbrev)) %>%
    relocate(region, 1) %>%
    arrange(region)
  summaryStats2 <- summaryStats1 %>%
    pivot_longer(-region) %>%
    group_by(name) %>%
    summarize(
      `Минимум` = min(value, na.rm = T),
      `25-й перцентиль` = quantile(value, 0.25, na.rm = T),
      `Среднее` = mean(value, na.rm = T),
      `Медиана` = median(value, na.rm = T),
      `75-й перцентиль` = quantile(value, 0.75, na.rm = T),
      `Максимум` = max(value, na.rm = T)
    ) %>%
    rename('Индикатор' = 'name')
  corstats <- transformData(
    econ,
    extra_vars = c('Abbrev', 'region', 'OKATO_id'),
    levels = meta$variable[meta$transf_model == 0],
    logLevels = meta$variable[meta$transf_model == 2],
    logLevels_by100 = meta$variable[meta$transf_model == 3],
    diffM = meta$variable[meta$transf_model == 1],
    diffQ = meta$variable[meta$transf_model == 4],
    diffY = meta$variable[meta$transf_model == 12],
    lDiffM = meta$variable[meta$transf_model == 10],
    lDiffQ = meta$variable[meta$transf_model == 40],
    lDiffY = meta$variable[meta$transf_model == 120]
  ) %>%
    filter(!is.na(value)) %>%
    group_by(region, Abbrev, name, year(date)) %>%
    summarize(value = mean(value), .groups = "drop") %>%
    ungroup()
  corstats <- corstats %>%
    left_join(
      corstats %>%
        filter(name == "grp") %>%
        select(region, Abbrev, `year(date)`, grp = `value`),
      by = c("region", "Abbrev", "year(date)")
    ) %>%
    filter(!is.na(grp)) %>%
    group_by(name, region) %>%
    summarize(cor = cor(value, grp), .groups = "drop")
  list(
    summaryStats1 = summaryStats1,
    summaryStats2 = summaryStats2,
    corstats = corstats,
    metadata = meta,
    renameVector = renameVector,
    renameVectorRu = renameVectorRu,
    renameVectorEng = renameVectorEng,
    renameVectorEngShort = renameVectorEngShort
  )
}
computeStdParams <- function(sourceDataRaw, noStdVars = NULL) {
  datRows <- 4:nrow(sourceDataRaw)
  varCols <- setdiff(colnames(sourceDataRaw), c("date", "year", "quarter"))
  skipMask <- if (!is.null(noStdVars)) {
    grepl(noStdVars, varCols)
  } else {
    rep(FALSE, length(varCols))
  }
  tibble(
    var = varCols,
    mean = sapply(seq_along(varCols), function(i) {
      if (skipMask[i]) {
        return(0)
      }
      mean(as.numeric(sourceDataRaw[[varCols[i]]][datRows]), na.rm = TRUE)
    }),
    sd = sapply(seq_along(varCols), function(i) {
      if (skipMask[i]) {
        return(1)
      }
      sd(as.numeric(sourceDataRaw[[varCols[i]]][datRows]), na.rm = TRUE)
    })
  )
}
applyStandardization <- function(sourceDataRaw, stdParams) {
  dat <- sourceDataRaw
  datRows <- 4:nrow(dat)
  for (i in seq_len(nrow(stdParams))) {
    v <- stdParams$var[i]
    vals <- as.numeric(dat[[v]][datRows])
    dat[[v]][datRows] <- (vals - stdParams$mean[i]) / stdParams$sd[i]
  }
  dat
}
addGdpToSourceData <- function(sourceDataRaw, gdpNat) {
  dataRows <- 4:nrow(sourceDataRaw)
  sourceDates <- as.Date(sourceDataRaw[dataRows, 1][[1]])
  gdpVals <- log(gdpNat$gdp / 100)
  gdpAligned <- gdpVals[match(sourceDates, gdpNat$date)]
  gdpColVec <- c(3, 2, 2, gdpAligned)
  varCols <- setdiff(colnames(sourceDataRaw), c("date", "year", "quarter"))
  aggs <- as.numeric(sourceDataRaw[2, varCols])
  lastQ <- max(which(aggs == 2))
  insertAfter <- which(colnames(sourceDataRaw) == varCols[lastQ])
  bind_cols(
    sourceDataRaw[, 1:insertAfter],
    tibble(gdp = gdpColVec),
    sourceDataRaw[, (insertAfter + 1):ncol(sourceDataRaw)]
  )
}
truncateGdpInSourceData <- function(sourceData, cutoffDate, gdpReleaseLag = 3) {
  if (!"gdp" %in% colnames(sourceData)) return(sourceData)
  gdpAvailableThrough <- cutoffDate - months(gdpReleaseLag)
  dataRows <- 4:nrow(sourceData)
  dateDat <- as.Date(sourceData[dataRows, 1][[1]])
  maskRows <- dataRows[dateDat > gdpAvailableThrough]
  if (length(maskRows) > 0) sourceData[maskRows, "gdp"] <- NA
  sourceData
}
buildEconSS <- function(econ, metadata, estimStart, sampleEnd) {
  meta <- metadata$meta %>%
    arrange(-freq)
  econSS <- transformData(
    econ,
    extra_vars = c('Abbrev', 'region', 'OKATO_id'),
    levels = meta$variable[meta$transf_model == 0],
    logLevels = meta$variable[meta$transf_model == 2],
    logLevels_by100 = meta$variable[meta$transf_model == 3],
    diffM = meta$variable[meta$transf_model == 1],
    diffQ = meta$variable[meta$transf_model == 4],
    diffY = meta$variable[meta$transf_model == 12],
    lDiffM = meta$variable[meta$transf_model == 10],
    lDiffQ = meta$variable[meta$transf_model == 40],
    lDiffY = meta$variable[meta$transf_model == 120]
  ) %>%
    mutate(value = ifelse(name == 'unemp', -value, value)) %>%
    arrange(Abbrev, match(name, meta$variable)) %>%
    ungroup() %>%
    select(-c(region, OKATO_id)) %>%
    pivot_wider(
      names_from = c('name', 'Abbrev'),
      names_sep = '.',
      names_vary = 'slowest'
    ) %>%
    filter(date >= estimStart & date <= sampleEnd) %>%
    mutate(
      year = as.numeric(month(date) == 12),
      quarter = as.numeric(month(date) %in% c(3, 6, 9, 12))
    ) %>%
    relocate(c(year, quarter), .after = date)
  econSS <- econSS[,
    -which(apply(econSS, 2, function(x) {
      all(is.na(x))
    }))
  ]
  freqs <- sapply(
    substr(names(econSS), 1, nchar(names(econSS)) - 3),
    function(x) {
      val <- meta$freq[meta$variable == x]
      ifelse(length(val) == 0, NA, val)
    }
  )
  names(freqs) <- names(econSS)
  freqs <- unlist(freqs)
  groups <- sapply(
    substr(names(econSS), 1, nchar(names(econSS)) - 3),
    function(x) {
      val <- meta$group[meta$variable == x]
      ifelse(length(val) == 0, NA, val)
    }
  )
  names(groups) <- names(econSS)
  groups <- unlist(groups)
  aggs <- sapply(
    substr(names(econSS), 1, nchar(names(econSS)) - 3),
    function(x) {
      val <- meta$agg_type[meta$variable == x]
      ifelse(length(val) == 0, NA, val)
    }
  )
  names(aggs) <- names(econSS)
  aggs <- unlist(aggs)
  rbind(freqs, aggs, groups, econSS)
}
buildWeightMatrix <- function(path, codes, estimStart, sampleEnd) {
  abbs <- codes$Abbrev
  regs <- codes$region
  w_raw <- read_excel(
    path,
    sheet = 'ВРПДоля (2005-2023)',
    range = cell_limits(c(3, 2), c(NA, NA))
  ) %>%
    rename(region = name_official) %>%
    left_join(codes %>% select(-OKATO_id), by = "region") %>%
    filter(!is.na(Abbrev)) %>%
    select(-region) %>%
    pivot_longer(-Abbrev, names_to = 'year') %>%
    mutate(
      year = as.integer(format(
        as.Date(as.numeric(year), origin = '1899-12-30'),
        "%Y"
      ))
    )
  dd <- seq.Date(estimStart, sampleEnd, by = 'month')
  date_df <- data.frame(date = dd) %>%
    mutate(year = year(date))
  w_expanded <- expand.grid(Abbrev = abbs, year = unique(date_df$year)) %>%
    left_join(w_raw, by = c("Abbrev", "year")) %>%
    arrange(Abbrev, year) %>%
    group_by(Abbrev) %>%
    fill(value, .direction = "downup") %>%
    ungroup() %>%
    mutate(value = value / 100)
  weight_matrix <- matrix(NA, length(regs), length(dd))
  for (i in seq_along(abbs)) {
    reg_weights <- w_expanded %>% filter(Abbrev == abbs[i])
    weight_matrix[i, ] <- date_df %>%
      left_join(reg_weights, by = "year") %>%
      pull(value)
  }
  weight_matrix
}

# Factor Initialization ----
factors_em <- function(x, kmax, jj, DEMEAN, NFAC) {
  if (sum(rowSums(is.na(x)) == ncol(x)) > 0) {
    stop("Input x contains entire row of missing values.")
  }
  if (sum(colSums(is.na(x)) == nrow(x)) > 0) {
    stop("Input x contains entire column of missing values.")
  }
  if (!((kmax <= ncol(x) && kmax >= 1 && floor(kmax) == kmax) || kmax == 99)) {
    stop("Input kmax is specified incorrectly.")
  }
  if (!(jj %in% c(1, 2, 3))) {
    stop("Input jj is specified incorrectly.")
  }
  if (!(DEMEAN %in% c(0, 1, 2, 3))) {
    stop("Input DEMEAN is specified incorrectly.")
  }
  maxit <- 10000
  T <- nrow(x)
  N <- ncol(x)
  err <- 999
  it <- 0
  x1 <- is.na(x)
  mut <- matrix(rep(colMeans(x, na.rm = TRUE), T), nrow = T, byrow = TRUE)
  x2 <- x
  x2[is.na(x)] <- mut[is.na(x)]
  transform_result <- transform_data(x2, DEMEAN)
  x3 <- transform_result$x22
  mut <- transform_result$mut
  sdt <- transform_result$sdt
  if (kmax != 99) {
    icstar <- baing(x3, kmax, jj)$ic1
  } else {
    icstar <- NFAC
  }
  pc_result <- pc2(x3, icstar)
  chat <- pc_result$chat
  Fhat <- pc_result$fhat
  lamhat <- pc_result$lambda
  ve2 <- pc_result$ss
  chat0 <- chat
  while (err > 0.000001 && it < maxit) {
    x2 <- chat * sdt + mut
    x2[!x1] <- x[!x1]
    transform_result <- transform_data(x2, DEMEAN)
    x3 <- transform_result$x22
    mut <- transform_result$mut
    sdt <- transform_result$sdt
    pc_result <- pc2(x3, icstar)
    chat <- pc_result$chat
    Fhat <- pc_result$fhat
    lamhat <- pc_result$lambda
    ve2 <- pc_result$ss
    err <- mean((chat - chat0)^2, na.rm = TRUE)
    chat0 <- chat
    it <- it + 1
  }
  x2 <- chat * sdt + mut
  x2[!x1] <- x[!x1]
  ehat <- x - x2
  return(list(ehat = ehat, Fhat = Fhat, lamhat = lamhat, ve2 = ve2, x2 = x2))
}
baing <- function(X, kmax, jj) {
  T <- nrow(X)
  N <- ncol(X)
  NT <- N * T
  NT1 <- N + T
  CT <- numeric(kmax)
  ii <- 1:kmax
  GCT <- min(N, T)
  if (jj == 1) {
    CT <- log(NT / NT1) * ii * NT1 / NT
  } else if (jj == 2) {
    CT <- (NT1 / NT) * log(min(N, T)) * ii
  } else if (jj == 3) {
    CT <- ii * log(GCT) / GCT
  }
  if (N < T) {
    ev_result <- eigen(crossprod(X))
    eigval <- ev_result$values
    Fhat0 <- X %*% ev_result$vectors
    Lambda0 <- ev_result$vectors
  } else {
    ev_result <- eigen(tcrossprod(X))
    eigval <- ev_result$values
    Lambda0 <- crossprod(X, ev_result$vectors) / sqrt(T)
    Fhat0 <- ev_result$vectors * sqrt(T)
  }
  Sigma <- numeric(kmax)
  for (i in 1:kmax) {
    Fhat <- Fhat0[, 1:i, drop = FALSE]
    lambda <- Lambda0[, 1:i, drop = FALSE]
    chat <- Fhat %*% t(lambda)
    ehat <- X - chat
    Sigma[i] <- mean(ehat^2)
  }
  IC1 <- log(Sigma) + CT
  ic1 <- which.min(IC1)
  ic1 <- ifelse(ic1 <= kmax, ic1, kmax)
  Fhat <- Fhat0[, 1:kmax, drop = FALSE]
  Lambda <- Lambda0[, 1:kmax, drop = FALSE]
  chat <- Fhat %*% t(Lambda)
  return(list(ic1 = ic1, chat = chat, Fhat = Fhat, eigval = eigval))
}
pc2 <- function(X, nfac) {
  N <- ncol(X)
  svd_result <- svd(crossprod(X))
  U <- svd_result$u
  S <- svd_result$d
  lambda <- U[, 1:nfac, drop = FALSE] * sqrt(N)
  fhat <- X %*% lambda / N
  chat <- fhat %*% t(lambda)
  ss <- S
  return(list(chat = chat, fhat = fhat, lambda = lambda, ss = ss))
}
transform_data <- function(x2, DEMEAN) {
  T <- nrow(x2)
  N <- ncol(x2)
  if (DEMEAN == 0) {
    mut <- matrix(0, T, N)
    sdt <- matrix(1, T, N)
    x22 <- x2
  } else if (DEMEAN == 1) {
    mut <- matrix(rep(colMeans(x2), T), nrow = T, byrow = TRUE)
    sdt <- matrix(1, T, N)
    x22 <- x2 - mut
  } else if (DEMEAN == 2) {
    mut <- matrix(rep(colMeans(x2), T), nrow = T, byrow = TRUE)
    sdt <- matrix(rep(apply(x2, 2, sd), T), nrow = T, byrow = TRUE)
    x22 <- (x2 - mut) / sdt
  } else if (DEMEAN == 3) {
    mut <- matrix(NA, nrow = T, ncol = N)
    for (t in 1:T) {
      mut[t, ] <- colMeans(x2[1:t, , drop = FALSE])
    }
    sdt <- matrix(rep(apply(x2, 2, sd), T), nrow = T, byrow = TRUE)
    x22 <- (x2 - mut) / sdt
  }
  return(list(x22 = x22, mut = mut, sdt = sdt))
}

# MCMC Samplers ----
GEN_PHI <- function(ztt, Sigma0, B0, sigma2, lag) {
  Y <- ztt
  X <- NULL
  for (i in 1:lag) {
    X <- cbind(X, lag0(matrix(Y, ncol = 1), i))
  }
  Y <- Y[(lag + 1):length(Y)]
  X <- X[(lag + 1):nrow(X), , drop = FALSE]
  XtX <- crossprod(X)
  XtY <- crossprod(X, Y)
  A <- solve(Sigma0) + (1 / sigma2) * XtX
  V <- solve(A)
  M <- V %*% (solve(Sigma0) %*% B0 + (1 / sigma2) * XtY)
  chck <- -1
  while (chck < 0) {
    B <- M + t(matrix(rnorm(lag), 1, lag) %*% safeChol(V))
    b <- rbind(t(B), cbind(diag(lag - 1), matrix(0, lag - 1, 1)))
    ee <- max(abs(eigen(b)$values))
    if (ee <= 1) {
      chck <- 1
    }
  }
  PHI <- t(B)
  return(PHI)
}
GEN_SIG2F <- function(ztt, phi, V0f, D0f, lag) {
  Y <- ztt[(lag + 1):length(ztt)]
  X <- NULL
  for (i in 1:lag) {
    X <- cbind(X, lag0(matrix(ztt, ncol = 1), i))
  }
  X <- X[(lag + 1):nrow(X), , drop = FALSE]
  resids <- Y - X %*% matrix(phi, ncol = 1)
  T_eff <- length(resids)
  sig2f <- 1 / rgamma(1, V0f + T_eff / 2, rate = D0f + sum(resids^2) / 2)
  return(sig2f)
}
GEN_PSI_AND_SIG_ROB <- function(
  Y0,
  sigma2,
  Sigma0,
  B0,
  T0,
  D0,
  lag,
  index,
  PSI_old
) {
  ind <- which(index == 1)[1]
  Y <- Y0[ind:length(Y0)]
  index1 <- index[ind:length(index)]
  X <- NULL
  for (i in 1:lag) {
    X <- cbind(X, lag0(matrix(Y, ncol = 1), i))
  }
  Y <- Y[(lag + 1):length(Y)]
  X <- X[(lag + 1):nrow(X), , drop = FALSE]
  index2 <- index1[(lag + 1):length(index1)]
  T <- nrow(X)
  XtX <- crossprod(X)
  XtY <- crossprod(X, Y)
  A <- solve(Sigma0) + (1 / sigma2) * XtX
  V <- solve(A)
  M <- V %*% (solve(Sigma0) %*% B0 + (1 / sigma2) * XtY)
  chck <- -1
  cou <- 0
  while (chck < 0) {
    B <- M + t(matrix(rnorm(lag), 1, lag) %*% safeChol(V))
    b <- rbind(t(B), cbind(diag(lag - 1), matrix(0, lag - 1, 1)))
    ee <- max(abs(eigen(b)$values))
    if (ee <= 1) {
      chck <- 1
    } else {
      cou <- cou + 1
      if (cou > 10000) {
        B <- matrix(PSI_old, ncol = 1)
        chck <- 1
      }
    }
  }
  PSI <- t(B)
  resids0 <- Y - tcrossprod(X, PSI)
  resids <- delif(resids0, 1 - index2)
  sigma2 <- 1 /
    rgamma(1, T0 + length(resids) / 2, rate = D0 + sum(resids^2) / 2)
  return(list(PSI = PSI, sigma2 = sigma2))
}
GEN_LAMDA <- function(yy0, zz, psi, sigma2, Sigma0, B0, kernel = 1, index0 = NULL) {
  yy0 <- t(yy0)
  Lu <- length(psi)
  TT0 <- length(yy0)
  YM <- yy0[(Lu + 1):TT0]
  X <- zz[(Lu + 1):length(zz)]
  for (j in 1:Lu) {
    YM <- YM - psi[j] * yy0[(Lu + 1 - j):(TT0 - j)]
    X <- X - psi[j] * zz[(Lu + 1 - j):(length(zz) - j)]
  }
  kLen <- length(kernel)
  if (kLen > 1) {
    TT <- length(X)
    ff <- matrix(X, ncol = 1)
    for (i in 1:(kLen - 1)) {
      ff <- cbind(ff, lag0(matrix(X, ncol = 1), i))
    }
    FM <- rep(NA_real_, TT)
    for (i in (kLen + 1):TT) {
      FM[i] <- sum(ff[i, ] * kernel)
    }
    X <- FM[(kLen + 1):TT]
    YM <- YM[(kLen + 1):length(YM)]
    if (!is.null(index0)) {
      index0 <- index0[(Lu + kLen + 1):length(index0)]
    }
  } else if (!is.null(index0)) {
    index0 <- index0[(Lu + 1):length(index0)]
  }
  if (!is.null(index0)) {
    ind <- which(index0 == 1)[1]
    X <- X[ind:length(X)]
    YM <- YM[ind:length(YM)]
  }
  Hpsid <- sigma2 * sum(kernel^2)
  fw_star <- X / sqrt(Hpsid)
  yw_star <- YM / sqrt(Hpsid)
  ftf <- crossprod(fw_star)
  fty <- crossprod(fw_star, yw_star)
  A <- solve(Sigma0) + ftf
  V <- solve(A)
  M <- V %*% (solve(Sigma0) %*% B0 + fty)
  M + t(matrix(rnorm(1), 1, 1) %*% safeChol(V))
}
# Precision Sampler Helpers ----


buildAggMatrix <- function(kernel, obsIndex, TT) {
  kLen <- length(kernel)
  nObs <- length(obsIndex)
  grid <- expand.grid(r = seq_len(nObs), k = seq_along(kernel))
  grid$col <- obsIndex[grid$r] - grid$k + 1L
  grid$x <- kernel[grid$k]
  valid <- grid$col >= 1L & grid$col <= TT
  sparseMatrix(
    i = grid$r[valid],
    j = grid$col[valid],
    x = grid$x[valid],
    dims = c(nObs, TT)
  )
}


precomputeObsCache <- function(dat, index0, aggs, kernels, N, TT_ext, sigma2z, pad) {
  invR <- 1 / sigma2z
  cache <- vector("list", N)
  for (i in 1:N) {
    kernel <- kernels[[aggs[i]]]
    obsIdx <- which(index0[, i] == 1)
    yObs <- as.numeric(dat[i, obsIdx])
    obsIdxShifted <- obsIdx + pad
    M <- buildAggMatrix(kernel, obsIdxShifted, TT_ext)
    MtM <- crossprod(M)
    MtY <- as.numeric(t(M) %*% yObs)
    MtM_g <- as(invR * MtM, "generalMatrix")
    sp <- summary(MtM_g)
    cache[[i]] <- list(
      MtY_invR = invR * MtY,
      sp_i = sp$i,
      sp_j = sp$j,
      sp_x = sp$x
    )
  }
  cache
}


computeYHAT <- function(f, tau, U, lamd, aggs, kernels, N, TT) {
  YHAT <- matrix(0, N, TT)
  for (i in 1:N) {
    kernel <- kernels[[aggs[i]]]
    kLen <- length(kernel)
    signal <- lamd[i] * f + U[, i]
    if (i == 1) signal <- signal + tau
    for (t in kLen:TT) {
      YHAT[i, t] <- sum(kernel * signal[(t - kLen + 1):t])
    }
    if (kLen > 1) {
      for (t in 1:(kLen - 1)) {
        jmax <- min(t, kLen) - 1
        YHAT[i, t] <- sum(kernel[1:(jmax + 1)] * signal[t:(t - jmax)])
      }
    }
  }
  YHAT
}

# Kalman Filter / State Space ----
GEN_STATE_VECTOR <- function(
  dat,
  phi,
  lamd,
  psi,
  sig2,
  sig2f,
  lag,
  lagu,
  index0,
  n,
  sig2_tau,
  a1 = NULL
) {
  tstar <- ncol(dat)
  NY <- n[1]
  NQ <- n[2]
  NM <- n[3]
  NM2 <- n[4]
  lamd_Y <- lamd[1:NY]
  lamd_Q <- lamd[(NY + 1):(NY + NQ)]
  lamd_M <- lamd[(NY + NQ + 1):(NY + NQ + NM)]
  lamd_M2 <- lamd[(NY + NQ + NM + 1):length(lamd)]
  y23s <- c(1:12, 11:1)
  y23a <- c(1:12, 11:1) / 12
  q14s <- c(1, 2, rep(3, 10), 2, 1)
  q14a <- c(1, 2, rep(3, 10), 2, 1) / 3
  m12 <- rep(1, 12)
  m1 <- 1
  yy <- y23a
  qq <- q14a
  mm <- m12
  mm2 <- m1
  NS <- length(yy) +
    NY * length(yy) +
    NQ * length(qq) +
    NM * length(mm) +
    NM2 * length(mm2) * lagu
  YY <- cbind(
    lamd_Y * matrix(rep(yy, each = NY), NY),
    kronecker(diag(NY), t(yy)),
    matrix(0, NY, length(qq) * NQ),
    matrix(0, NY, length(mm) * NM),
    matrix(0, NY, length(mm2) * lagu * NM2)
  )
  QQ <- cbind(
    lamd_Q * matrix(rep(qq, each = NQ), NQ),
    matrix(0, NQ, length(yy) - length(qq)),
    matrix(0, NQ, length(yy) * NY),
    kronecker(diag(NQ), t(qq)),
    matrix(0, NQ, length(mm) * NM),
    matrix(0, NQ, length(mm2) * lagu * NM2)
  )
  M0 <- cbind(
    lamd_M * matrix(rep(mm, each = NM), NM),
    matrix(0, NM, length(yy) - length(mm) + length(yy) * NY + length(qq) * NQ),
    kronecker(diag(NM), matrix(rep(1, length(mm)), nrow = 1)),
    matrix(0, NM, length(mm2) * lagu * NM2)
  )
  if (NM2 > 0) {
    M02 <- cbind(
      lamd_M2 * matrix(rep(mm2, each = NM2), NM2),
      matrix(
        0,
        NM2,
        length(yy) -
          length(mm2) +
          length(yy) * NY +
          length(qq) * NQ +
          length(mm) * NM
      ),
      kronecker(diag(NM2), t(c(1, rep(0, lagu - 1))))
    )
  } else {
    M02 <- matrix(ncol = NS, nrow = 0)
  }
  H <- rbind(YY, QQ, M0, M02)
  F_mat <- matrix(0, NS, NS)
  F_mat[1, 1:lag] <- phi
  if (length(yy) > 1) {
    F_mat[2:length(yy), 1:(length(yy) - 1)] <- diag(length(yy) - 1)
  }
  ini_y_r <- length(yy) + 1
  fin_y_r <- length(yy) * 2
  ini_y_c <- length(yy) + 1
  fin_y_c <- length(yy) * 2 - 1
  for (iy in 1:NY) {
    F_mat[ini_y_r:fin_y_r, ini_y_c:fin_y_c] <- rbind(
      c(psi[iy, ], rep(0, length(yy) - 1 - lagu)),
      diag(length(yy) - 1)
    )
    if (iy < NY) {
      ini_y_r <- ini_y_r + length(yy)
      fin_y_r <- fin_y_r + length(yy)
      ini_y_c <- ini_y_c + length(yy)
      fin_y_c <- fin_y_c + length(yy)
    }
  }
  ini_q_r <- fin_y_r + 1
  fin_q_r <- fin_y_r + length(qq)
  ini_q_c <- fin_y_c + 2
  fin_q_c <- fin_y_c + length(qq)
  for (iq in (NY + 1):(NY + NQ)) {
    F_mat[ini_q_r:fin_q_r, ini_q_c:fin_q_c] <- rbind(
      c(psi[iq, ], rep(0, length(qq) - 1 - lagu)),
      diag(length(qq) - 1)
    )
    if (iq < NY + NQ) {
      ini_q_r <- ini_q_r + length(qq)
      fin_q_r <- fin_q_r + length(qq)
      ini_q_c <- ini_q_c + length(qq)
      fin_q_c <- fin_q_c + length(qq)
    }
  }
  ini_m_r <- fin_q_r + 1
  fin_m_r <- fin_q_r + length(mm)
  ini_m_c <- fin_q_c + 2
  fin_m_c <- fin_q_c + length(mm)
  for (im in (NY + NQ + 1):(NY + NQ + NM)) {
    F_mat[ini_m_r:fin_m_r, ini_m_c:fin_m_c] <- rbind(
      c(psi[im, ], rep(0, length(mm) - 1 - lagu)),
      diag(length(mm) - 1)
    )
    if (im < NY + NQ + NM) {
      ini_m_r <- ini_m_r + length(mm)
      fin_m_r <- fin_m_r + length(mm)
      ini_m_c <- ini_m_c + length(mm)
      fin_m_c <- fin_m_c + length(mm)
    }
  }
  if (NM2 > 0) {
    ini_m2_r <- fin_m_r + 1
    fin_m2_r <- fin_m_r + lagu
    ini_m2_c <- fin_m_c + 2
    fin_m2_c <- fin_m_c + lagu + 1
    for (im in (NY + NQ + NM + 1):(NY + NQ + NM + NM2)) {
      F_mat[ini_m2_r:fin_m2_r, ini_m2_c:fin_m2_c] <- rbind(
        psi[im, ],
        cbind(diag(lagu - 1), matrix(0, lagu - 1, 1))
      )
      if (im < NY + NQ + NM + NM2) {
        ini_m2_r <- ini_m2_r + lagu
        fin_m2_r <- fin_m2_r + lagu
        ini_m2_c <- ini_m2_c + lagu
        fin_m2_c <- fin_m2_c + lagu
      }
    }
  }
  Q <- matrix(0, NS, NS)
  Q[1, 1] <- sig2f
  ENT <- 1
  entry_y <- length(yy) + 1
  for (iy in 1:NY) {
    Q[entry_y, entry_y] <- sig2[iy]
    ENT <- c(ENT, entry_y)
    if (iy < NY) {
      entry_y <- entry_y + length(yy)
    }
  }
  entry_q <- entry_y + length(yy)
  for (im in (NY + 1):(NY + NQ)) {
    Q[entry_q, entry_q] <- sig2[im]
    ENT <- c(ENT, entry_q)
    if (im < NY + NQ) {
      entry_q <- entry_q + length(qq)
    }
  }
  entry_m <- entry_q + length(qq)
  for (im in (NY + NQ + 1):(NY + NQ + NM)) {
    Q[entry_m, entry_m] <- sig2[im]
    ENT <- c(ENT, entry_m)
    if (im < NY + NQ + NM) {
      entry_m <- entry_m + length(mm)
    }
  }
  if (NM2 > 0) {
    entry_m2 <- entry_m + length(mm)
    for (im in (NY + NQ + NM + 1):(NY + NQ + NM + NM2)) {
      Q[entry_m2, entry_m2] <- sig2[im]
      ENT <- c(ENT, entry_m2)
      if (im < NY + NQ + NM + NM2) {
        entry_m2 <- entry_m2 + lagu
      }
    }
  }
  blockStarts <- 0L
  pos <- length(yy)
  for (iy in seq_len(NY)) {
    blockStarts <- c(blockStarts, pos)
    pos <- pos + length(yy)
  }
  for (iq in seq_len(NQ)) {
    blockStarts <- c(blockStarts, pos)
    pos <- pos + length(qq)
  }
  for (im in seq_len(NM)) {
    blockStarts <- c(blockStarts, pos)
    pos <- pos + length(mm)
  }
  if (NM2 > 0) {
    for (im in seq_len(NM2)) {
      blockStarts <- c(blockStarts, pos)
      pos <- pos + lagu
    }
  }
  NS_old <- NS
  nTau <- length(yy)
  NS <- NS + nTau
  H <- cbind(H, matrix(0, nrow(H), nTau))
  H[1, (NS_old + 1):NS] <- yy
  F_new <- matrix(0, NS, NS)
  F_new[1:NS_old, 1:NS_old] <- F_mat
  F_new[NS_old + 1, NS_old + 1] <- 1
  if (nTau > 1) {
    F_new[(NS_old + 2):NS, (NS_old + 1):(NS - 1)] <- diag(nTau - 1)
  }
  F_mat <- F_new
  Q_new <- matrix(0, NS, NS)
  Q_new[1:NS_old, 1:NS_old] <- Q
  Q_new[NS_old + 1, NS_old + 1] <- sig2_tau
  Q <- Q_new
  ENT <- c(ENT, NS_old + 1)
  blockStarts <- c(blockStarts, NS_old)
  if (is.null(a1)) {
    a1 <- rep(0, NS)
  }
  res <- dk_simulation_smoother_cpp(
    H, F_mat, Q, dat, index0, ENT, blockStarts, a1
  )
  res$F_mat <- F_mat
  res$H <- H
  res$NS_base <- NS_old
  return(res)
}

# Estimation ----
prepareRegionData <- function(sourceData, regionName, vargroups) {
  varnames <- colnames(sourceData)
  regCols <- grep(regionName, varnames, fixed = TRUE)
  metaCols <- c("date", "year", "quarter")
  natCols <- which(!grepl("\\.", varnames) & !varnames %in% metaCols)
  tabdat_cols <- sort(union(regCols, natCols))
  tabdat <- sourceData[, tabdat_cols]
  tabdat <- tabdat[, as.numeric(tabdat[3, ]) %in% vargroups]
  allNA <- sapply(tabdat[4:nrow(tabdat), ], function(col) {
    all(is.na(as.numeric(col)))
  })
  tabdat <- tabdat[, !allNA]
  tabdat <- tabdat[, order(as.numeric(tabdat[2, ]))]
  tabdat
}
settingsToTable <- function(settings) {
  desc <- c(
    hor = "Forecast horizon (months)",
    M0 = "Burn-in iterations",
    M = "Post-burn-in iterations",
    L = "Factor AR lags",
    Lu = "Idiosyncratic AR lags",
    T00 = "Loading prior mean",
    R00 = "Loading prior variance",
    T0 = "Factor AR prior mean",
    R0 = "Factor AR prior variance",
    T0p = "Idiosyncratic AR prior mean",
    R0p = "Idiosyncratic AR prior variance",
    V0 = "Idiosyncratic variance IG shape",
    D0 = "Idiosyncratic variance IG rate",
    V0f = "Factor variance IG shape",
    D0f = "Factor variance IG rate",
    fix_lam1 = "Fix first loading to 1",
    fix_sig2f = "Fix factor error variance",
    fix_psi_grp = "Fix GRP idiosyncratic AR (FALSE/TRUE/numeric)",
    sig2_tau = "Tau innovation variance (fixed)",
    smoother = "Simulation smoother method",
    estim_start = "Estimation start date",
    QQ = "Credible interval quantiles",
    progressEvery = "Progress report interval"
  )
  params <- intersect(names(desc), names(settings))
  tibble(
    parameter = params,
    description = desc[params],
    value = sapply(params, function(p) paste(settings[[p]], collapse = ", "))
  )
}
estimRegionTarget <- function(
  regionName,
  sourceData,
  estSettings,
  vargroups,
  folder,
  grpData,
  stdParams = NULL,
  regionOverrides = list()
) {
  if (regionName %in% names(regionOverrides)) {
    estSettings <- modifyList(estSettings, regionOverrides[[regionName]])
  }
  outpath <- paste0(folder, "/FACTORS", regionName, ".xlsx")
  statusFile <- paste0(folder, "/.status_", regionName)
  if (file.exists(outpath)) {
    if (file.exists(statusFile)) {
      file.remove(statusFile)
    }
    return(outpath)
  }
  writeLines(sprintf("%s|starting|0|0|%s", regionName, Sys.time()), statusFile)
  # WORKER_RCPP_CACHE: per-worker cache dir to avoid races between
  # parallel crew workers all writing to a shared .rcpp_cache.
  .cacheDir <- file.path(tempdir(), ".rcpp_cache")
  dir.create(.cacheDir, showWarnings = FALSE, recursive = TRUE)
  Rcpp::sourceCpp("R/functions.cpp", cacheDir = .cacheDir)
  Rcpp::sourceCpp("R/precision_sampler.cpp", cacheDir = .cacheDir)
  tabdat <- prepareRegionData(sourceData, regionName, vargroups)
  extraSheets <- list()
  if (!is.null(stdParams)) {
    regionVars <- colnames(tabdat)
    extraSheets$standardization <- stdParams %>% filter(var %in% regionVars)
  }
  extraSheets$settings <- settingsToTable(estSettings)
  estim_region(
    regionName,
    tabdat,
    estSettings,
    folder,
    extraSheets = extraSheets
  )
  file.remove(statusFile)
  outpath
}
estim_region <- function(
  sss,
  tabdat,
  settings,
  folder,
  extraSheets = list()
) {
  M0 <- settings$M0
  M <- settings$M
  capn <- M0 + M
  aggs <- as.numeric(tabdat[2, ])
  groups <- as.numeric(tabdat[3, ])
  data <- as.matrix(tabdat[4:nrow(tabdat), ])
  timewF <- seq(settings$estim_start, by = "month", length.out = nrow(data))
  yy <- t(data)
  NY <- sum(aggs == 1)
  NQ <- sum(aggs == 2)
  NM <- sum(aggs == 3)
  NM2 <- sum(aggs == 4)
  n <- c(NY, NQ, NM, NM2)
  N <- nrow(yy)
  T <- ncol(yy)
  L <- settings$L
  Lu <- settings$Lu
  hor <- settings$hor
  index <- 1 - is.na(data)
  QQ <- settings$QQ
  T00 <- settings$T00
  R00 <- settings$R00
  T0 <- c(settings$T0, rep(0, L - 1))
  R0 <- diag(settings$R0 / (1:L)^2)
  T0p <- rep(settings$T0p, Lu)
  R0p <- diag(settings$R0p / (1:Lu)^2)
  V0 <- settings$V0
  D0 <- settings$D0
  fix_lam1 <- settings$fix_lam1
  fix_sig2f <- isTRUE(settings$fix_sig2f)
  fix_psi_grp <- settings$fix_psi_grp
  if (is.null(fix_psi_grp) || identical(fix_psi_grp, FALSE)) {
    fix_psi_grp <- FALSE
  } else if (isTRUE(fix_psi_grp)) {
    fix_psi_grp <- 0
  }
  smoother <- if (!is.null(settings$smoother)) settings$smoother else "dk"
  V0f <- settings$V0f
  D0f <- settings$D0f
  sig2_tau <- settings$sig2_tau
  grpVals <- as.numeric(yy[1, ])
  grpVals <- grpVals[!is.na(grpVals)]
  iniTau <- mean(grpVals) / 12
  em_result <- factors_em(t(yy), 99, 1, 0, 1)
  ini_fac <- em_result$Fhat
  ini_load <- em_result$lamhat
  if (ini_load[1] < 0) {
    LAMDTT <- -ini_load
  } else {
    LAMDTT <- ini_load
  }
  PHITT <- rep(0.1, L)
  PSITT <- matrix(0.1, N, Lu)
  SIG2TT <- rep(0.1, N)
  if (!identical(fix_psi_grp, FALSE)) {
    PSITT[1, ] <- 0
    SIG2TT[1] <- fix_psi_grp
  }
  SIG2FTT <- 1
  ZTTMM <- matrix(0, T, M)
  LAMDMM <- matrix(0, N, M)
  GRP_DRAWS <- matrix(0, T + hor, M)
  PHITTMM <- matrix(0, L, M)
  PSITTMM <- matrix(0, N * Lu, M)
  SIG2TTMM <- matrix(0, N, M)
  SIG2FTTMM <- rep(NA_real_, M)
  TAUTTMM <- matrix(0, T, M)
  progressEvery <- settings$progressEvery
  if (is.null(progressEvery)) {
    progressEvery <- 0
  }
  statusFile <- paste0(folder, "/.status_", sss)
  writeLines(sprintf("%s|running|0|%d|%s", sss, capn, Sys.time()), statusFile)
  kernels <- list(
    `1` = c(1:12, 11:1) / 12, `2` = c(1, 2, rep(3, 10), 2, 1) / 3,
    `3` = rep(1, 12), `4` = 1
  )
  if (smoother == "precision") {
    maxKernelLen <- max(sapply(kernels[as.character(unique(aggs))], length))
    pad <- maxKernelLen - 1L
    TT_ext <- T + pad
    sigma2z <- 1e-8
    obsCache <- precomputeObsCache(
      yy, index, aggs, kernels, N, TT_ext, sigma2z, pad
    )
    sIdx <- (pad + 1):TT_ext
    y23a_prec <- kernels[["1"]]
  }
  for (itr in 1:capn) {
    if (progressEvery > 0 && itr %% progressEvery == 0) {
      writeLines(
        sprintf("%s|running|%d|%d|%s", sss, itr, capn, Sys.time()),
        statusFile
      )
    }
    if (smoother == "precision") {
      mT <- (2L + N) * TT_ext
      z <- precisionSampleCpp(
        PHITT, PSITT, SIG2TT, SIG2FTT,
        sig2_tau, TRUE, N, TT_ext,
        obsCache, LAMDTT, rnorm(mT), 1.0
      )
      ZTT <- z[sIdx]
      offIdx <- TT_ext
      tauTT <- z[offIdx + sIdx]
      offIdx <- offIdx + TT_ext
      UTT <- matrix(0, T, N)
      for (i in 1:N) {
        UTT[, i] <- z[offIdx + sIdx]
        offIdx <- offIdx + TT_ext
      }
      fFull <- z[1:TT_ext]
      tauFull <- z[(TT_ext + 1):(2 * TT_ext)]
      UFull <- matrix(0, TT_ext, N)
      offIdx2 <- 2L * TT_ext
      for (i in 1:N) {
        UFull[, i] <- z[(offIdx2 + 1):(offIdx2 + TT_ext)]
        offIdx2 <- offIdx2 + TT_ext
      }
      YHAT <- computeYHAT(
        fFull, tauFull, UFull, LAMDTT, aggs, kernels, N, TT_ext
      )
      YHAT <- YHAT[, sIdx, drop = FALSE]
    } else {
      nTau <- 23
      NS_base <- nTau + NY * nTau + NQ * 14 + NM * 12 + NM2 * Lu
      a1 <- rep(0, NS_base + nTau)
      a1[(NS_base + 1):(NS_base + nTau)] <- iniTau
      state_result <- GEN_STATE_VECTOR(
        yy, PHITT, LAMDTT, PSITT, SIG2TT, SIG2FTT,
        L, Lu, index, n, sig2_tau, a1
      )
      ZTTall <- state_result$Z_MAT
      YHAT <- state_result$YHAT
      km <- state_result$km
      F_mat <- state_result$F_mat
      H_mat <- state_result$H
      NS_base <- state_result$NS_base
      ZTT <- ZTTall[, 1]
      nKm <- length(km)
      tauTT <- ZTTall[, NS_base + 1]
      UTT <- ZTTall[, km[2:(nKm - 1)]]
    }
    for (i in 1:N) {
      if (!identical(fix_psi_grp, FALSE) && i == 1) {
        next
      }
      psi_result <- GEN_PSI_AND_SIG_ROB(
        UTT[, i],
        SIG2TT[i],
        R0p,
        T0p,
        V0,
        D0,
        Lu,
        index[, i],
        PSITT[i, ]
      )
      PSITT[i, ] <- psi_result$PSI
      SIG2TT[i] <- psi_result$sigma2
    }
    LAMDTT_old <- LAMDTT
    cou <- 0
    for (i in 1:N) {
      if (i == 1 && fix_lam1) {
        LAMDTT[i] <- 1
        next
      }
      idx0 <- if (aggs[i] %in% c(3, 4)) index[, i] else NULL
      if (i == 1) {
        chk <- -1
        while (chk < 0) {
          LAMDTT[i] <- GEN_LAMDA(
            YHAT[i, ], ZTT, PSITT[i, ], SIG2TT[i],
            R00, T00, kernels[[aggs[i]]], idx0
          )
          if (LAMDTT[i] > 0) {
            chk <- 1
          } else {
            cou <- cou + 1
            if (cou > 10000) {
              LAMDTT[i] <- LAMDTT_old[i]
              chk <- 1
            }
          }
        }
      } else {
        LAMDTT[i] <- GEN_LAMDA(
          YHAT[i, ], ZTT, PSITT[i, ], SIG2TT[i],
          R00, T00, kernels[[aggs[i]]], idx0
        )
      }
    }
    PHITT <- GEN_PHI(ZTT, R0, T0, SIG2FTT, L)
    if (!fix_sig2f) {
      SIG2FTT <- GEN_SIG2F(ZTT, PHITT, V0f, D0f, L)
    }
    if (itr > M0) {
      idx <- itr - M0
      ZTTMM[, idx] <- ZTT
      LAMDMM[, idx] <- LAMDTT
      if (smoother == "precision") {
        grpMonthly <- LAMDTT[1] * ZTT + UTT[, 1] + tauTT
        for (t in length(y23a_prec):T) {
          GRP_DRAWS[t, idx] <- sum(y23a_prec * grpMonthly[(t - length(y23a_prec) + 1):t])
        }
        fExt <- c(ZTT, rep(0, hor))
        tauExt <- c(tauTT, rep(tauTT[T], hor))
        for (h in 1:hor) {
          fExt[T + h] <- sum(PHITT * fExt[(T + h - 1):(T + h - L)])
        }
        grpExt <- LAMDTT[1] * fExt + tauExt
        grpExt[1:T] <- grpExt[1:T] + UTT[, 1]
        for (h in 1:hor) {
          t <- T + h
          if (t >= length(y23a_prec)) {
            GRP_DRAWS[t, idx] <- sum(y23a_prec * grpExt[(t - length(y23a_prec) + 1):t])
          }
        }
      } else {
        GRP_DRAWS[1:T, idx] <- as.numeric(ZTTall %*% H_mat[1, ])
        S_curr <- ZTTall[T, ]
        for (h in seq_len(hor)) {
          S_curr <- as.numeric(F_mat %*% S_curr)
          GRP_DRAWS[T + h, idx] <- sum(H_mat[1, ] * S_curr)
        }
      }
      PHITTMM[, idx] <- PHITT
      PSITTMM[, idx] <- c(PSITT)
      SIG2TTMM[, idx] <- c(SIG2TT)
      SIG2FTTMM[idx] <- SIG2FTT
      TAUTTMM[, idx] <- tauTT
    }
  }
  varnames <- colnames(tabdat)
  hpdProb <- QQ[2] - QQ[1]
  ufactor <- apply(ZTTMM, 1, median)
  ufactor_hpd <- apply(ZTTMM, 1, hpdInterval, prob = hpdProb)
  ufactor_ql <- sapply(ufactor_hpd, `[[`, "l")
  ufactor_qu <- sapply(ufactor_hpd, `[[`, "u")
  futureDates <- seq(max(timewF) + months(1), by = "month", length.out = hor)
  allDates <- c(timewF, futureDates)
  grpFitted <- apply(GRP_DRAWS, 1, median)
  grpHpd <- apply(GRP_DRAWS, 1, hpdInterval, prob = hpdProb)
  grpQl <- sapply(grpHpd, `[[`, "l")
  grpQu <- sapply(grpHpd, `[[`, "u")
  lambda <- apply(LAMDMM, 1, median)
  lambdaHpd <- apply(LAMDMM, 1, hpdInterval, prob = hpdProb)
  lambda_ql <- sapply(lambdaHpd, `[[`, "l")
  lambda_qu <- sapply(lambdaHpd, `[[`, "u")
  phi <- apply(PHITTMM, 1, median)
  phiHpd <- apply(PHITTMM, 1, hpdInterval, prob = hpdProb)
  phi_ql <- sapply(phiHpd, `[[`, "l")
  phi_qu <- sapply(phiHpd, `[[`, "u")
  psi <- apply(PSITTMM, 1, median)
  psiHpd <- apply(PSITTMM, 1, hpdInterval, prob = hpdProb)
  psi_ql <- sapply(psiHpd, `[[`, "l")
  psi_qu <- sapply(psiHpd, `[[`, "u")
  sig <- apply(SIG2TTMM, 1, median)
  sigHpd <- apply(SIG2TTMM, 1, hpdInterval, prob = hpdProb)
  sig_ql <- sapply(sigHpd, `[[`, "l")
  sig_qu <- sapply(sigHpd, `[[`, "u")
  factors_df <- data.frame(
    date = as.character(timewF),
    factor = ufactor,
    l = ufactor_ql,
    u = ufactor_qu
  )
  utau <- apply(TAUTTMM, 1, median)
  utauHpd <- apply(TAUTTMM, 1, hpdInterval, prob = hpdProb)
  utau_ql <- sapply(utauHpd, `[[`, "l")
  utau_qu <- sapply(utauHpd, `[[`, "u")
  factors_df$tau <- utau
  factors_df$tau_l <- utau_ql
  factors_df$tau_u <- utau_qu
  factors_df$grpMoM <- factors_df$factor + factors_df$tau
  factors_df$grpYoY <- zoo::rollsumr(factors_df$grpMoM, 12, fill = NA)
  factors_df$grpBase <- 100 * exp(cumsum(factors_df$grpMoM / 100))
  lambda_df <- data.frame(
    var = varnames,
    lambda = lambda,
    l = lambda_ql,
    u = lambda_qu
  )
  yhat_df <- data.frame(
    date = as.character(allDates),
    grp = grpFitted,
    l = grpQl,
    u = grpQu
  )
  colnames(yhat_df)[2] <- varnames[1]
  phi_df <- data.frame(lag = 1:L, coeff = phi, l = phi_ql, u = phi_qu)
  psi_df <- data.frame(
    var = rep(varnames, Lu),
    lag = rep(1:Lu, each = N),
    coeff = psi,
    l = psi_ql,
    u = psi_qu
  )
  sig_df <- data.frame(var = varnames, sigma2 = sig, l = sig_ql, u = sig_qu)
  if (!fix_sig2f) {
    sig2f_med <- median(SIG2FTTMM)
    sig2fHpd <- hpdInterval(SIG2FTTMM, prob = hpdProb)
    sig2f_ql <- sig2fHpd$l
    sig2f_qu <- sig2fHpd$u
    sig2f_df <- data.frame(
      var = "sig2f",
      sigma2 = sig2f_med,
      l = sig2f_ql,
      u = sig2f_qu
    )
    SIG2FTTMM_df <- data.frame(sig2f = SIG2FTTMM)
  }
  LAMDMM <- data.frame(t(LAMDMM))
  names(LAMDMM) <- varnames
  PHITTMM <- data.frame(t(PHITTMM))
  names(PHITTMM) <- as.character(1:L)
  PSITTMM <- data.frame(t(PSITTMM))
  names(PSITTMM) <- paste0(rep(varnames, times = Lu), '.', rep(1:Lu, each = N))
  SIG2TTMM <- data.frame(t(SIG2TTMM))
  names(SIG2TTMM) <- varnames
  GRP_DRAWS_df <- data.frame(t(GRP_DRAWS))
  names(GRP_DRAWS_df) <- as.character(allDates)
  sheets <- list(
    factors = factors_df,
    lambda = lambda_df,
    fitted = yhat_df,
    arcommon = phi_df,
    aridi = psi_df,
    sig = sig_df,
    data = cbind(date = as.character(timewF), tabdat[-c(1:3), ]),
    lambdamm = LAMDMM,
    phimm = PHITTMM,
    psimm = PSITTMM,
    sigmm = SIG2TTMM,
    grpmm = GRP_DRAWS_df
  )
  if (!fix_sig2f) {
    sheets$sig2f <- sig2f_df
    sheets$sig2fmm <- SIG2FTTMM_df
  }
  sheets$taumm <- data.frame(t(TAUTTMM))
  if (length(extraSheets) > 0) {
    sheets <- c(sheets, extraSheets)
  }
  writexl::write_xlsx(sheets, paste0(folder, "/FACTORS", sss, ".xlsx"))
  writeLines(
    sprintf("%s|done|%d|%d|%s", sss, capn, capn, Sys.time()),
    statusFile
  )
  return(factors_df)
}

# Monitoring ----
monitorProgress <- function(
  folder = "Results/EST_TEST",
  refresh = 5,
  once = FALSE
) {
  repeat {
    files <- list.files(
      folder,
      pattern = "^\\.status_",
      full.names = TRUE,
      all.files = TRUE
    )
    if (length(files) == 0) {
      cat("No status files found in", folder, "\n")
      if (once) {
        return(invisible(NULL))
      }
      Sys.sleep(refresh)
      next
    }
    rows <- lapply(files, function(f) {
      line <- tryCatch(readLines(f, n = 1, warn = FALSE), error = function(e) {
        NA_character_
      })
      if (is.na(line) || length(line) == 0) {
        return(NULL)
      }
      parts <- strsplit(line, "\\|")[[1]]
      if (length(parts) < 5) {
        return(NULL)
      }
      data.frame(
        region = parts[1],
        status = parts[2],
        iter = as.integer(parts[3]),
        total = as.integer(parts[4]),
        time = parts[5],
        stringsAsFactors = FALSE
      )
    })
    df <- do.call(rbind, Filter(Nonnull <- Negate(is.null), rows))
    if (is.null(df) || nrow(df) == 0) {
      cat("No readable status files.\n")
      if (once) {
        return(invisible(NULL))
      }
      Sys.sleep(refresh)
      next
    }
    df$pct <- ifelse(
      df$total > 0,
      sprintf("%5.1f%%", df$iter / df$total * 100),
      ""
    )
    running <- df[df$status == "running", ]
    done <- df[df$status == "done", ]
    other <- df[!df$status %in% c("running", "done"), ]
    if (!once) {
      cat("\033[2J\033[H")
    }
    cat(sprintf("=== Progress Monitor [%s] ===\n", Sys.time()))
    cat(sprintf(
      "Running: %d | Done: %d | Total: %d\n\n",
      nrow(running),
      nrow(done),
      nrow(df)
    ))
    if (nrow(running) > 0) {
      running <- running[order(running$region), ]
      for (i in seq_len(nrow(running))) {
        r <- running[i, ]
        bar <- paste0(
          strrep("#", round(r$iter / max(r$total, 1) * 30)),
          strrep(".", 30 - round(r$iter / max(r$total, 1) * 30))
        )
        cat(sprintf(
          "  %s [%s] %s  iter %d/%d\n",
          r$region,
          bar,
          r$pct,
          r$iter,
          r$total
        ))
      }
      cat("\n")
    }
    if (nrow(done) > 0) {
      cat("Done:", paste(sort(done$region), collapse = ", "), "\n")
    }
    if (nrow(other) > 0) {
      cat(
        "Other:",
        paste(sprintf("%s (%s)", other$region, other$status), collapse = ", "),
        "\n"
      )
    }
    if (once) {
      return(invisible(df))
    }
    if (nrow(running) == 0 && nrow(other) == 0) {
      cat("\nAll regions finished.\n")
      return(invisible(df))
    }
    Sys.sleep(refresh)
  }
}

# MCMC Diagnostics ----
computeInefficency <- function(folder, region) {
  path <- file.path(folder, paste0("FACTORS", region, ".xlsx"))
  sheets <- readxl::excel_sheets(path)

  lambdamm <- as.matrix(readxl::read_excel(path, sheet = "lambdamm"))
  phimm    <- as.matrix(readxl::read_excel(path, sheet = "phimm"))
  sigmm    <- as.matrix(readxl::read_excel(path, sheet = "sigmm"))

  draws <- list(lambda = lambdamm, phi = phimm, sigma = sigmm)

  if ("sig2fmm" %in% sheets) {
    sig2fmm <- as.matrix(readxl::read_excel(path, sheet = "sig2fmm"))
    draws$sig2f <- sig2fmm
  }
  if ("taumm" %in% sheets) {
    taumm <- as.matrix(readxl::read_excel(path, sheet = "taumm"))
    draws$tau <- taumm
  }

  M <- nrow(lambdamm)

  computeIF <- function(mat) {
    ess <- coda::effectiveSize(coda::mcmc(mat))
    M / ess
  }

  computeGeweke <- function(mat) {
    gd <- coda::geweke.diag(coda::mcmc(mat))
    gd$z
  }

  ifList <- lapply(draws, computeIF)
  gewekeList <- lapply(draws, computeGeweke)

  byParam <- purrr::imap_dfr(ifList, function(ifs, group) {
    tibble::tibble(group = group, parameter = names(ifs), IF = as.numeric(ifs))
  }) %>%
    dplyr::filter(is.finite(IF)) %>%
    dplyr::arrange(dplyr::desc(IF))

  byGroup <- byParam %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(mean_IF = mean(IF), max_IF = max(IF), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_IF))

  gewekeByParam <- purrr::imap_dfr(gewekeList, function(zs, group) {
    tibble::tibble(group = group, parameter = names(zs), geweke_z = as.numeric(zs))
  }) %>%
    dplyr::filter(is.finite(geweke_z))

  list(byParam = byParam, byGroup = byGroup, gewekeByParam = gewekeByParam)
}

computeAllDiagnostics <- function(folder) {
  files <- list.files(folder, pattern = "^FACTORS.*\\.xlsx$", full.names = TRUE)
  regions <- stringr::str_extract(basename(files), "(?<=FACTORS)[A-Z]+")
  purrr::map_dfr(regions, function(reg) {
    diag <- computeInefficency(folder, reg)
    ifDf <- diag$byParam %>% dplyr::mutate(region = reg, .before = 1)
    gwDf <- diag$gewekeByParam %>% dplyr::mutate(region = reg, .before = 1)
    dplyr::left_join(ifDf, gwDf, by = c("region", "group", "parameter"))
  })
}

# Forecast Evaluation ----

scoreEvalResults <- function(
  evalFolder, cutoffGrid, grpData, ivboYtdData, weightGrpFile,
  codes, cfg, districtMapPath = "Technical/district_map.csv"
) {
  modelForecasts <- collectEvalForecasts(cutoffGrid, evalFolder)
  actualGrp <- grpData %>%
    group_by(Abbrev) %>%
    arrange(date) %>%
    mutate(grpYoY = log(value / lag(value)) * 100) %>%
    filter(!is.na(grpYoY))
  benchmarks <- list(
    rw = computeRandomWalkBenchmark(grpData, cutoffGrid),
    ivbo = computeIvboBenchmark(ivboYtdData, cutoffGrid)
  )
  wMat <- readRDS(weightGrpFile)
  wMonths <- seq.Date(cfg$estim_start, cfg$sample_end, by = "month")
  decIdx <- which(month(wMonths) == 12)
  grpWeights <- tibble(
    region = rep(codes$Abbrev, length(decIdx)),
    year = rep(year(wMonths[decIdx]), each = nrow(codes)),
    weight = as.vector(wMat[, decIdx])
  )
  districtMap <- read_csv(districtMapPath, show_col_types = FALSE)
  evaluateForecasts(modelForecasts, actualGrp, benchmarks, grpWeights, districtMap)
}

# Build grid of cutoff dates and evaluation map for forecast assessment
defineCutoffGrid <- function(
  horizons,
  minGrpObs,
  grpReleaseLag,
  estimStart,
  sampleEnd
) {
  cutoffMonths <- sort(unique((12 + horizons) %% 12))
  cutoffMonths[cutoffMonths == 0] <- 12
  cutoffMonths <- sort(unique(cutoffMonths))
  firstGrpYear <- year(estimStart)
  minGrpThrough <- firstGrpYear + minGrpObs - 1
  earliestRelease <- as.Date(sprintf("%d-12-01", minGrpThrough)) +
    months(grpReleaseLag)
  lastCutoffMonth <- max(cutoffMonths[cutoffMonths <= month(sampleEnd)])
  lastCutoff <- as.Date(sprintf("%d-%02d-01", year(sampleEnd), lastCutoffMonth))
  allDates <- seq(
    as.Date(sprintf("%d-01-01", year(earliestRelease))),
    lastCutoff,
    by = "month"
  )
  cutoffDates <- allDates[
    month(allDates) %in% cutoffMonths & allDates >= earliestRelease
  ]
  grpYearAtCutoff <- function(d) {
    ref <- d - months(grpReleaseLag)
    year(ref) - as.integer(month(ref) < 12)
  }
  cutoffs <- tibble(
    cutoff_id = seq_along(cutoffDates),
    cutoff_date = cutoffDates,
    grp_available_through = sapply(cutoffDates, grpYearAtCutoff)
  )
  evalMap <- tidyr::crossing(cutoffs, tibble(horizon = horizons)) %>%
    mutate(
      cMonth = month(cutoff_date),
      target_year = as.integer(year(cutoff_date) - (horizon - cMonth + 12) / 12)
    ) %>%
    filter(
      (horizon - cMonth) %% 12 == 0,
      target_year > grp_available_through
    ) %>%
    select(cutoff_id, target_year, horizon)
  list(cutoffs = cutoffs, evalMap = evalMap)
}

# Truncate source data to a cutoff date, masking future GRP observations
truncateSourceData <- function(sourceDataRaw, cutoffDate, grpAvailableYear) {
  dates <- as.Date(sourceDataRaw[4:nrow(sourceDataRaw), 1][[1]])
  dat <- sourceDataRaw[c(1:3, 3 + which(dates <= cutoffDate)), ]
  aggs <- as.numeric(dat[2, -1])
  annualCols <- which(aggs == 1) + 1
  if (length(annualCols) > 0) {
    dataRows <- 4:nrow(dat)
    dateDat <- as.Date(dat[dataRows, 1][[1]])
    decRows <- dataRows[month(dateDat) == 12 & year(dateDat) > grpAvailableYear]
    if (length(decRows) > 0) dat[decRows, annualCols] <- NA
  }
  dat
}

# Collect model forecasts from per-cutoff estimation output folders
collectEvalForecasts <- function(cutoffGrid, evalFolder) {
  purrr::map_dfr(unique(cutoffGrid$evalMap$cutoff_id), function(cid) {
    ct <- cutoffGrid$cutoffs %>% filter(cutoff_id == cid)
    targets <- cutoffGrid$evalMap %>% filter(cutoff_id == cid)
    folder <- file.path(
      evalFolder,
      sprintf("CUT_%s", format(ct$cutoff_date, "%Y_%m"))
    )
    files <- list.files(
      folder,
      pattern = "^FACTORS.*\\.xlsx$",
      full.names = TRUE
    )
    purrr::map_dfr(files, function(f) {
      reg <- stringr::str_extract(basename(f), "(?<=FACTORS)[A-Z]+")
      fitted <- read_excel(f, sheet = "fitted")
      grpCol <- names(fitted)[grepl("^grp\\.", names(fitted))][1]
      if (is.na(grpCol)) {
        return(NULL)
      }
      fittedDates <- as.Date(fitted$date)
      grpDraws <- read_excel(f, sheet = "grpmm")
      purrr::map_dfr(seq_len(nrow(targets)), function(j) {
        targetDec <- as.Date(sprintf("%d-12-01", targets$target_year[j]))
        idx <- match(targetDec, fittedDates)
        if (is.na(idx)) {
          return(NULL)
        }
        tibble(
          region = reg,
          target_year = targets$target_year[j],
          horizon = targets$horizon[j],
          cutoff_id = cid,
          model_forecast = as.numeric(fitted[idx, grpCol]),
          draws = list(as.numeric(grpDraws[[as.character(targetDec)]]))
        )
      })
    })
  })
}

# Collect forecasts for a single cutoff (used as a targets branch)
collectCutoffForecasts <- function(cutoffId, cutoffGrid, evalFolder) {
  ct <- cutoffGrid$cutoffs %>% filter(cutoff_id == cutoffId)
  targets <- cutoffGrid$evalMap %>% filter(cutoff_id == cutoffId)
  folder <- file.path(
    evalFolder,
    sprintf("CUT_%s", format(ct$cutoff_date, "%Y_%m"))
  )
  files <- list.files(folder, pattern = "^FACTORS.*\\.xlsx$", full.names = TRUE)
  purrr::map_dfr(files, function(f) {
    reg <- stringr::str_extract(basename(f), "(?<=FACTORS)[A-Z]+")
    fitted <- read_excel(f, sheet = "fitted")
    grpCol <- names(fitted)[grepl("^grp\\.", names(fitted))][1]
    if (is.na(grpCol)) return(NULL)
    fittedDates <- as.Date(fitted$date)
    grpDraws <- read_excel(f, sheet = "grpmm")
    purrr::map_dfr(seq_len(nrow(targets)), function(j) {
      targetDec <- as.Date(sprintf("%d-12-01", targets$target_year[j]))
      idx <- match(targetDec, fittedDates)
      if (is.na(idx)) return(NULL)
      tibble(
        region = reg,
        target_year = targets$target_year[j],
        horizon = targets$horizon[j],
        cutoff_id = cutoffId,
        model_forecast = as.numeric(fitted[idx, grpCol]),
        draws = list(as.numeric(grpDraws[[as.character(targetDec)]]))
      )
    })
  })
}

# Score pre-collected forecasts against actuals and benchmarks
scoreFromForecasts <- function(
  modelForecasts, cutoffGrid, grpData, ivboYtdData, weightGrpFile,
  codes, cfg, districtMapPath = "Technical/district_map.csv"
) {
  actualGrp <- grpData %>%
    group_by(Abbrev) %>%
    arrange(date) %>%
    mutate(grpYoY = log(value / lag(value)) * 100) %>%
    filter(!is.na(grpYoY))
  benchmarks <- list(
    rw = computeRandomWalkBenchmark(grpData, cutoffGrid),
    ivbo = computeIvboBenchmark(ivboYtdData, cutoffGrid)
  )
  wMat <- readRDS(weightGrpFile)
  wMonths <- seq.Date(cfg$estim_start, cfg$sample_end, by = "month")
  decIdx <- which(month(wMonths) == 12)
  grpWeights <- tibble(
    region = rep(codes$Abbrev, length(decIdx)),
    year = rep(year(wMonths[decIdx]), each = nrow(codes)),
    weight = as.vector(wMat[, decIdx])
  )
  districtMap <- read_csv(districtMapPath, show_col_types = FALSE)
  evaluateForecasts(modelForecasts, actualGrp, benchmarks, grpWeights, districtMap)
}

# Compute random-walk benchmark forecasts (last observed GRP growth)
computeRandomWalkBenchmark <- function(grpData, cutoffGrid) {
  grpGrowth <- grpData %>%
    group_by(Abbrev) %>%
    arrange(date) %>%
    mutate(grpYoY = log(value / lag(value)) * 100, grpYear = year(date)) %>%
    filter(!is.na(grpYoY)) %>%
    ungroup()
  cutoffGrid$evalMap %>%
    left_join(
      cutoffGrid$cutoffs %>% select(cutoff_id, grp_available_through),
      by = "cutoff_id"
    ) %>%
    inner_join(
      grpGrowth %>% transmute(region = Abbrev, grpYear, bm_forecast = grpYoY),
      by = c("grp_available_through" = "grpYear")
    ) %>%
    select(region, target_year, horizon, cutoff_id, bm_forecast)
}

computeIvboBenchmark <- function(ivboYtdData, cutoffGrid) {
  cutoffGrid$evalMap %>%
    left_join(
      cutoffGrid$cutoffs %>% select(cutoff_id, cutoff_date),
      by = "cutoff_id"
    ) %>%
    mutate(
      ivboAvail = ifelse(
        month(cutoff_date) >= 3,
        year(cutoff_date) - 1,
        year(cutoff_date) - 2
      )
    ) %>%
    filter(target_year <= ivboAvail) %>%
    inner_join(
      ivboYtdData %>%
        transmute(region = Abbrev, year, bm_forecast = ivboYtdGrowth),
      by = c("target_year" = "year")
    ) %>%
    select(region, target_year, horizon, cutoff_id, bm_forecast)
}

dmTest <- function(e1, e2, h = 1) {
  if (length(e1) < 3) return(list(statistic = NA_real_, p.value = NA_real_))
  res <- tryCatch(
    forecast::dm.test(e1, e2, alternative = "less", h = h, power = 1),
    error = function(e) list(statistic = NA_real_, p.value = NA_real_)
  )
  list(statistic = as.numeric(res$statistic), p.value = res$p.value)
}

# Score model forecasts against actuals and benchmarks, aggregate by multiple dimensions
evaluateForecasts <- function(
  modelForecasts,
  actualGrp,
  benchmarks,
  grpWeights = NULL,
  districtMap = NULL
) {
  bmNames <- names(benchmarks)
  combined <- modelForecasts %>%
    inner_join(
      actualGrp %>%
        transmute(region = Abbrev, target_year = year(date), actual = grpYoY),
      by = c("region", "target_year")
    ) %>%
    mutate(error_model = model_forecast - actual)
  for (bm in bmNames) {
    fcCol <- paste0("forecast_", bm)
    errCol <- paste0("error_", bm)
    combined <- combined %>%
      left_join(
        benchmarks[[bm]] %>%
          transmute(region, target_year, horizon, cutoff_id,
                    !!fcCol := bm_forecast),
        by = c("region", "target_year", "horizon", "cutoff_id")
      ) %>%
      mutate(!!errCol := .data[[fcCol]] - actual)
  }
  combined <- combined %>%
    mutate(crps = purrr::map2_dbl(actual, draws, scoringRules::crps_sample))
  aggregate <- function(data, groupVars, w = NULL, dm = FALSE) {
    purrr::map_dfr(bmNames, function(bm) {
      errCol <- paste0("error_", bm)
      df <- data %>% filter(!is.na(.data[[errCol]]))
      if (is.null(w)) {
        res <- df %>%
          group_by(across(all_of(groupVars))) %>%
          summarise(
            benchmark = bm, n = n(),
            rmse_model = sqrt(mean(error_model^2)),
            mae_model = mean(abs(error_model)),
            crps = mean(crps),
            rmse_bm = sqrt(mean(.data[[errCol]]^2)),
            mae_bm = mean(abs(.data[[errCol]])),
            rel_rmse = sqrt(mean(error_model^2)) / sqrt(mean(.data[[errCol]]^2)),
            rel_mae = mean(abs(error_model)) / mean(abs(.data[[errCol]])),
            .groups = "drop"
          )
      } else {
        res <- df %>%
          group_by(across(all_of(groupVars))) %>%
          summarise(
            benchmark = bm, n = n(),
            wmae_model = weighted.mean(abs(error_model), .data[[w]]),
            wcrps = weighted.mean(crps, .data[[w]]),
            wmae_bm = weighted.mean(abs(.data[[errCol]]), .data[[w]]),
            rel_wmae = weighted.mean(abs(error_model), .data[[w]]) /
              weighted.mean(abs(.data[[errCol]]), .data[[w]]),
            .groups = "drop"
          )
      }
      if (dm) {
        dmRes <- df %>%
          group_by(across(all_of(groupVars))) %>%
          summarise(
            benchmark = bm,
            {
              r <- dmTest(error_model, .data[[errCol]], h = 1)
              tibble(dm_stat = r$statistic, dm_pvalue = r$p.value)
            },
            .groups = "drop"
          )
        res <- res %>% left_join(dmRes, by = c(groupVars, "benchmark"))
      }
      res
    }) %>%
      arrange(across(all_of(groupVars)), benchmark)
  }
  byHorizon <- aggregate(combined, "horizon")
  byRegion <- aggregate(combined, c("region", "horizon"), dm = TRUE)
  byYear <- aggregate(combined, c("target_year", "horizon"))
  result <- list(
    byHorizon = byHorizon, byRegion = byRegion,
    byYear = byYear, details = combined
  )
  if (!is.null(grpWeights)) {
    dw <- combined %>%
      inner_join(grpWeights %>% rename(target_year = year),
                 by = c("region", "target_year"))
    result$byHorizonW <- aggregate(dw, "horizon", w = "weight")
    byYearW <- aggregate(dw, c("target_year", "horizon"), w = "weight")
    result$byYear <- result$byYear %>%
      left_join(byYearW, by = c("target_year", "horizon", "benchmark"))
  }
  if (!is.null(districtMap)) {
    combinedDist <- combined %>%
      inner_join(districtMap %>% select(Abbrev, FD),
                 by = c("region" = "Abbrev")) %>%
      rename(district = FD)
    result$byDistrict <- aggregate(combinedDist, c("district", "horizon"))
    distDM <- purrr::map_dfr(bmNames, function(bm) {
      errCol <- paste0("error_", bm)
      combinedDist %>%
        filter(!is.na(.data[[errCol]])) %>%
        group_by(district, horizon, cutoff_id) %>%
        summarise(
          avg_err_model = mean(abs(error_model)),
          avg_err_bm = mean(abs(.data[[errCol]])),
          .groups = "drop"
        ) %>%
        group_by(district, horizon) %>%
        summarise(
          benchmark = bm,
          {
            r <- dmTest(avg_err_model, avg_err_bm, h = 1)
            tibble(dm_stat = r$statistic, dm_pvalue = r$p.value)
          },
          .groups = "drop"
        )
    })
    result$byDistrict <- result$byDistrict %>%
      left_join(distDM, by = c("district", "horizon", "benchmark"))
    if (!is.null(grpWeights)) {
      dwDist <- dw %>%
        inner_join(districtMap %>% select(Abbrev, FD),
                   by = c("region" = "Abbrev")) %>%
        rename(district = FD)
      wDist <- aggregate(dwDist, c("district", "horizon"), w = "weight")
      wDistDM <- purrr::map_dfr(bmNames, function(bm) {
        errCol <- paste0("error_", bm)
        dwDist %>%
          filter(!is.na(.data[[errCol]])) %>%
          group_by(district, horizon, cutoff_id) %>%
          summarise(
            wavg_err_model = weighted.mean(abs(error_model), weight),
            wavg_err_bm = weighted.mean(abs(.data[[errCol]]), weight),
            .groups = "drop"
          ) %>%
          group_by(district, horizon) %>%
          summarise(
            benchmark = bm,
            {
              r <- dmTest(wavg_err_model, wavg_err_bm, h = 1)
              tibble(wdm_stat = r$statistic, wdm_pvalue = r$p.value)
            },
            .groups = "drop"
          )
      })
      wDist <- wDist %>%
        left_join(wDistDM, by = c("district", "horizon", "benchmark"))
      result$byDistrict <- result$byDistrict %>%
        left_join(wDist, by = c("district", "horizon", "benchmark"))
    }
  }
  result
}
