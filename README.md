# DFM-Regions-RU -- Пакет для воспроизведения результатов

Код и данные для воспроизведения результатов динамической факторной модели
(DFM), применяемой для наукастинга валового регионального продукта (ВРП)
субъектов Российской Федерации.

## Структура

```
Data/РегДанные.xlsx      Исходная панель данных (региональные показатели Росстата)
R/                       Вычислительный конвейер (targets + Rcpp/RcppArmadillo)
  _targets.R             Определение конвейера
  functions.R            Подготовка данных, оценивание, расчёт мер качества
  functions.cpp          Процедуры на Rcpp
  precision_sampler.cpp  Сэмплер состояния на основе матрицы точности
  plot_eval_timeline.R   График динамики прогнозной точности
_targets.yaml            Конфигурация targets
Paper/Submission-1/      Рукопись и приложение (R Markdown)
Technical/               Стиль CSL, шаблоны Word, карта федеральных округов
Literature/DFMRU.bib     База библиографических ссылок
```

## Требования

- R (>= 4.2) и компилятор C++ для Rcpp (Rtools в Windows).
- Пакеты R: targets, crew, tidyverse, readxl, writexl, tsibble, foreach,
  seastests, deseats, seasonal, tempdisagg, zoo, xts, glue, flextable, huxtable,
  imputeTS, Rcpp, RcppArmadillo, scoringRules, Matrix, coda. Для сборки рукописи
  дополнительно используются rmarkdown, knitr, showtext и sysfonts.

## Воспроизведение результатов

1. Сохраните файл `Data/РегДанные.xlsx` на месте.
2. Из корня репозитория запустите вычислительный конвейер в R:

   ```r
   targets::tar_make()
   ```

   Расчёт ресурсоёмкий: для каждого региона параллельно выполняется процедура
   MCMC (по умолчанию 16 рабочих процессов; их число задаётся в вызове
   `crew_controller_local()` в файле `R/_targets.R`). В результате создаются
   файл `Data/monthly.xlsx`, матрицы весов и сводные статистики в формате
   `.RDS`, а также порегиональные результаты в каталоге `Results/`.

3. Соберите рукопись из каталога `Paper/Submission-1/`:

   ```r
   rmarkdown::render("Paper/Submission-1/main.Rmd")
   rmarkdown::render("Paper/Submission-1/appendix.Rmd")
   ```

## Примечания

- Рукопись загружает шрифт Times по пути, специфичному для Windows
  (`C:/Windows/Fonts/times.ttf`). В macOS и Linux измените вызов `font_add()`
  в файлах `main.Rmd` и `appendix.Rmd`, указав путь к локальному шрифту Times.
- Исходные данные получены из Росстата (Федеральная служба государственной
  статистики).
