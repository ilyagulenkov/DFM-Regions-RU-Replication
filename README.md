# DFM-Regions-RU

Код и данные для репликации статьи "Оперативная оценка региональной экономической активности в России с использованием динамических факторных моделей"

## Структура

```
Data/РегДанные.xlsx      Исходные данные
R/                       Код
  _targets.R             Основной pipeline
  functions.R            Вспомогательные функции на R
  functions.cpp          Вспомогательные функции на Rcpp
  precision_sampler.cpp  Векторизованный фильтр Калмана Antolín-Díaz et al (2024)  
  plot_eval_timeline.R   Визуализация вневыборочного прогноза
_targets.yaml            Конфигурация pipeline
Paper/Submission-1/      Рукопись и приложения (R Markdown)
Technical/               Технические файлы
Literature/DFMRU.bib     Библиография
```

## Требования

- R (>= 4.2) и компилятор C++ для Rcpp (Rtools в Windows).
- Пакеты R: targets, crew, tidyverse, readxl, writexl, tsibble, foreach,
  seastests, deseats, seasonal, tempdisagg, zoo, xts, glue, flextable, huxtable,
  imputeTS, Rcpp, RcppArmadillo, scoringRules, Matrix, coda. Для сборки рукописи
  дополнительно используются rmarkdown, knitr, showtext и sysfonts.

## Воспроизведение результатов

1. Сохранить локальную копию репозитория
2. В корневой директории репозитория исполнить команду:

   ```r
   targets::tar_make()
   ```

   В результате создаются файл `Data/monthly.xlsx`, матрицы весов и сводные статистики в формате
   `.RDS`, а также порегиональные результаты в каталоге `Results/`. Расчёт ресурсоёмкий: для каждого региона параллельно выполняется процедура
   MCMC (по умолчанию 16 рабочих процессов; их число задаётся в вызове
   `crew_controller_local()` в файле `R/_targets.R`). Ожидаемое время работы: несколько суток.

3. После завершения расчётов рукопись и приложения могут быть собраны с помощью команд:

   ```r
   rmarkdown::render("Paper/Submission-1/main.Rmd")
   rmarkdown::render("Paper/Submission-1/appendix.Rmd")
   ```
