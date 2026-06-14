library(tidyverse)
library(showtext)

font_add("Times", "C:/Windows/Fonts/times.ttf")
showtext_auto()
options(scipen = 9999)

PLOT_FONT_SIZE <- 18
BAR_THICKNESS  <- 2   # gray history-bar line width; lower this if rows look crowded
POINT_SIZE     <- 4   # target diamond size

theme_set(theme_minimal() +
            theme(text = element_text(size = PLOT_FONT_SIZE),
                  legend.text = element_text(size = PLOT_FONT_SIZE),
                  strip.placement = "outside",
                  strip.text.x = element_text(size = PLOT_FONT_SIZE, angle = 0),
                  strip.text.y.left = element_text(size = PLOT_FONT_SIZE, angle = 0),
                  strip.text.y.right = element_text(size = PLOT_FONT_SIZE, angle = 0),
                  legend.position = "bottom",
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  axis.ticks.y = element_line(),
                  axis.ticks.x = element_line(),
                  axis.line.x = element_line(),
                  axis.line.y = element_line(),
                  axis.text.y.right = element_text(margin = margin(l = 0.5, r = 0))))

cols      <- c("#595959", "#262626")
cols_1    <- c("#8cc5e3", "#3594cc")
targetDec <- as.Date("2020-12-01")
dataStart <- as.Date("2019-02-01")
xMin      <- as.Date("2019-02-01")
xMax      <- as.Date("2022-04-01")

# ypos: 7 top rows, then a wider gap, then 2 bottom rows.
# Edit this vector to add/remove rows or widen the gap.
horizons <- tibble(
  h      = c(6, 5, 4, 3, 2, 1, 0, -3, -6),
  cutoff = as.Date(c("2020-06-01", "2020-07-01", "2020-08-01", "2020-09-01",
                     "2020-10-01", "2020-11-01", "2020-12-01",
                     "2021-03-01", "2021-06-01")),
  ypos   = c(8, 7, 6, 5, 4, 3, 2, 0, -1),
  label  = case_when(h > 0 ~ sprintf("h = +%d", h),
                     h < 0 ~ sprintf("h = %d",  h),
                     TRUE  ~ "h = 0")
)

xBreaks <- seq(as.Date("2019-03-01"), as.Date("2022-03-01"), by = "3 months")
xLabels <- format(xBreaks, "%B\n%Y")

pEvalTimeline <- ggplot() +
  annotate("rect",
           xmin = as.Date("2020-01-01"), xmax = as.Date("2020-12-31"),
           ymin = min(horizons$ypos) - 0.45, ymax = max(horizons$ypos) + 0.45,
           fill = cols_1[1], alpha = 0.2) +
  annotate("segment",
           x = xMin + 30, xend = xMax - 30,
           y = 1, yend = 1,
           linewidth = 0.4, linetype = "dashed", color = "grey70") +
  geom_segment(data = horizons,
               aes(x = dataStart, xend = cutoff, y = ypos, yend = ypos),
               linewidth = BAR_THICKNESS, color = "grey78", lineend = "butt") +
  geom_segment(data = filter(horizons, h != 0),
               aes(x = cutoff, xend = targetDec + ifelse(h > 0, -15, 15),
                   y = ypos, yend = ypos),
               linewidth = 1, linetype = "dashed", color = cols[1],
               arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  geom_point(data = horizons,
             aes(x = targetDec, y = ypos),
             shape = 18, size = POINT_SIZE, color = cols_1[2]) +
  geom_segment(data = horizons,
               aes(x = cutoff, xend = cutoff, y = ypos - 0.28, yend = ypos + 0.28),
               linewidth = 1.2, color = cols[2]) +
  scale_x_date(
    limits = c(xMin, xMax),
    breaks = xBreaks,
    labels = xLabels,
    expand = expansion(mult = 0.01)) +
  scale_y_continuous(
    breaks = horizons$ypos,
    labels = horizons$label,
    expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = NULL) +
  theme(
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    axis.text.y = element_text(family = "Times", size = PLOT_FONT_SIZE * 0.9, hjust = 1),
    axis.text.x = element_text(family = "Times", size = PLOT_FONT_SIZE * 0.65),
    plot.margin = margin(6, 12, 6, 6))

pEvalTimeline
