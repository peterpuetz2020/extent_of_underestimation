# setup ----
# in case not installed, install pacman package by using the command
# install.packages("pacman", repos="https://cran.rstudio.com/")
# clean up environment
rm(list = ls())

# (install and) load packages
pacman::p_load(
  here,
  tidyverse,
  MASS,
  lubridate,
  MMWRweek,
  cowplot,
  segmented,
  fANCOVA,
  haven,
  padr,
  zoo
)

# set paths where data and results are located / stored
read_data_here <-
  here(here(), "In")
results_here <-
  here(here(), "Out")

# read in data (see Table 1 in the corresponding paper
# for a description of the variables)
df <- read_csv(here(read_data_here, "uaf_clean.csv"))

df <- df %>%
  # create date variable for the calendar week of the year
  mutate(KW = MMWRweek2Date(
    MMWRyear = Jahr,
    MMWRweek = Woche,
    MMWRday = 2
  ))

# get number of observations for variables that are smoothed
no_obs <- df %>%
  dplyr::select(GW_VPR, GW_SR, SC2_VL_WW) %>%
  summarise_all(~ sum(!is.na(.)))

# expand data set (for graphical purposes only)
df <- df %>%
  pad(by = "KW",
      interval = "day") %>%
  fill(Jahr, Woche, .direction = "down") %>%
  # also for graphical purposes: fill in incidence data
  mutate(GNS = na.spline(GNS))

# produce loess estimates
df <- df %>%
  mutate(
    obs = row_number(),
    GW_VPR_LOESS_037 = predict(
      loess.as(
        obs[!is.na(GW_VPR)],
        GW_VPR[!is.na(GW_VPR)],
        degree = 1,
        user.span = 8 / (no_obs %>% pull(GW_VPR))
      ),
      newdata = data.frame(x = obs)
    ),
    GW_SR_LOESS_091 = predict(
      loess.as(
        obs[!is.na(GW_SR)],
        GW_SR[!is.na(GW_SR)],
        degree = 1,
        user.span = 8 / (no_obs %>% pull(GW_SR))
      ),
      newdata = data.frame(x = obs)
    ),
    SC2_VL_WW_LOESS_09 = predict(
      loess.as(
        obs[!is.na(SC2_VL_WW)],
        SC2_VL_WW[!is.na(SC2_VL_WW)],
        degree = 1,
        user.span = 8 / (no_obs %>% pull(SC2_VL_WW))
      ),
      newdata = data.frame(x = obs)
    )
  )

# calculate underestimation factors based on smoothed values
df <- df %>%
  mutate(UEF_GW_VPR = GW_VPR_LOESS_037 / GNS,
         UEF_GW_SR = GW_SR_LOESS_091 / GNS)

# number of observations for variables that are smoothed
no_obs <- df %>%
  dplyr::select(UEF_GW_VPR, UEF_GW_SR) %>%
  summarise_all(~ sum(!is.na(.)))

# produce loess estimates for underestimation factors
df <- df %>%
  mutate(
    UEF_GW_VPR_LOESS_056 = predict(
      loess.as(
        obs[!is.na(UEF_GW_VPR)],
        UEF_GW_VPR[!is.na(UEF_GW_VPR)],
        degree = 1,
        user.span = 8 / (no_obs %>% pull(UEF_GW_VPR))
      ),
      newdata = data.frame(x = obs)
    ),
    UEF_GW_SR_LOESS_136 = predict(
      loess.as(
        obs[!is.na(UEF_GW_SR)],
        UEF_GW_SR[!is.na(UEF_GW_SR)],
        degree = 1,
        user.span = 8 / (no_obs %>% pull(UEF_GW_SR))
      ),
      newdata = data.frame(x = obs)
    )
  )

# plot A ------------------------------------------------------------------

# time frames for dominance of variants
# Alpha			 09/2021 - 23/2021
# Delta			31/2021 - 51/2021
# Omicron BA.1		52/2021 - 8/2022
# Omicron BA.2		9/2022 - 22/2022
# Omicron BA.5		23/2022 - 5/2023
# Recombinant Lineages		8/2023 - 47/2023
# Omicron BA.2		48/2023 - 4/2024 (current)

# for plotting purposes transfer the time frames above to a data frame
variant_data = tibble(
  x1 = c(
    df %>% filter(Woche == 9, Jahr == 2021) %>%  slice(1) %>% pull(KW),
    df %>% filter(Woche == 31, Jahr == 2021) %>%  slice(1)  %>% pull(KW),
    df %>% filter(Woche == 52, Jahr == 2021) %>%  slice(1) %>% pull(KW),
    df %>% filter(Woche == 9, Jahr == 2022) %>%  slice(1) %>% pull(KW),
    df %>% filter(Woche == 23, Jahr == 2022) %>%  slice(1) %>% pull(KW),
    df %>% filter(Woche == 8, Jahr == 2023) %>%  slice(1) %>% pull(KW),
    df %>% filter(Woche == 48, Jahr == 2023) %>%  slice(1)  %>% pull(KW)
  ),
  x2 = c(
    df %>% filter(Woche == 23, Jahr == 2021) %>%  slice(1)  %>% pull(KW),
    df %>% filter(Woche == 51, Jahr == 2021) %>%  slice(1) %>% pull(KW),
    df %>% filter(Woche == 8, Jahr == 2022) %>%  slice(1) %>% pull(KW),
    df %>% filter(Woche == 22, Jahr == 2022) %>%  slice(1) %>% pull(KW),
    df %>% filter(Woche == 5, Jahr == 2023) %>%  slice(1) %>% pull(KW),
    df %>% filter(Woche == 47, Jahr == 2023) %>%  slice(1) %>% pull(KW),
    df %>% filter(Woche == 4, Jahr == 2024) %>%  slice(1) %>% pull(KW)
  ),
  y1 = rep(2600, 7),
  y2 = rep(2900, 7),
  label = c(
    "Alpha",
    "Delta",
    "Omicron\nBA.1",
    "Omicron\nBA.2",
    "Omicron\nBA.5",
    "Recombinant Lineages",
    "Omicron\nBA.2"
  )
)

# get variables that determine the two axes for plot A
x = pull(df, GW_SR_LOESS_091)
z = pull(df, SC2_VL_WW_LOESS_09)

# compute transformation of axes (for plotting purposes)
ylim.prim = c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))
ylim.sec = c(min(z, na.rm = TRUE), max(z, na.rm = TRUE))
beta <- diff(ylim.sec) / diff(ylim.prim)
alpha = ylim.sec[1] - beta * ylim.prim[1]

# reshape data
plot_data <-
  df %>%
  dplyr::select(-contains("UEF"),-GW_SR,-SC2_VL_WW,-GW_VPR,-obs) %>%
  gather(criterion, value,-Jahr,-Woche,-KW)

# generate plot A
plot_1 <- df %>%
  ggplot() +
  scale_y_continuous(
    breaks = c(0, 500, 1000, 1500, 2000, 2500),
    limits = c(0, 2900),
    expand = c(0.02, 0),
    sec.axis = sec_axis(~ . * beta + alpha,
                        
                        name = "SARS-CoV-2 viral load in wastewater\nin 1,000 gene copies per liter")
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.line = element_line(),
    legend.title = element_blank(),
    text = element_text(size = 16),
    axis.title.y.right = element_text(vjust = +2),
    panel.background = element_rect(fill = 'white'),
    legend.key = element_blank()
  ) +
  
  scale_x_date(
    breaks = plot_data$KW[seq(1, length(unique(plot_data$KW)), 70)],
    labels =  paste0(plot_data$Woche[seq(1, length(unique(plot_data$KW)), 70)],
                     "\n",
                     plot_data$Jahr[seq(1, length(unique(plot_data$KW)), 70)]),
    expand = c(0.01, 0.01)
  ) +
  labs(x = "Calendar week and year", y = "COVID-19 cases\n per 100,000 inhabitants (adults)") +
  geom_rect(
    data = variant_data,
    aes(
      xmin = x1,
      xmax = x2,
      ymin = y1,
      ymax = y2
    ),
    color = "black",
    alpha = 0.5
  ) +
  geom_text(
    data = variant_data,
    aes(
      x = x1 + (x2 - x1) / 2,
      y = y1 + (y2 - y1) / 2,
      label = label
    ),
    size = 4.5,
    lineheight = .8
  )  +
  geom_line(aes(KW, GNS, color = "GNS-I (non-sm)"),  linewidth = 1.0) +
  geom_line(aes(KW, GW_VPR_LOESS_037, color = "GW-VPR-I (sm)"),  linewidth = 1.0) +
  geom_line(aes(KW, GW_SR_LOESS_091, color = "GW-SR-I (sm)"),  linewidth = 1.0) +
  geom_line(aes(KW, (SC2_VL_WW_LOESS_09 - alpha) / beta, color = "SC2-VL-WW (sm)"),
            linewidth = 1.0) +
  scale_color_manual(name = "",
                     values = c("black", "#E69F00", "#56B4E9", "#999999")) +
  guides(color = guide_legend(title = "", override.aes =
                                list(lwd = rep(0.7, 4))))

plot_1

# change point regression -------------------------------------------------

# transform underestimation factos
df <- df %>%
  mutate(across(contains("UEF"), ~ log10(.)))

# reshape data
df_gathered <- df %>%
  gather(criterion, UEF,-Jahr,-Woche,-obs)

# fit change point regression for 3 breakpoints with smoothed data
fit_segmented <- segmented(lm(UEF ~ obs, data = df_gathered %>%
                                filter(
                                  criterion %in% c("UEF_GW_VPR_LOESS_056",
                                                   "UEF_GW_SR_LOESS_136")
                                )),
                           seg.Z = ~ obs,
                           npsi = 3)
summary(fit_segmented)

# get change points (round them up)
infl_points <- ceiling(fit_segmented$psi[, 2])

df_exp <- df %>%
  # add fitted segmented regression line
  add_column(fitted = predict(fit_segmented,
                              newdata = data.frame(obs = df$obs))) %>%
  # add variable for start of new phases (after breakpoints)
  mutate(infl = (obs %in% infl_points) * 1)

# show breakpoints
df_exp %>%
  filter(obs %in% infl_points) %>%
  dplyr::select(Jahr, Woche, KW)

# plot B ------------------------------------------------------------------

# reshape and transform data
plot_data <-
  df_exp %>%
  dplyr::select(contains("UEF"), Jahr, Woche, KW, fitted, infl) %>%
  gather(criterion, UEF,-Jahr,-Woche,-KW,-contains("inf")) %>%
  mutate(across(c(contains("UEF")),
                ~ 10 ^ (.)))

# for plotting purposes transfer the time frames for the phases
# determined by the breakpoints above to a data frame
phase_data = tibble(
  x1 = c(
    df_exp %>% filter(KW == min(KW)) %>% pull(KW),
    df_exp %>% filter(infl == 1) %>% slice(1) %>% pull(KW) -
      0,
    df_exp %>% filter(infl == 1) %>% slice(2) %>% pull(KW) -
      0,
    df_exp %>% filter(infl == 1) %>% slice(3) %>% pull(KW) -
      0
  ),
  x2 = c(
    df_exp %>% filter(infl == 1) %>% slice(1) %>% pull(KW) - 0,
    df_exp %>% filter(infl == 1) %>% slice(2) %>% pull(KW) -
      0,
    df_exp %>% filter(infl == 1) %>% slice(3) %>% pull(KW) -
      0,
    df_exp %>% filter(KW == max(KW)) %>% pull(KW)
  ),
  y1 = rep(100, 4),
  y2 = rep(300, 4),
  label = c("phase 1",
            "phase 2",
            "phase 3",
            "phase 4")
)

# generate plot B
plot_2 <- plot_data %>%
  ggplot() +
  geom_rect(
    data = phase_data,
    aes(
      xmin = x1,
      xmax = x2,
      ymin = y1,
      ymax = y2
    ),
    color = "white",
    fill = "white",
    alpha = 0.5
  ) +
  geom_text(data = phase_data,
            aes(
              x = x1 + (x2 - x1) / 2,
              y = y1 + (y2 - y1) / 2,
              label = label
            ),
            size = 5) +
  geom_vline(
    xintercept = plot_data %>% filter(infl == 1) %>% pull(KW),
    color = "grey",
    linewidth = 1.2,
    show.legend = TRUE
  ) +
  geom_point(aes(KW, UEF, color = criterion),
             data = plot_data %>% filter(criterion %in%
                                           c("UEF_GW_VPR",
                                             "UEF_GW_SR"))) +
  geom_line(
    aes(KW, UEF, color = criterion, linetype = criterion),
    linewidth = 1.0,
    data = plot_data %>% filter(
      criterion %in%
        c("UEF_GW_VPR_LOESS_056",
          "UEF_GW_SR_LOESS_136",
          "fitted")
    )
  ) +
  scale_x_date(
    breaks = plot_data$KW[seq(1, length(unique(plot_data$KW)), 70)],
    labels =  paste0(plot_data$Woche[seq(1, length(unique(plot_data$KW)), 70)],
                     "\n",
                     plot_data$Jahr[seq(1, length(unique(plot_data$KW)), 70)]),
    expand = c(0.01, 0.01)
  ) +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    text = element_text(size = 16),
    axis.line = element_line(),
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.key = element_blank(),
    legend.position = "bottom"
  ) +
  scale_color_manual(
    name = "",
    breaks =  c(
      "UEF_GW_SR",
      "UEF_GW_SR_LOESS_136",
      "UEF_GW_VPR",
      "UEF_GW_VPR_LOESS_056",
      "fitted"
    ),
    labels =  c(
      "UEF-GW-SR-I",
      "UEF-GW-SR-I (sm)",
      "UEF-GW-VPR-I",
      "UEF-GW-VPR-I (sm)",
      "Fitted segmented regression line"
    ),
    values = c("#E69F00", "#E69F00", "#56B4E9", "#56B4E9",  "black")
  ) +
  scale_linetype_manual(
    values = c(
      "fitted" = "solid",
      "UEF_GW_VPR_LOESS_056" = "solid",
      "UEF_GW_SR_LOESS_136" = "solid"
    ),
    guide = "none"
  ) +
  guides(color = guide_legend(
    title = "",
    override.aes =
      list(
        linetype = c(0, 1, 0, 1, 1),
        shape = c(16, NA, 16, NA, NA),
        lwd = rep(1, 5)
      )
  )) +
  
  scale_shape(guide = "none") +
  labs(x = "Calendar week and year", y = "Underestimation factor (log10-scale)") +
  scale_y_continuous(
    trans = 'log10',
    expand = c(0.01, 0),
    breaks = c(.5, 1, 1.5, 2.5, 5 , 10, 25, 50, 100, 200),
    limits = c(0.1, 300)
  )
plot_2

# put both plots together
plot_grid(
  plot_1,
  plot_2,
  ncol = 1,
  labels = "AUTO",
  align = "v",
  label_x = 0.1,
  label_y = .95,
  label_size = 17
)

# save plot
ggsave(
  here(results_here, "plots.png"),
  width = 45,
  height = 30,
  units = "cm"
)