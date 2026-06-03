##################################################################################################################
#### Methods for Measuring Neural Activity During Voluntary Wheel Running                                      ####
#### Companion analysis script                                                                                 ####
####                                                                                                           ####
#### To adapt this workflow to new data:                                                                       ####
####   1. Place .rds files from FiPhA in `dir_path` defined in line 60                                         ####
####   2. Edit `mouse_lookup` to match the series IDs and labels of new mice                                   ####
####   3. Run the full script top to bottom                                                                    ####
####                                                                                                           ####
#### Nomenclature                                                                                              ####
####  # Event = a single running event and the 15 seconds before and after                                     ####
####  # Series = a single 30 minute recording that contains multiple running events                            ####
####  # Phase = acquisition (1st week of wheel exposure) or maintenance (after one month) series               ####
####                                                                                                           ####
#### All outputs save to a single `Analysis/` folder created next to the .rds files.                           ####
###############################################################################################################

#### Libraries ####
library(data.table)
library(ggplot2)
library(scico)
library(svglite)
library(ggridges)
library(nlme)
library(emmeans)
library(Hmisc)
library(MASS)
library(multcompView)
library(afex)

#### Configuration:  Mouse ID mapping ####
# To adapt this workflow to new mice, edit this single table.
# `series_id` is the numeric prefix for FiPhA series names of the original manuscript, a naming convention of the lab (ex: for 1120_RVDG_Aq1: 1120 is the mouse ID, RVDG is the right ventral dentate gyrus, and Aq1 is the first aquisition session)
# `label`     is the display name used throughout figures and tables.
# `color`     is the plotting color for that mouse.
mouse_lookup <- data.table(
  series_id = c("1120",     "1123",     "1124"),
  label     = c("Mouse1",   "Mouse2",   "Mouse3"),
  color     = c("#0E6251",  "#7D8B3A",  "#D4B96F")
)

# Extract series_id from any Start_/Stop_ object name, return mouse label. 
# Note, Start and Stop were defined in FiPhA during the event identification step shown in Figure 4e (ex: our preferred naming event naming convention: "Start_1120_RVDG_Aq1" or "Stop_1120_RVDG_Aq1")
.mouse_from <- function(nm) {
  id <- sub("^(?:Start|Stop)_(\\d+).*$", "\\1", nm, perl = TRUE)
  mouse_lookup$label[match(id, mouse_lookup$series_id)]
}

# Extract phase code/label from object name
.phase_code <- function(nm) ifelse(grepl("_Ma", nm), "Ma", "Aq")
.phase_long <- function(code) ifelse(code == "Ma", "Maintenance", "Acquisition")

# Mouse color palette, named by label in line 38 (use with scale_color_manual)
mouse_pal <- setNames(mouse_lookup$color, mouse_lookup$label)
# Mouse color palette, named by series_id (for plots that use raw IDs)
mouse_pal_id <- setNames(mouse_lookup$color, mouse_lookup$series_id)

#### Paths ####
# Input: directory containing .rds files from FiPhA
dir_path <- "path/to/data/folder"  # update to your local directory containing .rds files

# Output: all tables, figures, and modeling outputs nest under this root
out_root <- file.path(dir_path, "Analysis")
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

.fig_dir <- function(name) {
  d <- file.path(out_root, name)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

#### Import:         Event-level datasets (from FiPhA .rds files) ####
rds_files <- dir(dir_path, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(rds_files) > 0)

.pieces <- list()
.k <- 0L
for (r in rds_files) {
  .rds_obj <- readRDS(r)
  for (f in names(.rds_obj$events)) {
    for (s in names(.rds_obj$events[[f]])) {
      for (e in seq_along(.rds_obj$events[[f]][[s]])) {
        .k <- .k + 1L
        .pieces[[.k]] <- data.table(
          `[rds]` = basename(r),
          Series  = s,
          Event   = e,
          .rds_obj$events[[f]][[s]][[e]]
        )
      }
    }
  }
}
phot_data <- rbindlist(.pieces, fill = TRUE, use.names = TRUE)
rm(.pieces, .rds_obj, .k)

#### Import:         Complete session-level dataset (for unbiased off-wheel kinematics) ####
CompleteFiPhaDataset <- readRDS(
  file.path(dir_path, "RawFiPhaData.rds")
)
cat("CompleteFiPhaDataset loaded:", nrow(CompleteFiPhaDataset), "rows,",
    uniqueN(CompleteFiPhaDataset$Series), "sessions\n")

#### Transform:      Rename columns ####
setnames(phot_data,
         old = c("(time)", "(event time)", "(interval)", "(interval time)"),
         new = c("Absolute_time", "Event_time", "Interval_phase", "Interval_time"))

#### Transform:      Standardize times ####
phot_data[, `:=`(
  Absolute_time = round(Absolute_time, 2),
  Event_time    = round(Event_time, 2),
  Interval_time = round(Interval_time, 2)
)]

#### Transform:      Make a unique dataframe for each Series ####
series_list <- split(phot_data, by = "Series")
for (nm in names(series_list)) {
  dt <- series_list[[nm]]
  assign(paste0("Start_", nm),
         dt[Interval_phase %in% c("Before", "Running")], envir = .GlobalEnv)
  assign(paste0("Stop_", nm),
         dt[Interval_phase %in% c("Running", "After")],  envir = .GlobalEnv)
}

#### Transform:      Prefix event names by series ####
for (nm in ls(pattern = "^(Start_|Stop_)")) {
  dt <- get(nm, inherits = TRUE); setDT(dt)
  if (!"Event_raw" %in% names(dt)) dt[, Event_raw := Event]
  suf <- sub("^.*_((?:Aq|Ma)\\d+)$", "\\1", nm, perl = TRUE)
  dt[, Event := paste0(suf, "_", as.character(Event))]
  assign(nm, dt, envir = .GlobalEnv)
}

#### Transform:      For Stop intervals, reverse Running times ####
stop_series_names <- ls(pattern = "^Stop_")
for (nm in stop_series_names) {
  dt <- get(nm); setDT(dt)
  dt[, Event_time := {
    out     <- Event_time
    run_idx <- which(Interval_phase == "Running")
    if (length(run_idx) > 0) {
      vals <- out[run_idx]
      if (any(!is.na(vals))) {
        max_t <- max(vals, na.rm = TRUE)
        out[run_idx] <- round(vals - max_t, 2)
      }
    }
    out
  }, by = "Event"]
  assign(nm, dt, envir = .GlobalEnv)
}

#### Transform:      Convert After Event_time to Interval_time ####
for (nm in stop_series_names) {
  dt <- get(nm); setDT(dt)
  dt[Interval_phase == "After", Event_time := Interval_time]
  assign(nm, dt, envir = .GlobalEnv)
}
##################################################################################################################
#### Supplemental Figure 3 - Wheel Running event lengths and ACh amplitudes ridgeplots ####
#### Transform:      Calculate per-event wheel running lengths and z-scores ####

.start_names <- ls(pattern = "^Start_\\d{4}_", envir = .GlobalEnv)
.stop_names  <- ls(pattern = "^Stop_\\d{4}_",  envir = .GlobalEnv)
stopifnot(length(.start_names) > 0)

start_ev <- rbindlist(lapply(.start_names, function(nm) {
  dt_start <- get(nm, envir = .GlobalEnv)
  setDT(dt_start)
  stopifnot(all(c("Event","Event_time","Interval_phase","Ratio") %in% names(dt_start)))
  
  nm_stop <- sub("^Start_", "Stop_", nm)
  dt_stop <- if (nm_stop %in% .stop_names) {
    tmp <- get(nm_stop, envir = .GlobalEnv); setDT(tmp); tmp
  } else NULL
  
  # baseline rows: pre-event (-15 to -10) from Start_, post-event (10 to 15) from Stop_
  base_before <- dt_start[Interval_phase == "Before" & Event_time >= -15 & Event_time <= -10,
                          .(Event, Ratio)]
  base_after <- if (!is.null(dt_stop)) {
    dt_stop[Interval_phase == "After" & is.finite(Event_time) &
              Event_time >= 10 & Event_time <= 15, .(Event, Ratio)]
  } else data.table(Event = integer(), Ratio = numeric())
  base_all <- rbind(base_before, base_after, use.names = TRUE, fill = TRUE)
  
  # per-event baseline statistics
  base <- base_all[, .(mu = mean(Ratio, na.rm = TRUE),
                       sd = sd(Ratio,   na.rm = TRUE),
                       n  = .N),
                   by = Event][is.finite(mu) & is.finite(sd) & sd > 0 & n > 1]
  if (!nrow(base)) return(NULL)
  
  # z-score Start_ rows using their baseline
  dt <- merge(dt_start, base[, .(Event, mu, sd)], by = "Event", all.x = TRUE)
  dt[, z := (Ratio - mu) / sd]
  
  run <- dt[Interval_phase == "Running" & is.finite(Event_time) &
              is.finite(Ratio) & is.finite(z)]
  if (!nrow(run)) return(NULL)
  
  ans <- run[, .(
    event_length_seconds = max(Event_time, na.rm = TRUE) - min(Event_time, na.rm = TRUE),
    median_running_z     = median(z, na.rm = TRUE),
    n_samples_running    = .N
  ), by = .(Event)]
  
  ans[, `:=`(
    series_name = nm,
    phase       = .phase_code(nm),
    mouse       = .mouse_from(nm)
  )][]
}), fill = TRUE)

start_ev <- start_ev[is.finite(event_length_seconds) & is.finite(median_running_z)]
start_ev[, phase_long := .phase_long(phase)]
start_ev[, `:=`(
  phase = factor(phase, levels = c("Aq","Ma")),
  mouse = factor(as.character(mouse), levels = mouse_lookup$label)
)]
setorder(start_ev, phase, mouse, Event)

.summarize_lengths_and_z <- function(x) {
  x[, .(
    n_events                = .N,
    min_event_length_s      = min(event_length_seconds,    na.rm = TRUE),
    median_event_length_s   = median(event_length_seconds, na.rm = TRUE),
    mean_event_length_s     = mean(event_length_seconds,   na.rm = TRUE),
    max_event_length_s      = max(event_length_seconds,    na.rm = TRUE),
    min_median_running_z    = min(median_running_z,        na.rm = TRUE),
    median_median_running_z = median(median_running_z,     na.rm = TRUE),
    mean_median_running_z   = mean(median_running_z,       na.rm = TRUE),
    max_median_running_z    = max(median_running_z,        na.rm = TRUE)
  )]
}

desc_phase_mouse <- start_ev[, .summarize_lengths_and_z(.SD), by = .(phase_long, mouse)]

.round3 <- function(dt) {
  numcols <- names(dt)[vapply(dt, is.numeric, TRUE)]
  copy(dt)[, (numcols) := lapply(.SD, \(v) ifelse(is.finite(v), round(v, 3), v)),
           .SDcols = numcols][]
}
desc_phase_mouse_fmt <- .round3(desc_phase_mouse)
print(desc_phase_mouse_fmt)

#### Plot:           SF3_EventLengthMedianZRidgeplots####
p_lengths_by_mouse_phase <- ggplot(start_ev, aes(x = event_length_seconds, y = mouse, fill = mouse)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 0.5,
                      alpha = 0.85, scale = 1.10, rel_min_height = 0.01, linewidth = 0.3) +
  geom_point(aes(x = event_length_seconds, y = mouse),
             inherit.aes = FALSE, color = "grey50", size = 0.8, alpha = 0.85) +
  facet_wrap(~ phase_long, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = mouse_pal, guide = "none") +
  labs(title = "Wheel Running Event Lengths",
       x = "Wheel Running Event Length (s)", y = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

p_medianz_by_mouse_phase <- ggplot(start_ev, aes(x = median_running_z, y = mouse, fill = mouse)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 0.5,
                      alpha = 0.85, scale = 1.10, rel_min_height = 0.01, linewidth = 0.3) +
  geom_point(aes(x = median_running_z, y = mouse),
             inherit.aes = FALSE, color = "grey50", size = 0.8, alpha = 0.85) +
  facet_wrap(~ phase_long, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = mouse_pal, guide = "none") +
  labs(title = "Wheel Running Event Signal (Median z-score)",
       x = "Wheel Running Event Signal (Median z-score)", y = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p_lengths_by_mouse_phase)
print(p_medianz_by_mouse_phase)

#### Save:           SF3_EventLengthMedianZRidgeplots ####
fig_dir <- .fig_dir("SF3_EventLength_MedianZ_Ridgeplots")
fwrite(desc_phase_mouse_fmt,
       file.path(fig_dir, "WheelRunning_EventLength_and_MedianZ_by_Mouse_and_Phase.csv"))
ggsave(file.path(fig_dir, "WheelRunning_EventLengths_by_Mouse_and_Phase.tiff"),
       p_lengths_by_mouse_phase, width = 6.5, height = 6.0, units = "in",
       dpi = 600, bg = "white", device = "tiff", compression = "lzw")
ggsave(file.path(fig_dir, "WheelRunning_EventLengths_by_Mouse_and_Phase.svg"),
       p_lengths_by_mouse_phase, width = 6.5, height = 6.0, units = "in", bg = "white")
ggsave(file.path(fig_dir, "WheelRunning_MedianZ_by_Mouse_and_Phase.tiff"),
       p_medianz_by_mouse_phase, width = 6.5, height = 6.0, units = "in",
       dpi = 600, bg = "white", device = "tiff", compression = "lzw")
ggsave(file.path(fig_dir, "WheelRunning_MedianZ_by_Mouse_and_Phase.svg"),
       p_medianz_by_mouse_phase, width = 6.5, height = 6.0, units = "in", bg = "white")
##################################################################################################################
#### Figure 5a - Wheel running ACh amplitude z-score versus event length correlations ####
#### Transform:      Build OLS models ####

.per_mouse_corr <- function(ev_phase) {
  ev_phase[, {
    ok <- sum(is.finite(event_length_seconds) & is.finite(median_running_z)) >= 3 &&
      length(unique(event_length_seconds)) > 1 &&
      length(unique(median_running_z))     > 1
    pearson_r  <- if (ok) suppressWarnings(cor(event_length_seconds, median_running_z, method = "pearson"))  else NA_real_
    spearman_r <- if (ok) suppressWarnings(cor(event_length_seconds, median_running_z, method = "spearman")) else NA_real_
    .(n_events = .N, pearson_r = pearson_r, spearman_r = spearman_r)
  }, by = mouse][order(mouse)]
}

.fit_fe_ols <- function(ev_phase) {
  fit <- lm(median_running_z ~ event_length_seconds + mouse, data = ev_phase)
  cf  <- summary(fit)$coefficients
  est <- unname(cf["event_length_seconds","Estimate"])
  se  <- unname(cf["event_length_seconds","Std. Error"])
  df  <- max(1, df.residual(fit))
  t   <- est / se
  p   <- 2 * (1 - pt(abs(t), df))
  tcrit <- qt(0.975, df)
  data.table(model = "FE (mouse intercepts)",
             beta_length = est, se = se, df = df, t = t, p = p,
             ci95_lo = est - tcrit*se, ci95_hi = est + tcrit*se)
}

.fit_pooled_ols <- function(ev_phase) {
  fit <- lm(median_running_z ~ event_length_seconds, data = ev_phase)
  cf  <- summary(fit)$coefficients
  est <- unname(cf["event_length_seconds","Estimate"])
  se  <- unname(cf["event_length_seconds","Std. Error"])
  df  <- max(1, df.residual(fit))
  t   <- est / se
  p   <- 2 * (1 - pt(abs(t), df))
  tcrit <- qt(0.975, df)
  b0 <- unname(coef(fit)["(Intercept)"])
  data.table(model = "Pooled (black line)",
             beta_length = est, se = se, df = df, t = t, p = p,
             ci95_lo = est - tcrit*se, ci95_hi = est + tcrit*se,
             intercept = b0)
}

.per_mouse_ols <- function(ev_phase) {
  mice <- levels(droplevels(ev_phase$mouse))
  rbindlist(lapply(mice, function(m) {
    d <- ev_phase[mouse == m]
    if (nrow(d) < 2L || length(unique(d$event_length_seconds)) < 2L) {
      return(data.table(mouse = m, slope = NA_real_, intercept = NA_real_))
    }
    fm <- lm(median_running_z ~ event_length_seconds, data = d)
    data.table(mouse = m,
               slope = unname(coef(fm)["event_length_seconds"]),
               intercept = unname(coef(fm)["(Intercept)"]))
  }), use.names = TRUE, fill = TRUE)
}

.build_phase_results <- function(phase_label) {
  evp <- start_ev[phase_long == phase_label]
  list(
    per_mouse_corr  = .per_mouse_corr(evp),
    fe_results      = .fit_fe_ols(evp),
    pooled_results  = .fit_pooled_ols(evp),
    per_mouse_coefs = .per_mouse_ols(evp),
    data            = evp
  )
}

res_Aq <- .build_phase_results("Acquisition")
res_Ma <- .build_phase_results("Maintenance")

per_mouse_corr_table <- rbindlist(list(
  cbind(phase = "Acquisition",  res_Aq$per_mouse_corr),
  cbind(phase = "Maintenance", res_Ma$per_mouse_corr)
), use.names = TRUE)

model_results_table <- rbindlist(list(
  cbind(phase = "Acquisition",  rbind(res_Aq$fe_results, res_Aq$pooled_results, fill = TRUE)),
  cbind(phase = "Maintenance", rbind(res_Ma$fe_results, res_Ma$pooled_results, fill = TRUE))
), use.names = TRUE)

print(per_mouse_corr_table)
print(model_results_table)

#### Plot:           F5a_BoutLengthvsMedianZCorrelations####
.plot_permouse_with_pooled <- function(res_obj, title_text, ellipse_level = 0.90) {
  evp  <- res_obj$data
  pm   <- res_obj$per_mouse_coefs
  pool <- res_obj$pooled_results
  
  xr <- range(evp$event_length_seconds, na.rm = TRUE)
  
  pm_lines  <- pm[is.finite(slope) & is.finite(intercept),
                  .(x = xr, y = intercept + slope * xr), by = mouse]
  pooled_df <- data.table(x = xr, y = pool$intercept + pool$beta_length * xr)
  
  subtitle_txt <- sprintf("Pooled OLS (black): β = %.4f (95%% CI %.4f, %.4f), p = %.3g",
                          pool$beta_length, pool$ci95_lo, pool$ci95_hi, pool$p)
  
  ggplot(evp, aes(event_length_seconds, median_running_z,
                  color = mouse, shape = mouse)) +
    stat_ellipse(aes(group = mouse), type = "norm", level = ellipse_level,
                 linewidth = 0.6, alpha = 0.35, show.legend = FALSE) +
    geom_point(alpha = 0.9, size = 3) +
    geom_line(data = pm_lines, aes(x, y, color = mouse),
              linewidth = 0.9, inherit.aes = FALSE) +
    geom_line(data = pooled_df, aes(x, y),
              linewidth = 2, color = "black", inherit.aes = FALSE) +
    scale_color_manual(values = mouse_pal) +
    labs(title = title_text, subtitle = subtitle_txt,
         x = "Wheel Running Length (s)",
         y = "Median z-score (in-running)") +
    coord_cartesian(xlim = c(0, NA)) +
    theme_classic(base_size = 12) +
    theme(plot.title    = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
}

p_len_z_Aq <- .plot_permouse_with_pooled(res_Aq, "Acquisition")
p_len_z_Ma <- .plot_permouse_with_pooled(res_Ma, "Maintenance")
print(p_len_z_Aq); print(p_len_z_Ma)

#### Save:           F5a_BoutLengthvsMedianZCorrelations ####
fig_dir <- .fig_dir("F5a_BoutLengthvsMedianZCorrelations")
fwrite(per_mouse_corr_table,
       file.path(fig_dir, "BoutLength_vs_MedianZ_PerMouse_Pearson_Spearman.csv"))
fwrite(model_results_table,
       file.path(fig_dir, "BoutLength_vs_MedianZ_ModelResults_FE_and_Pooled.csv"))
ggsave(file.path(fig_dir, "Acquisition_BoutLength_vs_MedianZ_perMouseOLS_plusPooled.tiff"),
       p_len_z_Aq, width = 6, height = 4.5, units = "in",
       dpi = 600, bg = "white", device = "tiff", compression = "lzw")
ggsave(file.path(fig_dir, "Acquisition_BoutLength_vs_MedianZ_perMouseOLS_plusPooled.svg"),
       p_len_z_Aq, width = 6, height = 4.5, units = "in", bg = "white")
ggsave(file.path(fig_dir, "Maintenance_BoutLength_vs_MedianZ_perMouseOLS_plusPooled.tiff"),
       p_len_z_Ma, width = 6, height = 4.5, units = "in",
       dpi = 600, bg = "white", device = "tiff", compression = "lzw")
ggsave(file.path(fig_dir, "Maintenance_BoutLength_vs_MedianZ_perMouseOLS_plusPooled.svg"),
       p_len_z_Ma, width = 6, height = 4.5, units = "in", bg = "white")


##################################################################################################################
#### Figure 5bc - Wheel running event length retention curves ####
#### Transform:      Build event-length retention table ####

L_max  <- floor(max(start_ev$event_length_seconds, na.rm = TRUE))
L_grid <- 0:L_max

decay_curve <- start_ev[
  is.finite(event_length_seconds),
  {
    pct <- 100 * vapply(L_grid, function(L) mean(event_length_seconds >= L), numeric(1))
    data.table(L = L_grid, retained_pct = pct)
  },
  by = .(mouse, phase_long)
]
decay_curve[, `:=`(
  L            = as.numeric(L),
  retained_pct = as.numeric(retained_pct),
  mouse        = factor(as.character(mouse), levels = mouse_lookup$label)
)]
setorder(decay_curve, mouse, L)

#### Plot:           F5bc_BoutLengthRetentionCurves####
p_decay <- ggplot(decay_curve, aes(L, retained_pct, color = mouse)) +
  geom_step(linewidth = 1) +
  facet_wrap(~ phase_long, ncol = 1, scales = "free_x") +
  scale_color_manual(values = mouse_pal, name = "Mouse") +
  scale_x_continuous(breaks = seq(0, L_max, 60),
                     minor_breaks = seq(0, L_max, 15),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100), minor_breaks = NULL) +
  labs(title = "Event-length decay (retention) curves",
       subtitle = "% of running events with length >= L",
       x = "Event length L (s)", y = "Events retained (%)") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
        panel.grid.minor.x = element_line(color = "grey93", linewidth = 0.3),
        panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
        strip.background = element_rect(fill = "white", color = NA))
print(p_decay)

p_decay_zoom <- ggplot(decay_curve, aes(L, retained_pct, color = mouse)) +
  geom_step(linewidth = 1) +
  facet_wrap(~ phase_long, ncol = 1, scales = "free_x") +
  scale_color_manual(values = mouse_pal, name = "Mouse") +
  scale_x_continuous(breaks = seq(0, 60, 15), limits = c(0, 60), minor_breaks = NULL) +
  scale_y_continuous(breaks = seq(50, 100, 25), limits = c(50, 100), minor_breaks = NULL) +
  labs(title = "Event-length decay (retention) curves",
       subtitle = "% of running events with length >= L",
       x = "Event length L (s)", y = "Events retained (%)") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
        panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
        strip.background = element_rect(fill = "white", color = NA))
print(p_decay_zoom)

#### Save:           F5bc_BoutLengthRetentionCurves ####
fig_dir <- .fig_dir("F5bc_BoutLengthRetentionCurves")
ggsave(file.path(fig_dir, "EventLength_RetentionCurves.tiff"),
       p_decay, width = 6.5, height = 5.5, units = "in",
       dpi = 600, compression = "lzw", bg = "white")
ggsave(file.path(fig_dir, "EventLength_RetentionCurves.svg"),
       p_decay, width = 6.5, height = 5.5, units = "in",
       device = svglite::svglite, bg = "white")
ggsave(file.path(fig_dir, "EventLength_RetentionCurves_Zoom0_60s.tiff"),
       p_decay_zoom, width = 6.5, height = 5.5, units = "in",
       dpi = 600, compression = "lzw", bg = "white")
ggsave(file.path(fig_dir, "EventLength_RetentionCurves_Zoom0_60s.svg"),
       p_decay_zoom, width = 6.5, height = 5.5, units = "in",
       device = svglite::svglite, bg = "white")
##################################################################################################################
#### Supplemental Figure 4 - Wheel running Mouse x Phase Trends ####
#### Helpers:        Grid enforcement, short-run drop, baseline-z ####

enforce_grid_25hz <- function(dt, win, hz = 25, value_col = "Ratio") {
  step <- 1 / hz
  grid <- round(seq(win[1], win[2], by = step), 2)
  setDT(dt)
  dt <- dt[is.finite(Event_time) & is.finite(get(value_col))]
  dt[, Event_time := round(as.numeric(Event_time), 2)]
  dt <- dt[Event_time >= win[1] & Event_time <= win[2]]
  setorder(dt, Event, Event_time)
  dt <- unique(dt, by = c("Event", "Event_time"))
  need_n <- length(grid)
  keep_events <- dt[, {
    times <- sort(unique(Event_time))
    .(ok = length(times) == need_n && identical(times, grid))
  }, by = Event][ok == TRUE, Event]
  dt[Event %in% keep_events][]
}

drop_short_running <- function(dt, min_run = 25) {
  setDT(dt)
  runlen <- dt[Interval_phase == "Running",
               .(run_len = max(Event_time) - min(Event_time)), by = Event]
  keep <- runlen[is.finite(run_len) & run_len >= min_run, Event]
  dt[Event %in% keep][]
}

z_event_then_center <- function(dt, baseline_win, value_col = "Ratio",
                                out_col = "Ratio_z_event") {
  setDT(dt)
  dt[, (out_col) := {
    v      <- get(value_col)
    mu_all <- mean(v, na.rm = TRUE)
    sd_all <- sd(v,  na.rm = TRUE)
    z_all  <- if (is.finite(sd_all) && sd_all > 0) (v - mu_all) / sd_all else 0 * v
    base_m <- mean(z_all[Event_time >= baseline_win[1] & Event_time <= baseline_win[2]],
                   na.rm = TRUE)
    z_all - base_m
  }, by = Event]
  dt[]
}

process_series <- function(nm) {
  dt <- get(nm, inherits = TRUE); setDT(dt)
  is_start <- grepl("^Start_", nm)
  if (is_start) { win <- c(-15, 30); base <- c(-15, -10) }
  else          { win <- c(-30, 15); base <- c(10,   15) }
  dt <- enforce_grid_25hz(dt, win = win, hz = 25, value_col = "Ratio")
  dt <- drop_short_running(dt, min_run = 25)
  dt <- z_event_then_center(dt, baseline_win = base, value_col = "Ratio",
                            out_col = "Ratio_z_event")
  assign(nm, dt, envir = .GlobalEnv)
  invisible(dt[])
}

for (nm in ls(pattern = "^(Start_|Stop_)")) process_series(nm)

#### Plot function ####
plot_series_trend <- function(dt, subtitle = "", lw_scale = 1, base_size = 16) {
  setDT(dt)
  dt  <- dt[is.finite(Event_time) & is.finite(Ratio_z_event)]
  xlo <- floor(min(dt$Event_time,   na.rm = TRUE))
  xhi <- ceiling(max(dt$Event_time, na.rm = TRUE))
  
  ggplot() +
    geom_line(data = dt,
              aes(Event_time, Ratio_z_event, group = Event, color = factor(Event)),
              linewidth = 0.6 * lw_scale, alpha = 0.85, show.legend = FALSE) +
    stat_summary(data = dt, aes(Event_time, Ratio_z_event),
                 fun.data = mean_cl_normal, geom = "ribbon", inherit.aes = FALSE,
                 fill = "grey70", alpha = 0.8, na.rm = TRUE) +
    stat_summary(data = dt, aes(Event_time, Ratio_z_event),
                 fun = mean, geom = "line", inherit.aes = FALSE,
                 color = "white", linewidth = 2.2 * lw_scale,
                 lineend = "round", na.rm = TRUE) +
    stat_summary(data = dt, aes(Event_time, Ratio_z_event),
                 fun = mean, geom = "line", inherit.aes = FALSE,
                 color = "black", linewidth = 1.25 * lw_scale,
                 lineend = "round", na.rm = TRUE) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed",
               linewidth = 1.2 * lw_scale) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8 * lw_scale) +
    coord_cartesian(xlim = c(xlo, xhi), ylim = c(-5, 5), expand = FALSE) +
    scale_x_continuous(breaks = seq(xlo, xhi, 5)) +
    scale_y_continuous(breaks = seq(-5, 5, 1)) +
    scale_color_scico_d(palette = "bamako", begin = 0.1, end = 0.8,
                        direction = 1, guide = "none") +
    labs(title    = "Acetylcholine - Ventral Dentate Gyrus",
         subtitle = subtitle,
         x        = "Event Time (s)",
         y        = "z-score") +
    theme_classic(base_size = base_size) +
    theme(plot.title    = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
}

#### SF4a: Mouse x Session (diagnostic only) ####
series_names <- ls(pattern = "^(Start_|Stop_)")
fig_dir      <- .fig_dir("SF4a_TrendsMouseSession")

for (nm in series_names) {
  dt <- get(nm, inherits = TRUE)
  if (!("Ratio_z_event" %in% names(dt)) || !nrow(dt)) next
  p <- plot_series_trend(dt, subtitle = nm)
  ggsave(file.path(fig_dir, paste0(nm, ".tiff")), p,
         device = "tiff", width = 6.5, height = 5.0, units = "in",
         dpi = 600, compression = "lzw", bg = "white")
}

#### Table:          Event counts per Mouse x Phase ####
event_counts <- rbindlist(lapply(ls(pattern = "^(Start_|Stop_)\\d{4}_"), function(nm) {
  dt <- get(nm, inherits = TRUE)
  if (!"Ratio_z_event" %in% names(dt) || !nrow(dt)) return(NULL)
  data.table(kind     = sub("^(Start|Stop)_.*$", "\\1", nm),
             series   = nm,
             mouse    = .mouse_from(nm),
             phase    = .phase_code(nm),
             n_events = uniqueN(dt$Event))
}))
event_counts[, phase_long := .phase_long(phase)]
setorder(event_counts, kind, mouse, phase, series)

event_counts_summary <- event_counts[, .(
  n_series     = .N,
  total_events = sum(n_events)
), by = .(kind, mouse, phase_long)]
print(event_counts_summary)

#### Transform:      Mouse x Phase groups ####
mouse_phase_data <- list()
for (kind in c("Start", "Stop")) {
  nms <- ls(pattern = paste0("^", kind, "_\\d{4}_"), envir = .GlobalEnv)
  nms <- nms[sapply(nms, function(nm) {
    d <- get(nm, inherits = TRUE)
    "Ratio_z_event" %in% names(d) && nrow(d) > 0
  })]
  if (!length(nms)) next
  
  combos <- unique(data.table(
    nm    = nms,
    mouse = .mouse_from(nms),
    phase = .phase_code(nms)
  ))[, key := paste(kind, mouse, phase, sep = "_")]
  
  for (k in unique(combos$key)) {
    grp <- combos[key == k]
    dt_combined <- rbindlist(lapply(grp$nm, function(nm) {
      d <- get(nm, inherits = TRUE)
      d[is.finite(Event_time) & is.finite(Ratio_z_event),
        .(Event         = paste0(nm, "_", Event),
          Event_time    = Event_time,
          Ratio_z_event = Ratio_z_event)]
    }), use.names = TRUE, fill = TRUE)
    
    n_ev       <- uniqueN(dt_combined$Event)
    phase_long <- .phase_long(grp$phase[1])
    mouse_phase_data[[k]] <- list(
      dt         = dt_combined,
      mouse      = grp$mouse[1],
      phase      = grp$phase[1],
      phase_long = phase_long,
      kind       = kind,
      n_events   = n_ev,
      fname      = paste0(kind, "_", grp$mouse[1], "_", grp$phase[1]),
      subtitle   = paste0(phase_long, " - ", grp$mouse[1],
                          " (n = ", n_ev, " events)")
    )
  }
}

#### Save:           SF4_TrendsMousexPhase  ####
fig_dir <- .fig_dir("SF4b_TrendsMousePhase")
fwrite(event_counts,         file.path(fig_dir, "EventCounts_per_Series.csv"))
fwrite(event_counts_summary, file.path(fig_dir, "EventCounts_per_MousePhase.csv"))

for (k in names(mouse_phase_data)) {
  fname <- mouse_phase_data[[k]]$fname
  p_full <- plot_series_trend(mouse_phase_data[[k]]$dt,
                              subtitle = mouse_phase_data[[k]]$subtitle)
  p_bare <- plot_series_trend(mouse_phase_data[[k]]$dt,
                              subtitle  = "",
                              lw_scale  = 0.15,
                              base_size = 5) +
    labs(title = NULL, subtitle = NULL) +
    theme(plot.title    = element_blank(),
          plot.subtitle = element_blank(),
          axis.title    = element_text(size = 6, face = "bold"))
  
  # Full size
  ggsave(file.path(fig_dir, paste0(fname, ".tiff")), p_full,
         device = "tiff", width = 6.5, height = 5.0, units = "in",
         dpi = 600, compression = "lzw", bg = "white")
  ggsave(file.path(fig_dir, paste0(fname, ".svg")), p_full,
         device = svglite::svglite, width = 6.5, height = 5.0, units = "in",
         bg = "white")
  
  # Bare (4 x 4 cm)
  ggsave(file.path(fig_dir, paste0(fname, "_bare.tiff")), p_bare,
         device = "tiff", width = 4, height = 4, units = "cm",
         dpi = 1200, compression = "lzw", bg = "white")
  ggsave(file.path(fig_dir, paste0(fname, "_bare.svg")), p_bare,
         device = svglite::svglite, width = 4, height = 4, units = "cm",
         bg = "white")
}



##################################################################################################################
#### Figure 6 - Wheel running nested phase event trends (nested avg mouse x phase trends) ####
#### Transform:      Group Start_/Stop_ objects by Mouse × Phase ####

start_objs <- ls(pattern = "^Start_\\d{4}_", envir = .GlobalEnv)
stop_objs  <- ls(pattern = "^Stop_\\d{4}_",  envir = .GlobalEnv)

groups_start <- if (length(start_objs))
  split(start_objs, paste("Start", .mouse_from(start_objs), .phase_code(start_objs), sep = "_"))

groups_stop  <- if (length(stop_objs))
  split(stop_objs,  paste("Stop",  .mouse_from(stop_objs),  .phase_code(stop_objs),  sep = "_"))


#### Transform:      Collapse each Mouse × Phase group to a per-time mouse trace ####
.combine_group_to_mouse_phase <- function(key, members, kind = c("Start","Stop")) {
  kind <- match.arg(kind)
  pieces <- lapply(members, function(nm) {
    dt <- as.data.table(get(nm, inherits = TRUE))
    if (!all(c("Event","Event_time","Ratio_z_event") %in% names(dt))) return(NULL)
    dt[is.finite(Event_time) & is.finite(Ratio_z_event),
       .(Event_time    = as.numeric(Event_time),
         Event         = as.character(Event),
         Ratio_z_event = as.numeric(Ratio_z_event))]
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(NULL)
  
  combo <- rbindlist(pieces, use.names = TRUE, fill = TRUE)
  out   <- combo[, .(z_mean = mean(Ratio_z_event, na.rm = TRUE)),
                 by = Event_time][order(Event_time)]
  
  mouse <- sub(paste0("^", kind, "_(Mouse\\d+)_.*$"), "\\1", key)
  phase <- sub(paste0("^", kind, "_Mouse\\d+_(Aq|Ma)$"), "\\1", key)
  out[, `:=`(mouse = mouse, phase = phase,
             phase_long = .phase_long(phase), kind = kind)]
  out[]
}

start_lines <- rbindlist(lapply(names(groups_start),
                                function(k) .combine_group_to_mouse_phase(k, groups_start[[k]], "Start")),
                         use.names = TRUE, fill = TRUE)
stop_lines  <- rbindlist(lapply(names(groups_stop),
                                function(k) .combine_group_to_mouse_phase(k, groups_stop[[k]],  "Stop")),
                         use.names = TRUE, fill = TRUE)

#### Plot:           F6_TrendsPhase ####
.plot_phase_trends <- function(lines_dt, phase_code = c("Aq","Ma"),
                               kind = c("Start","Stop")) {
  phase_code <- match.arg(phase_code)
  kind       <- match.arg(kind)
  
  df_m   <- lines_dt[phase == phase_code]
  df_avg <- df_m[, .(z_mean = mean(z_mean, na.rm = TRUE)), by = Event_time][order(Event_time)]
  
  xlim         <- if (kind == "Start") c(-15, 30) else c(-30, 15)
  title_txt    <- if (phase_code == "Aq") "Acquisition Phase" else "Maintenance Phase"
  subtitle_txt <- if (kind == "Start")  "Running Initiation" else "Running Termination"
  
  ggplot() +
    geom_line(data = df_m,
              aes(Event_time, z_mean, color = mouse), linewidth = 1) +
    geom_line(data = df_avg, aes(Event_time, z_mean),
              color = "white", linewidth = 2.4, lineend = "round") +
    geom_line(data = df_avg, aes(Event_time, z_mean),
              color = "black", linewidth = 1.35, lineend = "round") +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
    coord_cartesian(xlim = xlim, ylim = c(-2, 5), expand = FALSE) +
    scale_x_continuous(breaks = seq(xlim[1], xlim[2], 5)) +
    scale_y_continuous(breaks = seq(-2, 5, 1)) +
    scale_color_manual(values = mouse_pal, name = "Mouse Average") +
    labs(title = title_txt, subtitle = subtitle_txt,
         x = "Event Time (s)", y = "z-score") +
    theme_classic(base_size = 12) +
    theme(plot.title    = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
}

p_start_aq <- if (nrow(start_lines)) .plot_phase_trends(start_lines, "Aq", "Start")
p_start_ma <- if (nrow(start_lines)) .plot_phase_trends(start_lines, "Ma", "Start")
p_stop_aq  <- if (nrow(stop_lines))  .plot_phase_trends(stop_lines,  "Aq", "Stop")
p_stop_ma  <- if (nrow(stop_lines))  .plot_phase_trends(stop_lines,  "Ma", "Stop")

if (!is.null(p_start_aq)) print(p_start_aq)
if (!is.null(p_start_ma)) print(p_start_ma)
if (!is.null(p_stop_aq))  print(p_stop_aq)
if (!is.null(p_stop_ma))  print(p_stop_ma)

#### Save:           F6_TrendsPhase ####
fig_dir <- .fig_dir("F6_TrendsPhase")
.save_tiff_svg <- function(p, fname, w = 6.5, h = 4.8) {
  ggsave(file.path(fig_dir, paste0(fname, ".tiff")),
         p, width = w, height = h, units = "in",
         dpi = 600, bg = "white", device = "tiff", compression = "lzw")
  ggsave(file.path(fig_dir, paste0(fname, ".svg")),
         p, width = w, height = h, units = "in",
         device = svglite::svglite, bg = "white")
}
if (!is.null(p_start_aq)) .save_tiff_svg(p_start_aq, "Start_Trends_Phase_Acquisition")
if (!is.null(p_start_ma)) .save_tiff_svg(p_start_ma, "Start_Trends_Phase_Maintenance")
if (!is.null(p_stop_aq))  .save_tiff_svg(p_stop_aq,  "Stop_Trends_Phase_Acquisition")
if (!is.null(p_stop_ma))  .save_tiff_svg(p_stop_ma,  "Stop_Trends_Phase_Maintenance")


##################################################################################################################
#### Supplemental Figure 5 - Per-event slope across peri-event timeline ####
#### Transform:      Define 5-s windows spanning -15 through +15 around event ####

start_window_bounds <- seq(-15, 30, by = 5)
stop_window_bounds  <- seq(-30, 15, by = 5)

start_windows <- data.table(
  win_start = head(start_window_bounds, -1),
  win_end   = tail(start_window_bounds, -1)
)[, win_mid := (win_start + win_end) / 2]

stop_windows <- data.table(
  win_start = head(stop_window_bounds, -1),
  win_end   = tail(stop_window_bounds, -1)
)[, win_mid := (win_start + win_end) / 2]

#### Transform:      Compute per-event slope in each window ####
.slope_per_event <- function(series_names, windows, kind_label) {
  rbindlist(lapply(series_names, function(nm) {
    dt <- as.data.table(get(nm, inherits = TRUE))
    if (!all(c("Event", "Event_time", "Ratio_z_event") %in% names(dt))) return(NULL)
    if (!nrow(dt)) return(NULL)
    
    mouse_id <- .mouse_from(nm)
    phase_id <- .phase_code(nm)
    
    rbindlist(lapply(seq_len(nrow(windows)), function(i) {
      w <- windows[i]
      sub <- dt[is.finite(Event_time) & is.finite(Ratio_z_event) &
                  Event_time >= w$win_start & Event_time < w$win_end]
      if (!nrow(sub)) return(NULL)
      sub[, {
        if (.N < 3 || length(unique(Event_time)) < 2) {
          .(slope = NA_real_, n_tp = .N)
        } else {
          fm <- lm(Ratio_z_event ~ Event_time)
          .(slope = unname(coef(fm)["Event_time"]), n_tp = .N)
        }
      }, by = Event][, `:=`(
        win_start  = w$win_start,
        win_end    = w$win_end,
        win_mid    = w$win_mid,
        series     = nm,
        mouse      = mouse_id,
        phase      = phase_id,
        phase_long = .phase_long(phase_id),
        kind       = kind_label
      )][]
    }), fill = TRUE)
  }), fill = TRUE)
}

start_names <- ls(pattern = "^Start_\\d{4}_", envir = .GlobalEnv)
stop_names  <- ls(pattern = "^Stop_\\d{4}_",  envir = .GlobalEnv)

slope_start <- .slope_per_event(start_names, start_windows, "Start")
slope_stop  <- .slope_per_event(stop_names,  stop_windows,  "Stop")
slope_all   <- rbindlist(list(slope_start, slope_stop), fill = TRUE)
slope_all   <- slope_all[!is.na(slope)]

#### Transform:      Per-mouse and across-mouse summaries ####
slope_mouse <- slope_all[, .(mean_slope = mean(slope, na.rm = TRUE),
                             n_events   = .N),
                         by = .(mouse, win_mid, kind, phase_long)]

slope_pool <- slope_mouse[, .(pooled_slope = mean(mean_slope, na.rm = TRUE),
                              se_slope     = sd(mean_slope, na.rm = TRUE) / sqrt(.N)),
                          by = .(win_mid, kind, phase_long)]

for (d in list(slope_all, slope_mouse, slope_pool)) {
  d[, `:=`(phase_long = factor(phase_long, levels = c("Acquisition","Maintenance")),
           kind       = factor(kind,       levels = c("Start","Stop")))]
}

rect_baseline <- data.table(
  kind       = factor(c("Start","Start","Stop","Stop"), levels = c("Start","Stop")),
  phase_long = factor(c("Acquisition","Maintenance","Acquisition","Maintenance"),
                      levels = c("Acquisition","Maintenance")),
  xmin       = c(-15, -15, 10, 10),
  xmax       = c(-10, -10, 15, 15)
)

#### Plot:           SF5_BaselineSelection ####
p_slope_timeline <- ggplot() +
  geom_rect(data = rect_baseline,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey60", alpha = 0.2, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8, color = "black") +
  geom_point(data = slope_all,
             aes(x = win_mid, y = slope, color = mouse),
             size = 0.7, alpha = 0.25,
             position = position_jitter(width = 0.4, height = 0)) +
  geom_line(data = slope_mouse,
            aes(x = win_mid, y = mean_slope, color = mouse, group = mouse),
            linewidth = 0.9, alpha = 0.85) +
  geom_line(data = slope_pool,
            aes(x = win_mid, y = pooled_slope, group = 1),
            color = "black", linewidth = 1.4) +
  geom_point(data = slope_pool,
             aes(x = win_mid, y = pooled_slope),
             color = "black", size = 2.5) +
  facet_grid(phase_long ~ kind, scales = "free_x") +
  scale_x_continuous(breaks = function(x) seq(floor(min(x)/5)*5,
                                              ceiling(max(x)/5)*5,
                                              by = 5)) +
  scale_color_manual(values = mouse_pal, name = "Mouse") +
  labs(title    = "Per-event slope across peri-event timeline (5-s windows)",
       subtitle = "Grey rectangle = baseline window used in analysis",
       x        = "Window midpoint (s)",
       y        = "Slope of z-score ~ time (z/s)") +
  theme_classic(base_size = 12) +
  theme(plot.title        = element_text(hjust = 0.5),
        plot.subtitle     = element_text(hjust = 0.5),
        strip.background  = element_rect(fill = "white", color = NA),
        panel.spacing.x   = unit(0.75, "cm"))

print(p_slope_timeline)

#### Save:           SF5_BaselineSelection ####
fig_dir <- .fig_dir("SF5_BaselineSelection")
fwrite(slope_all,   file.path(fig_dir, "SlopeTimeline_PerEvent.csv"))
fwrite(slope_mouse, file.path(fig_dir, "SlopeTimeline_PerMouse.csv"))
fwrite(slope_pool,  file.path(fig_dir, "SlopeTimeline_Pooled.csv"))
ggsave(file.path(fig_dir, "SlopeTimeline.tiff"),
       p_slope_timeline, width = 9.0, height = 6.0, units = "in",
       dpi = 600, device = "tiff", compression = "lzw", bg = "white")
ggsave(file.path(fig_dir, "SlopeTimeline.svg"),
       p_slope_timeline, width = 9.0, height = 6.0, units = "in",
       device = svglite::svglite, bg = "white")
##################################################################################################################
#### Figure 6 Stats - Paired t-tests for Wheel running phase event trends (baseline window vs running window) ####
#### Transform:      Compute per-mouse window means (average across events) ####

.mouse_window_means <- function(dt, base_win, run_win) {
  setDT(dt)
  dt <- dt[is.finite(Event_time)]
  
  by_ev <- dt[, .(
    base_z   = mean(Ratio_z_event[Event_time >= base_win[1] & Event_time <= base_win[2]],
                    na.rm = TRUE),
    run_z    = mean(Ratio_z_event[Event_time >= run_win[1]  & Event_time <= run_win[2]],
                    na.rm = TRUE),
    base_raw = mean(Ratio[Event_time >= base_win[1] & Event_time <= base_win[2]],
                    na.rm = TRUE),
    run_raw  = mean(Ratio[Event_time >= run_win[1]  & Event_time <= run_win[2]],
                    na.rm = TRUE)
  ), by = Event][is.finite(base_z) & is.finite(run_z) &
                   is.finite(base_raw) & is.finite(run_raw)]
  if (!nrow(by_ev)) return(NULL)
  
  data.table(
    base_z   = mean(by_ev$base_z),
    run_z    = mean(by_ev$run_z),
    base_raw = mean(by_ev$base_raw),
    run_raw  = mean(by_ev$run_raw),
    n_events = nrow(by_ev)
  )
}

.collect_rows <- function(kind = c("Start","Stop"), phase = c("Aq","Ma")) {
  kind  <- match.arg(kind); phase <- match.arg(phase)
  glist <- if (kind == "Start") groups_start else groups_stop
  keys  <- names(glist)
  keys  <- keys[grepl(paste0("^", kind, "_Mouse\\d+_", phase, "$"), keys)]
  if (!length(keys)) return(NULL)
  
  if (kind == "Start") { base_win <- c(-15, -10); run_win <- c(0, 30)  }
  else                 { base_win <- c(10,  15);  run_win <- c(-30, 0) }
  
  rows <- rbindlist(lapply(keys, function(key) {
    members <- glist[[key]]
    dt <- rbindlist(lapply(members, function(nm) as.data.table(get(nm, inherits = TRUE))),
                    use.names = TRUE, fill = TRUE)
    mm <- .mouse_window_means(dt, base_win, run_win)
    if (is.null(mm)) return(NULL)
    data.table(mouse = sub(paste0("^", kind, "_(Mouse\\d+)_.*$"), "\\1", key), mm)
  }), use.names = TRUE, fill = TRUE)
  
  if (!nrow(rows)) return(NULL)
  rows[, `:=`(
    diff_z  = run_z - base_z,
    pct_raw = ifelse(abs(base_raw) > .Machine$double.eps,
                     100 * (run_raw - base_raw) / abs(base_raw), NA_real_)
  )][]
}

.run_test <- function(kind, phase, rows) {
  if (is.null(rows) || nrow(rows) < 2L) return(NULL)
  d  <- rows$diff_z
  n  <- length(d)
  tt <- t.test(rows$run_z, rows$base_z, paired = TRUE)
  dz <- mean(d) / sd(d)
  g  <- if (n > 2) dz * (1 - 3/(4*n - 9)) else NA_real_
  
  data.table(
    kind                = kind,
    phase_code          = phase,
    phase_long          = .phase_long(phase),
    n_mice              = n,
    mean_base_z         = mean(rows$base_z),
    mean_run_z          = mean(rows$run_z),
    mean_diff_z         = mean(d),
    t                   = unname(tt$statistic),
    df                  = unname(tt$parameter),
    p_value             = unname(tt$p.value),
    ci_lo               = unname(tt$conf.int[1]),
    ci_hi               = unname(tt$conf.int[2]),
    cohen_dz            = dz,
    hedges_g            = g,
    mean_pct_change_raw = mean(rows$pct_raw, na.rm = TRUE),
    min_pct_change_raw  = min(rows$pct_raw,  na.rm = TRUE),
    max_pct_change_raw  = max(rows$pct_raw,  na.rm = TRUE)
  )
}

rows_Start_Aq <- .collect_rows("Start","Aq")
rows_Start_Ma <- .collect_rows("Start","Ma")
rows_Stop_Aq  <- .collect_rows("Stop","Aq")
rows_Stop_Ma  <- .collect_rows("Stop","Ma")

sum_Start_Aq <- .run_test("Start","Aq", rows_Start_Aq)
sum_Start_Ma <- .run_test("Start","Ma", rows_Start_Ma)
sum_Stop_Aq  <- .run_test("Stop","Aq",  rows_Stop_Aq)
sum_Stop_Ma  <- .run_test("Stop","Ma",  rows_Stop_Ma)

#### Print:          F6_TrendsPhaseStats####
.print_block <- function(tag, rows, summary) {
  cat("\n", tag, " — per-mouse means (event-averaged)\n", sep = "")
  if (is.null(rows)) { cat("  (insufficient mice)\n"); return() }
  print(rows[, .(mouse, n_events, base_z, run_z, diff_z, base_raw, run_raw, pct_raw)])
  cat("\n", tag, " — paired t-test across mice (z; in-running vs baseline)\n", sep = "")
  print(summary)
}
.print_block("START / Acquisition", rows_Start_Aq, sum_Start_Aq)
.print_block("START / Maintenance", rows_Start_Ma, sum_Start_Ma)
.print_block("STOP  / Acquisition", rows_Stop_Aq,  sum_Stop_Aq)
.print_block("STOP  / Maintenance", rows_Stop_Ma,  sum_Stop_Ma)

#### Save:           F6_TrendsPhase ####
fig_dir <- .fig_dir("F6_TrendsPhase")
.save_pair <- function(rows, summary, stem) {
  if (is.null(rows) || is.null(summary)) return(invisible())
  fwrite(rows,    file.path(fig_dir, paste0(stem, "_PerMouseMeans.csv")))
  fwrite(summary, file.path(fig_dir, paste0(stem, "_PairedT_Summary.csv")))
}
.save_pair(rows_Start_Aq, sum_Start_Aq, "PhaseStats_START_Acquisition")
.save_pair(rows_Start_Ma, sum_Start_Ma, "PhaseStats_START_Maintenance")
.save_pair(rows_Stop_Aq,  sum_Stop_Aq,  "PhaseStats_STOP_Acquisition")
.save_pair(rows_Stop_Ma,  sum_Stop_Ma,  "PhaseStats_STOP_Maintenance")


##################################################################################################################
#### Supplemental Figure 7 and Stats - RMANOVA of phase trends with varying bin size  ####
#### Transform:      Harmonize inputs for binning ####

.norm_cols <- function(dt) {
  setDT(dt)
  nm <- names(dt)
  if ("Event_time" %in% nm) setnames(dt, "Event_time", "t")
  if (!"z_mean" %in% nm && "Ratio_z_event" %in% nm) setnames(dt, "Ratio_z_event", "z_mean")
  dt[, .(t = as.numeric(t), z_mean = as.numeric(z_mean), mouse = as.factor(mouse))]
}

Start_Aq_adj <- .norm_cols(start_lines[phase == "Aq",
                                       .(t = Event_time, z_mean = z_mean, mouse = factor(mouse))])
Start_Ma_adj <- .norm_cols(start_lines[phase == "Ma",
                                       .(t = Event_time, z_mean = z_mean, mouse = factor(mouse))])
Stop_Aq_adj  <- .norm_cols(stop_lines[phase == "Aq",
                                      .(t = Event_time, z_mean = z_mean, mouse = factor(mouse))])
Stop_Ma_adj  <- .norm_cols(stop_lines[phase == "Ma",
                                      .(t = Event_time, z_mean = z_mean, mouse = factor(mouse))])

#### Transform:      Bin summary per mouse ####
bin_summarize <- function(dt, bin_width = 5, range = c(-15, 30),
                          stat = c("mean", "auc")) {
  stat <- match.arg(stat)
  dt <- .norm_cols(dt)[is.finite(t) & is.finite(z_mean)]
  dt <- dt[t >= range[1] & t <= range[2]]
  setorder(dt, mouse, t)
  
  step <- median(diff(sort(unique(dt$t))))
  brks <- seq(range[1], range[2], by = bin_width)
  if (tail(brks, 1) < range[2]) brks <- c(brks, range[2])
  
  dt[, bin := cut(t, breaks = brks, right = FALSE, include.lowest = TRUE)]
  dt <- dt[!is.na(bin)]
  
  dt[, .(value = if (stat == "mean") mean(z_mean, na.rm = TRUE)
         else sum(z_mean, na.rm = TRUE) * step),
     by = .(mouse, bin)][, bin := droplevels(bin)][]
}

#### Stats:          RMANOVA with GG correction ####
rm_anova <- function(bins) {
  bins <- as.data.frame(bins)
  bins$mouse <- factor(bins$mouse)
  bins$bin   <- factor(bins$bin)
  fit <- afex::aov_ez(id = "mouse", dv = "value", data = bins,
                      within = "bin", type = 3,
                      anova_table = list(correction = "GG", es = "ges"))
  list(kind = "afex", nice = afex::nice(fit, es = "ges", correction = "GG"))
}

#### Stats:          Unadjusted paired t-tests for all bin pairs ####
param_posthoc_allpairs <- function(bins_dt) {
  bins <- as.data.table(bins_dt)
  bins[, bin := factor(bin, levels = levels(bins$bin))]
  wide <- dcast(bins, mouse ~ bin, value.var = "value")
  lev  <- levels(bins$bin)
  cmb  <- t(combn(lev, 2))
  
  rbindlist(lapply(seq_len(nrow(cmb)), function(i) {
    b1 <- cmb[i, 1]; b2 <- cmb[i, 2]
    x <- wide[[b1]]; y <- wide[[b2]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 2) {
      return(data.table(bin1 = b1, bin2 = b2, n = sum(ok),
                        mean_diff = NA_real_, sd_diff = NA_real_,
                        cohen_dz = NA_real_, t = NA_real_,
                        df = NA_real_, p_raw = NA_real_))
    }
    d  <- y[ok] - x[ok]
    tt <- t.test(y[ok], x[ok], paired = TRUE)
    data.table(bin1 = b1, bin2 = b2, n = sum(ok),
               mean_diff = mean(d), sd_diff = sd(d),
               cohen_dz  = mean(d) / sd(d),
               t  = unname(tt$statistic),
               df = unname(tt$parameter),
               p_raw = unname(tt$p.value))
  }), use.names = TRUE, fill = TRUE)[order(bin1, bin2)]
}

cld_from_posthoc <- function(posthoc) {
  if (is.null(posthoc) || !nrow(posthoc))
    return(data.table(period = character(), letters = character()))
  bins_all <- sort(unique(c(as.character(posthoc$bin1), as.character(posthoc$bin2))))
  id_map   <- setNames(sprintf("b%02d", seq_along(bins_all)), bins_all)
  nm <- paste(id_map[as.character(posthoc$bin1)],
              id_map[as.character(posthoc$bin2)], sep = "-")
  pv <- setNames(posthoc$p_raw, nm)
  pv <- pv[is.finite(pv)]
  if (!length(pv)) return(data.table(period = bins_all, letters = "A"))
  letters_id  <- multcompView::multcompLetters(pv)$Letters
  letters_vec <- letters_id[id_map[bins_all]]
  data.table(period = bins_all, letters = unname(letters_vec))
}

#### Plot:           Function for creating bar plots and CLD annotations####
plot_bin_bars <- function(bins, label, cld = NULL) {
  sm <- bins[, .(m = mean(value), sd = sd(value), n = .N), by = bin][order(bin)]
  sm[, se    := sd / sqrt(pmax(n, 1L))]
  sm[, tcrit := qt(0.975, df = pmax(n - 1L, 1L))]
  sm[, `:=`(lo = m - tcrit * se, hi = m + tcrit * se)]
  bins[, bin := factor(bin, levels = levels(sm$bin))]
  
  ann <- NULL
  if (!is.null(cld) && nrow(cld)) {
    ann <- merge(data.table(bin = levels(sm$bin)), cld,
                 by.x = "bin", by.y = "period", all.x = TRUE)
    ann <- merge(ann, sm[, .(bin, hi)], by = "bin", all.x = TRUE)
    lift <- diff(range(c(sm$hi, sm$lo), finite = TRUE)) * 0.08
    ann[, y := hi + lift]
  }
  
  ggplot() +
    geom_col(data = sm, aes(x = bin, y = m),
             fill = "grey80", color = "black", width = 0.8) +
    geom_errorbar(data = sm, aes(x = bin, ymin = lo, ymax = hi), width = 0.16) +
    geom_point(data = bins, aes(x = bin, y = value, color = mouse),
               position = position_jitter(width = 0.06, height = 0),
               size = 2.2, alpha = 0.9) +
    geom_line(data = bins, aes(x = bin, y = value, group = mouse, color = mouse),
              linewidth = 1.4, alpha = 0.55) +
    { if (!is.null(ann)) geom_text(data = ann, aes(x = bin, y = y, label = letters), vjust = 0) } +
    scale_color_manual(values = mouse_pal, guide = "none") +
    labs(title = label, x = "Event Time (bins)", y = "Mean z (baseline-adj)") +
    theme_classic(base_size = 12) +
    theme(plot.title  = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

run_block_parametric <- function(dt, label, range, stem, bin_width = 5) {
  bins <- bin_summarize(dt, bin_width = bin_width, stat = "mean", range = range)
  an   <- rm_anova(bins)
  post <- param_posthoc_allpairs(bins)
  cld  <- cld_from_posthoc(post)
  p    <- plot_bin_bars(bins, label, cld)
  list(bins = bins, anova = an, posthoc = post, cld = cld, plot = p, stem = stem)
}

#### Stats:          Run all bin widths × all four conditions ####
bin_widths <- c(2.5, 5, 10)
bin_sensitivity_results <- list()

for (bw in bin_widths) {
  stem_bw <- gsub("\\.", "p", as.character(bw))
  bin_sensitivity_results[[stem_bw]] <- list(
    Start_Aq = run_block_parametric(
      Start_Aq_adj,
      paste0("Acquisition — Running Initiation (", bw, "s bins)"),
      c(-15, 30), paste0("Start_Acquisition_", stem_bw, "s"),
      bin_width = bw),
    Start_Ma = run_block_parametric(
      Start_Ma_adj,
      paste0("Maintenance — Running Initiation (", bw, "s bins)"),
      c(-15, 30), paste0("Start_Maintenance_", stem_bw, "s"),
      bin_width = bw),
    Stop_Aq = run_block_parametric(
      Stop_Aq_adj,
      paste0("Acquisition — Running Termination (", bw, "s bins)"),
      c(-30, 15), paste0("Stop_Acquisition_", stem_bw, "s"),
      bin_width = bw),
    Stop_Ma = run_block_parametric(
      Stop_Ma_adj,
      paste0("Maintenance — Running Termination (", bw, "s bins)"),
      c(-30, 15), paste0("Stop_Maintenance_", stem_bw, "s"),
      bin_width = bw)
  )
}

#### Plot:           SF7_TrendsRMANOVA ####
for (bw_res in bin_sensitivity_results) {
  print(bw_res$Start_Aq$plot); print(bw_res$Start_Ma$plot)
  print(bw_res$Stop_Aq$plot);  print(bw_res$Stop_Ma$plot)
}

#### Save:           SF7_TrendsRMANOVA ####
fig_dir <- .fig_dir("SF7_TrendsRMANOVA")

.write_block <- function(res) {
  fwrite(res$bins,    file.path(fig_dir, paste0(res$stem, "_BinsPerMouse.csv")))
  fwrite(res$posthoc, file.path(fig_dir, paste0(res$stem, "_Posthoc_PairedT_Raw.csv")))
  fwrite(res$cld,     file.path(fig_dir, paste0(res$stem, "_CLD_Raw.csv")))
  capture.output(print(res$anova$nice),
                 file = file.path(fig_dir, paste0(res$stem, "_RMANOVA_GG.txt")))
  ggsave(file.path(fig_dir, paste0(res$stem, ".tiff")),
         res$plot, width = 7.2, height = 4.8, units = "in",
         dpi = 600, device = "tiff", compression = "lzw", bg = "white")
  ggsave(file.path(fig_dir, paste0(res$stem, ".svg")),
         res$plot, width = 7.2, height = 4.8, units = "in",
         device = svglite::svglite, bg = "white")
}
for (bw_res in bin_sensitivity_results) {
  lapply(list(bw_res$Start_Aq, bw_res$Start_Ma, bw_res$Stop_Aq, bw_res$Stop_Ma),
         .write_block)
}
##################################################################################################################
#### Figure 7a - ACh event contrast threshold sensitivity ####
#### Transform:     Per-mouse Cohen's dz at each threshold (session-level aggregation) ####

thresholds <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

behavior_labs <- c(
  "Run - None"     = "Run On Wheel",
  "Stumble - None" = "Stumble On Wheel",
  "Move - None"    = "Move Under Wheel"
)

mouse_thresh_results <- rbindlist(lapply(thresholds, function(thr) {
  cat("Computing per-mouse dz at threshold:", thr, "\n")
  
  dt_thr <- as.data.table(phot_data)
  dt_thr[, Ratio_z_sess := (Ratio - mean(Ratio, na.rm = TRUE)) / sd(Ratio, na.rm = TRUE),
         by = Series]
  dt_thr[, Mouse := factor(sub("^([0-9]{4}).*$", "\\1", Series))]
  
  dt_thr[, `:=`(
    RunFlag  = as.integer(is.finite(Run_on_wheel_probability)     & Run_on_wheel_probability     > thr),
    StumFlag = as.integer(is.finite(Stumble_on_wheel_probability) & Stumble_on_wheel_probability > thr),
    MoveFlag = as.integer(is.finite(Move_under_wheel_probability) & Move_under_wheel_probability > thr)
  )]
  dt_thr[, state := fifelse(RunFlag  == 1, "Run",
                            fifelse(StumFlag == 1, "Stumble",
                                    fifelse(MoveFlag == 1, "Move",
                                            fifelse(On_wheel == 0, "None", NA_character_))))]
  dt_thr <- dt_thr[!is.na(state)]
  
  rbindlist(lapply(c("Run","Stumble","Move"), function(beh) {
    # session-level means for behavior and None
    sess_means <- dt_thr[state %in% c(beh, "None"), .(
      beh_mean  = mean(Ratio_z_sess[state == beh],    na.rm = TRUE),
      none_mean = mean(Ratio_z_sess[state == "None"], na.rm = TRUE)
    ), by = .(Mouse, Series)]
    sess_means <- sess_means[is.finite(beh_mean) & is.finite(none_mean)]
    if (!nrow(sess_means)) return(NULL)
    
    # mouse-level means
    mouse_means <- sess_means[, .(beh_mean  = mean(beh_mean,  na.rm = TRUE),
                                  none_mean = mean(none_mean, na.rm = TRUE)),
                              by = Mouse]
    mouse_means[, diff := beh_mean - none_mean]
    if (nrow(mouse_means) < 2) return(NULL)
    
    # Cohen's dz across mice
    sd_diff <- sd(mouse_means$diff)
    mouse_means[, `:=`(threshold = thr,
                       behavior  = beh,
                       dz        = diff / sd_diff)]
    mouse_means[, .(Mouse, threshold, behavior,
                    beh_mean, none_mean, diff, dz)]
  }), fill = TRUE)
}), fill = TRUE)

mouse_thresh_results[, contrast := paste0(behavior, " - None")]

#### Transform:     Pooled summary across mice ####
pooled_thresh <- mouse_thresh_results[, {
  n     <- .N
  m     <- mean(dz, na.rm = TRUE)
  s     <- sd(dz,   na.rm = TRUE)
  se    <- s / sqrt(n)
  tcrit <- qt(0.975, df = n - 1)
  .(mean_dz  = m, se_dz = se,
    ci_lo = m - tcrit * se, ci_hi = m + tcrit * se,
    n_mice = n)
}, by = .(threshold, contrast, behavior)]

pooled_thresh[,        behavior_lab := behavior_labs[contrast]]
mouse_thresh_results[, behavior_lab := behavior_labs[contrast]]

# rename Mouse to label for plotting
mouse_thresh_results[, Mouse := mouse_lookup$label[match(as.character(Mouse), mouse_lookup$series_id)]]
mouse_thresh_results[, Mouse := factor(Mouse, levels = mouse_lookup$label)]

#### Print:          ####
thresh_table2 <- pooled_thresh[order(contrast, threshold), .(
  threshold, contrast,
  mean_dz = round(mean_dz, 3),
  ci_lo   = round(ci_lo,   3),
  ci_hi   = round(ci_hi,   3),
  n_mice
)]
print(thresh_table2)

#### Plot:           ####
y_range <- range(c(pooled_thresh$ci_lo, pooled_thresh$ci_hi,
                   mouse_thresh_results$dz), na.rm = TRUE)
y_pad   <- diff(y_range) * 0.1
ylims   <- c(y_range[1] - y_pad, y_range[2] + y_pad)

p_thresh_dz <- ggplot() +
  annotate("rect", xmin = 0.775, xmax = 0.825,
           ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_point(data = mouse_thresh_results[!is.na(dz)],
             aes(x = threshold, y = dz, color = Mouse),
             size = 2.2, alpha = 0.85,
             position = position_jitter(width = 0.008, height = 0)) +
  geom_ribbon(data = pooled_thresh[!is.na(mean_dz)],
              aes(x = threshold, ymin = ci_lo, ymax = ci_hi, group = 1),
              fill = "grey40", alpha = 0.25) +
  geom_line(data = pooled_thresh[!is.na(mean_dz)],
            aes(x = threshold, y = mean_dz, group = 1),
            color = "black", linewidth = 1.4) +
  geom_point(data = pooled_thresh[!is.na(mean_dz)],
             aes(x = threshold, y = mean_dz),
             color = "black", size = 3) +
  facet_wrap(~ behavior_lab, ncol = 3) +
  scale_color_manual(values = mouse_pal, name = "Mouse") +
  scale_x_continuous(breaks = thresholds) +
  coord_cartesian(ylim = ylims) +
  labs(title    = "Threshold sensitivity — Cohen's dz by classifier probability threshold",
       subtitle = "Mouse-level session means; None = off-wheel frames only",
       x        = "Classifier probability threshold",
       y        = "Cohen's dz (behavior vs None)") +
  theme_classic(base_size = 12) +
  theme(plot.title       = element_text(hjust = 0.5),
        plot.subtitle    = element_text(hjust = 0.5),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "white", color = NA))
print(p_thresh_dz)

#### Save:           F7a_EventContrasts ####
fig_dir <- .fig_dir("F7a_EventContrasts")

saveRDS(mouse_thresh_results, file.path(fig_dir, "mouse_thresh_results.rds"))
saveRDS(pooled_thresh,        file.path(fig_dir, "pooled_thresh.rds"))
fwrite(thresh_table2,         file.path(fig_dir, "ThresholdSensitivity_CohendZ_MouseLevel.csv"))

ggsave(file.path(fig_dir, "Fig7a_EventContrasts.tiff"),
       p_thresh_dz, width = 10.0, height = 4.5, units = "in",
       dpi = 600, device = "tiff", compression = "lzw", bg = "white")
ggsave(file.path(fig_dir, "Fig7a_EventContrasts.svg"),
       p_thresh_dz, width = 10.0, height = 4.5, units = "in",
       device = svglite::svglite, bg = "white")

p_thresh_bare <- p_thresh_dz + theme_void()
ggsave(file.path(fig_dir, "Fig7a_EventContrasts_bare.svg"),
       p_thresh_bare, width = 10.0, height = 4.5, units = "in", bg = "white")
###########################################################################################################



#### Figure 7b - Pooled correlations of ACh and kinematics ####
#### Transform:      Session z-score and mouse factor ####

dt <- as.data.table(CompleteFiPhaDataset)
setnames(dt, "(time)", "time", skip_absent = TRUE)
dt[, Ratio_z_sess := (Ratio - mean(Ratio, na.rm = TRUE)) / sd(Ratio, na.rm = TRUE),
   by = Series]
dt[, Mouse := factor(sub("^([0-9]{4}).*$", "\\1", Series))]

#### Transform:      Define variables, trims, and off-wheel/all-sample splits ####
kin_vars <- c("Center_movement","Mouse_area","Mouse_length","Mouse_width","Mouse_rotation")

#Set cutpoints to remove implausible data
rngs <- list(Center_movement = c(0,   5),
             Mouse_area      = c(50,  200),
             Mouse_length    = c(15,  90),
             Mouse_width     = c(15,  30),
             Mouse_rotation  = c(0,   100))  # absolute scale as rotation data is -100 to 100 for counterwise and clockwise directions

prob_var <- "Run_on_wheel_probability"

off_all <- rbindlist(lapply(kin_vars, \(v) {
  dd <- dt[On_wheel == 0 & is.finite(Ratio_z_sess) & is.finite(get(v)),
           .(Series, Mouse, Ratio_z_sess, x = get(v))]
  if (v == "Mouse_rotation")  dd[, x := abs(x)]
  rng <- rngs[[v]]
  if (!is.null(rng)) dd <- dd[x >= rng[1] & x <= rng[2]]
  dd[, var := v][]
}), use.names = TRUE)

prob_all <- dt[is.finite(Ratio_z_sess) & is.finite(get(prob_var)),
               .(Series, Mouse, Ratio_z_sess, x = get(prob_var))][, var := prob_var][]

#### Transform:      Per-session ρ → per-mouse Fisher-z → across-mouse pooled ####
fisher_z     <- function(r) 0.5 * log((1 + r) / (1 - r))
fisher_z_inv <- function(z) (exp(2*z) - 1) / (exp(2*z) + 1)
corr_by_session <- function(D) suppressWarnings(cor(D$Ratio_z_sess, D$x, method = "spearman"))

per_sess <- rbindlist(list(
  off_all[,  .(rho = corr_by_session(.SD)), by = .(var, Mouse, Series)],
  prob_all[, .(rho = corr_by_session(.SD)), by = .(var, Mouse, Series)]
), use.names = TRUE)

per_mouse <- per_sess[, {
  z    <- fisher_z(rho)
  zbar <- mean(z[is.finite(z)], na.rm = TRUE)
  .(zbar = zbar, n_sessions = sum(is.finite(z)), rho_mean = fisher_z_inv(zbar))
}, by = .(var, Mouse)]

pooled <- per_mouse[is.finite(zbar), {
  n_m    <- .N; zbar_m <- mean(zbar, na.rm = TRUE)
  se_z   <- if (n_m > 1) sd(zbar, na.rm = TRUE) / sqrt(n_m) else NA_real_
  tcrit  <- if (n_m > 1) qt(0.975, df = n_m - 1) else NA_real_
  ci     <- if (is.finite(se_z) && is.finite(tcrit)) tcrit * se_z else NA_real_
  .(rho_pooled = fisher_z_inv(zbar_m), n_mice = n_m,
    rho_lo = if (is.finite(ci)) fisher_z_inv(zbar_m - ci) else NA_real_,
    rho_hi = if (is.finite(ci)) fisher_z_inv(zbar_m + ci) else NA_real_)
}, by = var][order(rho_pooled)]

#### Plot:           F7b_KinematicCorrelations ####
nice_names <- c(
  Center_movement          = "Mouse_movement (OFF-wheel)",
  Mouse_length             = "Mouse_length (OFF-wheel)",
  Mouse_width              = "Mouse_width (OFF-wheel)",
  Mouse_area               = "Mouse_area (OFF-wheel)",
  Mouse_rotation           = "Mouse_rotation abs (OFF-wheel)",
  Run_on_wheel_probability = "Run_on_wheel_probability (ALL samples)"
)
pooled[, var_lab := nice_names[var]]

p_forest <- ggplot(pooled[is.finite(rho_pooled)],
                   aes(x = rho_pooled, y = reorder(var_lab, rho_pooled))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = rho_lo, xmax = rho_hi),
                 height = 0.18, linewidth = 0.8, color = "black") +
  geom_point(size = 4.0, color = "black", shape = 21,
             fill = "white", stroke = 1.2) +
  labs(title = NULL, x = "Spearman ρ (pooled across mice)", y = NULL) +
  theme_classic(base_size = 12)
print(p_forest)

#### Save:           F7b_KinematicCorrelations ####
fig_dir <- .fig_dir("F7b_KinematicCorrelations")
fwrite(per_sess,  file.path(fig_dir, "corr_per_session.csv"))
fwrite(per_mouse, file.path(fig_dir, "corr_per_mouse.csv"))
fwrite(pooled,    file.path(fig_dir, "corr_pooled_equal_mouse.csv"))
ggsave(file.path(fig_dir, "Fig7b_pooled_corr.tiff"),
       p_forest, width = 12, height = 5, units = "cm", 
       device = "tiff", compression = "lzw", bg = "white")
ggsave(file.path(fig_dir, "Fig7b_pooled_corr.svg"),
       p_forest, width = 12, height = 5, units = "cm", bg = "white")

p_bare <- p_forest + theme(axis.title  = element_blank(), axis.text   = element_blank(),)
ggsave(file.path(fig_dir, "Fig7b_pooled_corr_bare.svg"), p_bare, width = 12, height = 8, units = "cm", bg = "white")

##################################################################################################################
#### Supplemental Figure 8 - Scatter plots of ACh vs kinematics ####
#### Transform:      Data thinning for scatterplots (per var, per session) ####

set.seed(1)
thin_per_series <- 400

off_thin  <- off_all[, if (.N <= thin_per_series) .SD else .SD[sample(.N, thin_per_series)],
                     by = .(var, Series)]
prob_thin <- prob_all[, if (.N <= thin_per_series) .SD else .SD[sample(.N, thin_per_series)],
                      by = Series]

#### Plot:           S8_KinematicScatterPlots ####
plot_one_off <- function(var_name) {
  dd      <- off_thin[var == var_name]
  limx    <- rngs[[var_name]]
  rho_row <- pooled[var == var_name]
  ann_txt <- if (nrow(rho_row)) {
    sprintf("ρ = %.2f (95%% CI %.2f, %.2f)",
            rho_row$rho_pooled, rho_row$rho_lo, rho_row$rho_hi)
  } else NA_character_
  x_label <- var_name
  
  p <- ggplot(dd, aes(x = x, y = Ratio_z_sess, color = Mouse)) +
    geom_point(alpha = 0.25, size = 0.7) +
    stat_ellipse(level = 0.90, linewidth = 0.6, alpha = 0.35) +
    geom_smooth(method = MASS::rlm, se = FALSE, linewidth = 1.0) +
    geom_smooth(aes(color = NULL), method = MASS::rlm, se = FALSE,
                linewidth = 1.4, color = "black") +
    coord_cartesian(xlim = limx) +
    scale_color_manual(values = mouse_pal_id, drop = FALSE) +
    labs(x = x_label, y = "z-score") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")
  if (!is.na(ann_txt))
    p <- p + annotate("text", x = Inf, y = Inf, label = ann_txt,
                      hjust = 1.05, vjust = 1.5, size = 3.5, fontface = "plain")
  p
}

plot_prob <- function() {
  dd      <- prob_thin
  rho_row <- pooled[var == "Run_on_wheel_probability"]
  ann_txt <- if (nrow(rho_row)) {
    sprintf("ρ = %.2f (95%% CI %.2f, %.2f)",
            rho_row$rho_pooled, rho_row$rho_lo, rho_row$rho_hi)
  } else NA_character_
  
  p <- ggplot(dd, aes(x = x, y = Ratio_z_sess, color = Mouse)) +
    geom_point(alpha = 0.25, size = 0.7) +
    stat_ellipse(level = 0.90, linewidth = 0.6, alpha = 0.35) +
    geom_smooth(method = MASS::rlm, se = FALSE, linewidth = 1.0) +
    geom_smooth(aes(color = NULL), method = MASS::rlm, se = FALSE,
                linewidth = 1.4, color = "black") +
    coord_cartesian(xlim = c(0, 1)) +
    scale_color_manual(values = mouse_pal_id, drop = FALSE) +
    labs(x = "Run_on_wheel_probability", y = "z-score") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")
  if (!is.na(ann_txt))
    p <- p + annotate("text", x = Inf, y = Inf, label = ann_txt,
                      hjust = 1.05, vjust = 1.5, size = 3.5, fontface = "plain")
  p
}

p_list <- lapply(kin_vars, \(v) { p <- plot_one_off(v); print(p); p })
pprob  <- plot_prob(); print(pprob)

#### Save:           S8_KinematicScatterPlots ####
fig_dir <- .fig_dir("S8_KinematicScatterPlots")
for (v in kin_vars) {
  p <- plot_one_off(v)
  ggsave(file.path(fig_dir, paste0("S8_scatter_", v, ".tiff")),
         p, width = 6.2, height = 4.3, units = "in", dpi = 600,
         device = "tiff", compression = "lzw", bg = "white")
  ggsave(file.path(fig_dir, paste0("S8_scatter_", v, ".svg")),
         p, width = 6.2, height = 4.3, units = "in", bg = "white")
}
ggsave(file.path(fig_dir, "S8_scatter_Run_on_wheel_probability.tiff"),
       pprob, width = 6.2, height = 4.3, units = "in", dpi = 600,
       device = "tiff", compression = "lzw", bg = "white")
ggsave(file.path(fig_dir, "S8_scatter_Run_on_wheel_probability.svg"),
       pprob, width = 6.2, height = 4.3, units = "in", bg = "white")
##################################################################################################################
#### Session info - for reproducibility ####
capture.output(sessionInfo(),
               file = file.path(out_root, "sessionInfo.txt"))
cat("\nAnalysis complete. Outputs in:\n", out_root, "\n")