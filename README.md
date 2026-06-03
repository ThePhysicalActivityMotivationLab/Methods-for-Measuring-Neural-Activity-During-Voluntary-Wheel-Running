# Methods for Measuring Neural Activity During Voluntary Wheel Running

Companion data and analysis code for:

> Letsinger et al. "Methods for Measuring Neural Activity During Voluntary Wheel Running." *Journal of Neuroscience Methods* 2026.

---

## Repository Contents

| File | Description |
|------|-------------|
| `wheel_running_neural_activity_analysis.R` | Full analysis script reproducing all figures and statistics in the manuscript |
| `[RawFiPhaData.rds](https://osf.io/5jry9/files/uc4mn)` | Complete session-level dataset exported from FiPhA (all timepoints, all sessions) |
| `[1120_RVDG_Events.rds](https://osf.io/5jry9/files/s3dhe)` | FiPhA event-level export for Mouse 1 (right ventral dentate gyrus) |
| `[1123_RVDG_Events.rds](https://osf.io/5jry9/files/rbjxk)` | FiPhA event-level export for Mouse 2 (right ventral dentate gyrus) |
| `[1124_RVDG_Events.rds](https://osf.io/5jry9/files/k26c4)` | FiPhA event-level export for Mouse 3 (right ventral dentate gyrus) |
| `BORIS_annotations/` | Example annotated behavioral video and BORIS annotation file for classifier training |

---

## Pipeline Overview

This repository covers the R-based analysis stage of a three-step pipeline:

1. **DeepLabCut v3** — pose estimation of anatomical landmarks from behavioral video
2. **SimBA v3.8** — supervised behavioral classification using pose-derived features and user-defined ROIs
3. **FiPhA** — spectral unmixing, detrending, and event-based export of fiber photometry data
4. **R** (this repository) — peri-event signal processing, statistics, and figure generation

For step-by-step documentation of each tool:
- DeepLabCut: https://deeplabcut.github.io/DeepLabCut
- BORIS: https://boris.readthedocs.io
- SimBA: https://github.com/sgoldenlab/simba
- FiPhA: https://mfbridge.github.io/FiPhA

---

## Data Structure

### Nomenclature

| Term | Definition |
|------|------------|
| Event | A single wheel running bout plus 15 seconds before initiation and 15 seconds after termination |
| Series | A single 60-minute recording session containing multiple running events |
| Phase | Acquisition (week 12, first week of wheel exposure) or Maintenance (week 15, after one month of daily running) |

### FiPhA Export Format

The per-mouse RDS files (`1120_RVDG_Events.rds`, etc.) contain event-level data exported from FiPhA. Each file is a named list with the following structure:

```
$events
  $<filter_name>
    $<series_name>
      [[event_index]]  →  data.frame of timepoints for that event
```

Series names follow the convention `<MouseID>_<Region>_<Phase><Session>` (e.g., `1120_RVDG_Aq1` = Mouse 1120, right ventral dentate gyrus, acquisition session 1).

`RawFiPhaData.rds` contains the complete session-level dataframe with all timepoints across all sessions, used for off-wheel kinematic analyses.

---

## Software Requirements

| Software | Version |
|----------|---------|
| R | 4.5.1 |
| data.table | 1.17.8 |
| ggplot2 | 4.0.2 |
| nlme | 3.1-168 |
| emmeans | 2.0.3 |
| afex | 1.5-1 |
| lme4 | 2.0-1 |
| multcompView | 0.1-11 |
| Hmisc | 5.2-5 |
| MASS | 7.3-65 |
| ggridges | 0.5.7 |
| svglite | 2.2.2 |
| scico | 1.5.0 |
| cowplot | 1.2.0 |

Install all required packages:

```r
install.packages(c(
  "data.table", "ggplot2", "scico", "svglite", "ggridges",
  "nlme", "emmeans", "Hmisc", "MASS", "multcompView",
  "afex", "lme4", "cowplot"
))
```

---

## Usage

### Running the analysis on the provided data

1. Clone or download this repository
2. Open `wheel_running_neural_activity_analysis.R`
3. On line 60, set `dir_path` to the directory containing the `.rds` files:
   ```r
   dir_path <- "path/to/your/local/repository"
   ```
4. Run the script top to bottom

All figures and statistical outputs will save to an `Analysis/` folder created inside `dir_path`.

### Adapting to new data

1. Export event-level `.rds` files from FiPhA using the same event window parameters (15 seconds before initiation, 15 seconds after termination)
2. Export the complete session-level dataset from FiPhA as `RawFiPhaData.rds`
3. Place all `.rds` files in a single directory and set `dir_path` accordingly
4. Edit the `mouse_lookup` table at lines 38-42 to match your mouse IDs and labels:
   ```r
   mouse_lookup <- data.table(
     series_id = c("XXXX", "XXXX"),   # numeric prefix of your FiPhA series names
     label     = c("Mouse1", "Mouse2"),
     color     = c("#hex1",  "#hex2")
   )
   ```
5. Run the script top to bottom

### Important notes for adaptation

- Framerate, contrast, brightness, color spectrum, camera distance, mouse coat color, attached recording apparatus, and wheel orientation must be matched to maintain DeepLabCut and SimBA performance
- The classifier must be retrained on study-specific data; at least one annotated video per mouse from a minimum of five mice is recommended
- See the manuscript Discussion for full implementation guidance

---

## Figures Produced

| Figure | Script Section |
|--------|---------------|
| Figure 5a | Bout length vs. median z-score correlations |
| Figure 5b,c | Event length retention curves |
| Figure 6 | Peri-event acetylcholine traces by phase |
| Figure 7a | Cohen's dz threshold sensitivity analysis |
| Figure 7b | Pooled Spearman correlations with off-wheel kinematics |
| Supplemental Figure 3 | Event length and amplitude distributions |
| Supplemental Figure 4 | Per-mouse peri-event traces |
| Supplemental Figure 5 | Baseline window slope analysis |
| Supplemental Figure 7 | RMANOVA across bin widths |
| Supplemental Figure 8 | Acetylcholine vs. kinematic scatterplots |

---

## Contact

Ayland Letsinger  
Department of Kinesiology and Health Education  
University of Texas at Austin  
ayland.letsinger@austin.utexas.edu

Funding sources: This research was supported [in part] by the Intramural Research Program of the National Institutes of Health (NIH; Z1ES90998), a NIH Pathway to Independence Award from The National Institute of Drug Abuse (K99DA058974-01A1), and contract with Social & Scientific Systems, a DLH Holdings Corp. Company (GS-00F-173CA-75N96022F00055). The contributions of the NIH author(s) are considered Works of the United States Government. The findings and conclusions presented in this paper are those of the author(s) and do not necessarily reflect the views of the NIH or the U.S. Department of Health and Human Services.
