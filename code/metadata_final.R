## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## metadata for sharks
## GP Sept 2026
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list =ls())
library(dplyr)
library(lubridate)

# output directory for metadata products
dir.create(
  "../results/metadata",
  showWarnings = FALSE,
  recursive = TRUE)

#read in data
results_df <- read.csv("../results/results_annual_kde_coa.csv")
capture_df <- read.csv("../data/raw/acoustic_metadata/lemon_summary_read_updated.csv")
gap_df <- read.csv("../results/metadata/shark_detection_gaps.csv")

#--- aggregate results_df to one row per shark ---
metadata_summary <- results_df %>%
  mutate(max_date = ymd_hms(max_date), min_date = ymd_hms(min_date)) %>%
  group_by(shark_ID) %>%
  summarise(
    sex = tolower(first(sex)),
    
    final_date  = max(max_date, na.rm = TRUE), # latest detection date across all years
    
    n_tracking_days_total = sum(n_tracking_days),      #total days detected, summed across years
    total_detections = sum(num_detections),
    years_with_detections = n(),                       #number of rows = number of years with any data
    
    final_size = size[which.max(max_date)],             #size estimate for that final year
    .groups = "drop")

#--- pull tagging info from capture_df ---
capture_clean <- capture_df %>%
  transmute(
    shark_ID = FishID,
    sex = tolower(Sex),
    tagging_date = dmy_hm(Date.released),
    tagging_size = pcl,
    tag_type = Tag_type,
    estimated_tag_life_days = Est_life)

#--- join and compute derived columns ---
metadata_summary <- metadata_summary %>%
  left_join(capture_clean %>% dplyr::select(shark_ID, tagging_date, tagging_size, tag_type,
                                            estimated_tag_life_days), by = "shark_ID") %>%
  mutate(
    track_length = as.integer(as.Date(final_date) - as.Date(tagging_date)+1),
    tracking_ratio = round(n_tracking_days_total / track_length, 3),
    average_days_per_year = as.integer(round(n_tracking_days_total / years_with_detections)),
    average_detections_per_year = as.integer(round(total_detections / years_with_detections))) %>%
  dplyr::select(
    shark_ID, sex,
    tagging_date, tagging_size,
    tag_type, estimated_tag_life_days,
    final_date, final_size,
    n_tracking_days_total, total_detections,
    track_length, tracking_ratio,
    years_with_detections, average_days_per_year, average_detections_per_year)

#add maturity at tagging date and last detection
metadata_summary <- metadata_summary %>%
  mutate(
    maturity_at_tagging = case_when(
      is.na(tagging_size) ~ NA_character_,
      tagging_size >= 166 ~ "mature",
      TRUE ~ "immature"),
    maturity_at_last_detection = case_when(
      is.na(final_size) ~ NA_character_,
      final_size >= 166 ~ "mature",
      TRUE ~ "immature"))

metadata_summary <- metadata_summary %>%
  left_join(
    gap_df %>% rename(shark_ID = FishID),
    by = "shark_ID")

metadata_summary <- metadata_summary %>%
  mutate(
    expected_tag_expiry =
      as.Date(tagging_date) + estimated_tag_life_days,
    
    days_to_expected_expiry =
      as.integer(expected_tag_expiry - as.Date(final_date)))

write.csv(
  metadata_summary,
  "../results/metadata/shark_metadata_full.csv",
  row.names = FALSE)

saveRDS(
  metadata_summary,
  "../results/metadata/shark_metadata_full.rds")

metadata_final <- metadata_summary %>%
  transmute(
    `Shark ID` = shark_ID, Sex = toupper(sex),
    `Tag type` = tag_type,
    `Days before expected tag expiry` = days_to_expected_expiry,
    `Tagging date` = as.Date(tagging_date),
    `Size at tagging (cm PCL)` = round(tagging_size),
    `Last detection` = as.Date(final_date),
    `Estimated size at last detection (PCL, cm)` = round(final_size),
    `Days detected` = n_tracking_days_total,
    `Tracking span (days)` = track_length,
    `Days detected (%)` = round(tracking_ratio * 100, 1),
    `Years detected` = years_with_detections,
    `Life stage at tagging` = maturity_at_tagging,
    `Life stage at last detection` = maturity_at_last_detection,
    `Maximum no detection gap (days)` = max_detection_gap_days)

write.csv(
  metadata_final,
  "../results/metadata/shark_metadata_final.csv",
  row.names = FALSE)

# ---------------- General results paragraph --------------------------- 
report <- metadata_summary %>%
  summarise(
    n_sharks = n(),
    n_females = sum(sex == "f"),
    n_males = sum(sex == "m"),
    total_detections = sum(total_detections, na.rm = TRUE),
    
    mean_tracking_span = mean(track_length, na.rm = TRUE),
    median_tracking_span = median(track_length, na.rm = TRUE),
    min_tracking_span = min(track_length, na.rm = TRUE),
    max_tracking_span = max(track_length, na.rm = TRUE),
    
    mean_tagging_size = mean(tagging_size, na.rm = TRUE),
    min_tagging_size = min(tagging_size, na.rm = TRUE),
    max_tagging_size = max(tagging_size, na.rm = TRUE),
    
    earliest_tagging = min(tagging_date, na.rm = TRUE),
    latest_detection = max(final_date, na.rm = TRUE),
    
    mean_days_detected_pct = mean(tracking_ratio * 100, na.rm = TRUE),
    median_days_detected_pct = median(tracking_ratio * 100, na.rm = TRUE),
    
    median_max_gap = median(max_detection_gap_days, na.rm = TRUE),
    max_gap = max(max_detection_gap_days, na.rm = TRUE))
report

metadata_summary %>%
  slice_max(track_length, n = 1) %>%
  dplyr::select(shark_ID, track_length)

metadata_summary %>%
  count(maturity_at_tagging)

metadata_summary %>%
  count(maturity_at_last_detection)

metadata_summary %>%
  count(maturity_at_tagging, maturity_at_last_detection)

range(metadata_summary$tagging_date, na.rm = TRUE)
range(metadata_summary$final_date, na.rm = TRUE)

#------------ Tag life and cessation of detections -----------------
metadata_tags <- metadata_summary %>%
  mutate(
    proportion_tag_life_observed =
      track_length / estimated_tag_life_days)

summary(metadata_tags$days_to_expected_expiry)

metadata_tags %>%
  dplyr::select(
    shark_ID,
    tag_type,
    estimated_tag_life_days,
    tagging_date,
    final_date,
    expected_tag_expiry,
    days_to_expected_expiry,
    proportion_tag_life_observed) %>%
  arrange(days_to_expected_expiry)

metadata_tags %>%
  summarise(
    median_days_before_expiry =
      median(days_to_expected_expiry, na.rm = TRUE),
    
    n_before_expiry =
      sum(days_to_expected_expiry > 0, na.rm = TRUE),
    
    n_after_expiry =
      sum(days_to_expected_expiry <= 0, na.rm = TRUE),
    
    median_prop_tag_life =
      median(proportion_tag_life_observed, na.rm = TRUE))

metadata_tags %>%
  group_by(tag_type) %>%
  summarise(
    n = n(),
    median_tag_life = median(estimated_tag_life_days, na.rm = TRUE),
    median_days_to_expiry = median(days_to_expected_expiry, na.rm = TRUE),
    median_prop_life_observed = median(proportion_tag_life_observed, na.rm = TRUE))

study_end <- as.Date("2021-12-31")

metadata_tags <- metadata_tags %>%
  mutate(
    days_before_study_end =
      as.integer(study_end - as.Date(final_date)))

# How many sharks were still being detected close to study end?
metadata_tags %>%
  summarise(
    within_30d  = sum(days_before_study_end <= 30),
    within_90d  = sum(days_before_study_end <= 90),
    within_180d = sum(days_before_study_end <= 180))

metadata_tags %>%
filter(days_before_study_end > 90) %>%
  summarise(
    n = n(),
    median_days_before_expiry =
      median(days_to_expected_expiry, na.rm = TRUE),
    n_near_expiry =
      sum(abs(days_to_expected_expiry) <= 90, na.rm = TRUE),
    n_well_before_expiry =
      sum(days_to_expected_expiry > 90, na.rm = TRUE),
    median_prop_tag_life_elapsed =
      median(proportion_tag_life_observed, na.rm = TRUE))
# 36 sharks for which detections ceased before the end of monitoring
# 32 of them > 90 days before tag expiry  # the four left are the really short tags
# battery exhaustion cannot explain most cessation of detections

metadata_tags %>%
  filter(days_before_study_end > 90) %>%
  count(maturity_at_last_detection) 
# among the 36 sharks, 17 were immature and 19 mature at last detection

metadata_tags %>%
  filter(days_before_study_end > 90) %>%
  summarise(
    mean_final_size = mean(final_size, na.rm = TRUE),
    median_final_size = median(final_size, na.rm = TRUE),
    min_final_size = min(final_size, na.rm = TRUE),
    max_final_size = max(final_size, na.rm = TRUE)) 
# final estimated sizes spanning 95.1–214 cm PCL and a median of 173 cm PCL