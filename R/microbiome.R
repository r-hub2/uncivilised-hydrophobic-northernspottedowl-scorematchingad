#' 16s Microbiome Data for Soil-Transmitted Helminths
#'
#' The `microbiome` data contains paired DNA samples from before treatment and 21 months after treatment for helminth infections \insertCite{martin2019mi}{scorematchingad}.
#' This data was analysed by \insertCite{martin2019mi;textual}{scorematchingad} and a further subset was studied by \insertCite{scealy2023sc;textual}{scorematchingad}.
#' The data are from a study into the effect of helminth infections on the course of malaria infections (ImmunoSPIN-Malaria) in the Nangapanda subdistrict, Indonesia \insertCite{wiria2010do}{scorematchingad}.
#' As part of the study, some participants were given 400mg of albendazole every three months for 1.5 years,
#' remaining participants were given a placebo
#'  \insertCite{wiria2010do}{scorematchingad}.
#' 
#' @details 
#' The measurements in the data come from stool samples before and after treatment.
#' Gut microbiome prevalence was measured using 16s rRNA 454 sequencing \insertCite{martin2019mi}{scorematchingad}.
#' Helminth infections were detected by PCR or microscopy \insertCite{martin2019mi}{scorematchingad}.
#' 
#' The subset studied by  \insertCite{scealy2023sc;textual}{scorematchingad} contained only the measurements from before treatment, and only those individuals with a helminth infection.
#' These measurements can be obtained by running
#' ```
#' microbiome[(microbiome$Year == 2008) & microbiome$Helminth, ]
#' ```
#' Two further individuals (`IndividualID` of `2079` and `2280`) were deemed outliers by \insertCite{scealy2023sc;textual}{scorematchingad}.
#'
#' # Modifications from the Source
#' The `microbiome` data was created from the file `S1_Table.xlsx` hosted on [Nematode.net](http://nematode.net/Data/environmental_interaction/S1_Table.xlsx) using the below code.
#'
#' ```
#' microbiome <- readxl::read_excel("S1_Table.xlsx",
#'   range = "A3:AE303") #avoids the genus data, keeping - only phyla
#' metacolnames <- readxl::read_excel("S1_Table.xlsx",
#'   range = "A2:J2", 
#'   col_names = FALSE)
#' colnames(microbiome)[1:ncol(metacolnames)] <- metacolnames[1, ]
#' colnames(microbiome)[2] <- "Year"
#' microbiome[, 11] <- (microbiome$ct_Al <= 30) | (microbiome$ct_Na <= 30) |
#'   (microbiome$ct_Ad <= 30) | (microbiome$ct_St <= 30) |
#'   (microbiome$micr_Tt == 1)
#' colnames(microbiome)[11] <- "Helminth"
#' microbiome <- microbiome |>
#'   dplyr::mutate(across(c(1,2,3,12:31), as.integer)) |>
#'   dplyr::mutate(micr_Tt = as.logical(micr_Tt),
#'                 Treatment = as.logical(Treatment)) |>
#'   dplyr::rename(IndividualID = `Individual ID`)
#' microbiome <- as.data.frame(microbiome)
#save(microbiome, file = "../data/microbiome.rda")
#' ```
#' @format 
#' A dataframe with 300 rows (two rows per individual) and 31 columns:
#' \describe{
#'   \item{IndividualID}{An integer uniquely specifying the individual.}
#'   \item{Year}{The collection year for the sample. `2008` for before treatment. `2010` for after treatment.}
#'   \item{Sex}{`1` if female, `0` otherwise.}
#'   \item{Treatment}{`TRUE` if individual given 400mg of albendazole every three months for 1.5 years, `FALSE` otherwise.}
#'   \item{Age}{Age at first sample.}
#'   \item{ct_Al}{A Helminth measurement: The qPCR cycle threshold (CT) for *Ascaris lumbricoides* (large roundworm). *Ascaris lumbricoides* can be considered present if the value is 30 or less.}
#'   \item{ct_Na}{A Helminth measurement: The qPCR cycle threshold (CT) for *Necator americanus* (a hookworm). *Necator americanus* can be considered present if the value is 30 or less.}
#'   \item{ct_Ad}{A Helminth measurement: The qPCR cycle threshold (CT) for *Ancylostoma duodenale* (a hookworm). *Ancylostoma duodenale* can be considered present if the value is 30 or less.}
#'   \item{micr_Tt}{A Helminth measurement: The presence of *Trichuris trichiura* as determined by microscopy. A value of `TRUE` means *Trichuris trichiura* was detected.}
#'   \item{Helminth}{A Helminth measurement: If any of the above helminths were detected then `TRUE`, otherwise `FALSE`.}
#'   \item{Remaining columns}{Count prevalence of 18 bacterial phyla and 2 unclassified columns.}
#' }
#' 
#'
#' @source [S1_Table.xlsx](http://nematode.net/Data/environmental_interaction/S1_Table.xlsx) from [http://nematode.net](http://nematode.net/NN3_frontpage.cgi?navbar_selection=environmental_interaction&subnav_selection=Indonesia_16s). 
#' [S1_Table.xlsx](http://nematode.net/Data/environmental_interaction/S1_Table.xlsx) was created by Dr. Bruce A Rosa for \insertCite{martin2019mi;textual}{scorematchingad}. Permission to share this data was obtained from Dr. Bruce Rosa and Dr. Ivonne Martin.
#' @references
#' \insertAllCited{}
"microbiome"
