This repo gathers all the input files and scripts related to our study entitled "Leveraging of SARS-CoV-2 PCR cycle thresholds values (Ct) to forecast COVID-19 trends" (Yin *et al*. *submitted*): R script to prepare, analyse, and visualise the Ct data.

The Ct data are all gathered in the file "Data_LHUB-ULB_200521.csv". The R script "R_scripts_for_Ct_analyses.r" is divided into different sections allowing to perform the following analytical steps:

(1) estimating the median and mean Ct values through time

(2) co-plotting the mean Ct values on phase diagrams*

(3) plotting the median and mean Ct values through time

(*) the phase diagrams are generated according to Hens et al. (2021, Archives of Public Health)

In addition, we also gathered within the subdirectory "Ct_values_Analyses_by_age" the R script and files required to run specific analyses by age category.

Summary of the study: We assessed the usefulness of SARS-CoV-2 RT-PCR cycle thresholds (Ct) values trends produced by the LHUB-ULB (a consolidated microbiology laboratory located in Brussels, Belgium) for monitoring the epidemic’s dynamics at local and national levels and for improving forecasting models. SARS-CoV-2 RT-PCR Ct values produced from April 1, 2020, to May 15, 2021, were compared with national COVID-19 confirmed cases notifications according to their geographical and time distribution. These Ct values were evaluated against both a phase diagram predicting the number of COVID-19 patients requiring intensive care and an age-structured model estimating COVID-19 prevalence in Belgium. Over 155,811 RT-PCR performed, 12,799 were positive and 7,910 Ct values were available for analysis. The 14-day median Ct values were negatively correlated with the 14-day mean daily positive tests with a lag of 17 days. In addition, the 14-day mean daily positive tests in LHUB-ULB were strongly correlated with the 14-day mean confirmed cases in the Brussels-Capital and in Belgium with coinciding start, peak and end of the different waves of the epidemic. Ct values decreased concurrently with the forecasted phase-shifts of the diagram. Similarly, the evolution of 14-day median Ct values was negatively correlated with daily estimated prevalence for all age-classes. We provide preliminary evidence that trends of Ct values can help to both follow and predict the epidemic’s trajectory at local and national levels, underlining that consolidated microbiology laboratories can act as epidemic sensors as they gather data that are representative of the geographical area they serve.
