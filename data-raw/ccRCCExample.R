## code to prepare `ccRCCExample` dataset goes here
ccRCCExample = readRDS("data-raw/RCC7_16.rds")
usethis::use_data(ccRCCExample, overwrite = TRUE)
