## bSims simulations

fl <- c(list.files("_tmp/est_conv_mc", full.names = TRUE))
xb <- do.call(rbind, lapply(fl, readRDS))
arrow::write_parquet(xb, "estimates_bsims_mc.parquet")

## field data

SPP <- c("CHSP", "DEJU", "HETH", "MAWA", "MYWA", "OVEN", "RCKI", "REVI", "SWTH", "TEWA", "WIWR", "WTSP")
xf <- do.call(rbind, lapply(paste0("_tmp/paired_mc/estimates_field-data_", SPP, ".rds"), readRDS))
arrow::write_parquet(xf, "estimates_field_data_mc.parquet")
