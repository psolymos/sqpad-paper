## bSims simulations

fl <- c(list.files("_tmp/est_conv", full.names = TRUE),
    list.files("_tmp/est_n1000", full.names = TRUE))
# fl <- c(list.files("_tmp/est_conv", full.names = TRUE))
xb <- do.call(rbind, lapply(fl, readRDS))
arrow::write_parquet(xb, "estimates_bsims.parquet")

fl <- c(list.files("_tmp/est_conv_mc", full.names = TRUE))
xb <- do.call(rbind, lapply(fl, readRDS))
arrow::write_parquet(xb, "estimates_bsims_mc.parquet")

## field data

SPP <- c("CHSP", "DEJU", "HETH", "MAWA", "MYWA", "OVEN", "RCKI", "REVI", "SWTH", "TEWA", "WIWR", "WTSP")
xf <- do.call(rbind, lapply(paste0("_tmp/paired/estimates_field-data_", SPP, ".rds"), readRDS))
arrow::write_parquet(xf, "estimates_field_data.parquet")

xf <- do.call(rbind, lapply(paste0("_tmp/paired_mc/estimates_field-data_", SPP, ".rds"), readRDS))
arrow::write_parquet(xf, "estimates_field_data_mc.parquet")
