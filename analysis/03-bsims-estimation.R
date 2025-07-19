#!/usr/bin/env Rscript

# Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 1

t0 <- proc.time()
library(rconfig)
library(bSims)
library(detect)
# library(unmarked)
library(parallel)
library(qs)
source("00-functions-sqpad.R")
source("01-functions-helpers.R")

CONFIG <- rconfig::rconfig()
# str(CONFIG)

DIR <- rconfig::value(CONFIG[["dir"]], "_tmp/est_conv_mc")
B <- rconfig::value(CONFIG[["B"]], 200)
n <- rconfig::value(CONFIG[["n"]], 200)
set <- rconfig::value(CONFIG[["set"]], 1)
replace <- rconfig::value(CONFIG[["replace"]], FALSE)
notify <- rconfig::value(CONFIG[["notify"]], FALSE)

str(list(dir=DIR, B=B, n=n, set=set, replace=replace, notify=notify))


fl <- c(
    "sqpad-bsims_n-1000_D-0.5_phi-0.25_tau-1_M-0_E-0_C-0_J-1.qs", # 1 DONE
    "sqpad-bsims_n-1000_D-0.5_phi-0.5_tau-1_M-0_E-0_C-0_J-1.qs", # 2 DONE
    "sqpad-bsims_n-1000_D-1_phi-0.25_tau-1_M-0_E-0_C-0_J-1.qs", # 3 DONE
    "sqpad-bsims_n-1000_D-1_phi-0.5_tau-1_M-0_E-0_C-0_J-1.qs", # 4 DONE

    "sqpad-bsims_n-1000_D-0.5_phi-0.25_tau-1_M-1_E-1_C-1_J-1.qs", # 5 DONE
    "sqpad-bsims_n-1000_D-0.5_phi-0.5_tau-1_M-1_E-1_C-1_J-1.qs", # 6 DONE
    "sqpad-bsims_n-1000_D-1_phi-0.25_tau-1_M-1_E-1_C-1_J-1.qs", # 7 DONE
    "sqpad-bsims_n-1000_D-1_phi-0.5_tau-1_M-1_E-1_C-1_J-1.qs", # 8 DONE

    "sqpad-bsims_n-1000_D-0.5_phi-0.25_tau-1_M-1_E-1_C-0_J-1.qs", # 9 DONE
    "sqpad-bsims_n-1000_D-0.5_phi-0.5_tau-1_M-1_E-1_C-0_J-1.qs", # 10 DONE
    "sqpad-bsims_n-1000_D-1_phi-0.25_tau-1_M-1_E-1_C-0_J-1.qs", # 11 DONE
    "sqpad-bsims_n-1000_D-1_phi-0.5_tau-1_M-1_E-1_C-0_J-1.qs", # 12 DONE

    "sqpad-bsims_n-1000_D-0.5_phi-0.25_tau-1_M-1_E-0_C-1_J-1.qs", # 13 DONE
    "sqpad-bsims_n-1000_D-0.5_phi-0.5_tau-1_M-1_E-0_C-1_J-1.qs", # 14 DONE
    "sqpad-bsims_n-1000_D-1_phi-0.25_tau-1_M-1_E-0_C-1_J-1.qs", # 15 DONE
    "sqpad-bsims_n-1000_D-1_phi-0.5_tau-1_M-1_E-0_C-1_J-1.qs") # 16 DONE

fn <- fl[set]
sup <- gsub("\\.qs$", "", gsub("sqpad-bsims_n-1000_", "", fn))

s <- strsplit(strsplit(gsub("\\.qs$", "", fn), "_")[[1L]], "-")
s <- structure(as.numeric(sapply(s, "[[", 2)[-1]), names = sapply(s, "[[", 1)[-1])

f <- file.path("_tmp/bsims", fn)
all_sims <- qs::qread(f)

message("setup: ", set, ", n = ", n, ", file = ", fn)

res <- NULL
for (b in 1:B) {
    message("\n--->>> ", "setup: ", set, ", run: ", b, " / ", B, " [", sup, "]\n")
    r <- est_fun(n=n, all_sims=all_sims,
        replace = replace,
        verbose = TRUE)

    res <- cbind(
        res, 
        c(s[names(s) != "n"], r))
}
res <- data.frame(t(res))
res$source_file <- fn
res$run <- 1:B
rownames(res) <- NULL

f <- file.path(DIR, 
    paste0(gsub("\\.qs", "", gsub("sqpad-bsims_n-1000", "estimates_bsims", fn)), "_", 
        paste0(sample(letters,12), collapse=""),
        ".rds"))
make_dir(f)
saveRDS(res, f)

if (notify) {
    topic <- "a8m_bsims_alerts"
    details <- paste0(
        "set=", paste0(set, collapse = ","), 
        ", n=", n, ", B=", B,
        collapse = "")
    tm <- as.numeric(proc.time()[3] - t0[3])/60
    tm <- if (tm < 60) {
        paste0(round(tm, 1), " mins")
    } else {
        if (tm < 24*60) {
            paste0(round(tm/60, 1), " hrs")
        } else {
            paste0(tm %/% (60*24), " days & ", 
                round(tm %% (60*24) / 60), " hrs")
        }
    }
    msg <- paste("Finished job in", tm,
        "@", format(.POSIXct(Sys.time(), "America/Edmonton")), "\nJob details:", details)
    system2("curl", c(
        "-H", "\"X-Tags: white_check_mark\"", 
        # "-H", "\"X-Priority: 4\"", 
        "-d", sprintf("\"%s\"", msg), 
        sprintf("ntfy.sh/%s", topic)))
}

quit(save = "no", status = 0)
