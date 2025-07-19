#!/usr/bin/env Rscript

# TODO:
# Use 1:5 random pick for SAVED & 3x1min for MV & QPAD
# Use c(0.5,1,1.5,2)^2*pi random pick for SAVED and 1.5 as maxdis for MV & QPAD (pi*1.5^2=7.07)

# test:
# Rscript --vanilla 02-bsims-simulations.R --dir _tmp/test --n 10
# Rscript --vanilla 02-bsims-simulations.R --dir _tmp/test --n 25 --all
# Rscript --vanilla 02-bsims-simulations.R --dir _tmp/test --n 25 --move --phi 0.5

# real runs:

# cd sqpad-paper/analysis

# Rscript --vanilla 02-bsims-simulations.R --D 0.5 --phi 0.25 --n 1000
# Rscript --vanilla 02-bsims-simulations.R --D 0.5 --phi 0.25 --jointdist --n 1000
# Rscript --vanilla 02-bsims-simulations.R --D 0.5 --phi 0.25 --all --n 1000
# Rscript --vanilla 02-bsims-simulations.R --D 0.5 --phi 0.25 --jointdist --move --disterr --n 1000
# Rscript --vanilla 02-bsims-simulations.R --D 0.5 --phi 0.25 --jointdist --move --doublecount --n 1000

# Rscript --vanilla 02-bsims-simulations.R --D 0.5 --phi 0.5 --n 1000
# Rscript --vanilla 02-bsims-simulations.R --D 0.5 --phi 0.5 --jointdist --n 1000
# Rscript --vanilla 02-bsims-simulations.R --D 0.5 --phi 0.5 --all --n 1000
# Rscript --vanilla 02-bsims-simulations.R --D 0.5 --phi 0.5 --jointdist --move --disterr --n 1000
# Rscript --vanilla 02-bsims-simulations.R --D 0.5 --phi 0.5 --jointdist --move --doublecount --n 1000

# Rscript --vanilla 02-bsims-simulations.R --D 1 --phi 0.25 --n 1000
# Rscript --vanilla 02-bsims-simulations.R --D 1 --phi 0.25 --jointdist --n 1000
# Rscript --vanilla 02-bsims-simulations.R --D 1 --phi 0.25 --all --n 1000
# Rscript --vanilla 02-bsims-simulations.R --D 1 --phi 0.25 --jointdist --move --disterr --n 1000
# Rscript --vanilla 02-bsims-simulations.R --D 1 --phi 0.25 --jointdist --move --doublecount --n 1000

# Rscript --vanilla 02-bsims-simulations.R --D 1 --phi 0.5 --n 1000
# Rscript --vanilla 02-bsims-simulations.R --D 1 --phi 0.5 --jointdist --n 1000
# Rscript --vanilla 02-bsims-simulations.R --D 1 --phi 0.5 --all --n 1000
# Rscript --vanilla 02-bsims-simulations.R --D 1 --phi 0.5 --jointdist --move --disterr --n 1000
# Rscript --vanilla 02-bsims-simulations.R --D 1 --phi 0.5 --jointdist --move --doublecount --n 1000

t0 <- proc.time()

library(bSims)
library(rconfig)
library(qs)
source("01-functions-helpers.R")

CONFIG <- rconfig::rconfig()
str(CONFIG)

DIR <- rconfig::value(CONFIG$dir, "_tmp/bsims")

ALL <- rconfig::value(CONFIG$all, FALSE)
if (ALL) {
    MOVE <- TRUE
    DISTERR <- TRUE
    DOUBLECOUNT <- TRUE
    JOINTDIST <- TRUE
} else {
    MOVE <- rconfig::value(CONFIG$move, FALSE)
    DISTERR <- rconfig::value(CONFIG$disterr, FALSE)
    DOUBLECOUNT <- rconfig::value(CONFIG$doublecount, FALSE)
    JOINTDIST <- rconfig::value(CONFIG$jointdist, FALSE)
}

n <- rconfig::value(CONFIG$n, 10000)

D <- rconfig::value(CONFIG$D, 1)
phi <- rconfig::value(CONFIG$phi, 0.25)
tau <- rconfig::value(CONFIG$tau, 1)

ncl <- NULL
# ncl <- 4
tbr <- 1:10
rbr <- c(0.25, 0.5, 1, 1.5, 2, Inf)

s <- list(
    density = D,
    duration = max(tbr)*2,
    vocal_rate = phi,
    tau = tau,
    tint = tbr,
    rint = rbr,
    abund_fun = NULL,
    move_rate = 0,
    movement = 0,
    condition = "event1",
    error = 0,
    bias = 1,
    perception = NULL
)
if (MOVE) {
    s$move_rate <- 0.5
    s$movement <- 0.25
}
if (DISTERR) {
    s$error <- 0
    s$bias <- 0.8
}
if (DOUBLECOUNT) {
    s$perception <- 1.2
}
if (JOINTDIST) {
    s$condition <- "det1"
}

str(s)

message("Run with settings ",
    paste(c("n", "D", "phi", "tau"), c(n, D, phi, tau), sep = "=", collapse=" "),
    " [", paste0(c("M", "E", "C", "J"), as.integer(c(MOVE, DISTERR, DOUBLECOUNT, JOINTDIST)), collapse="+"), "]")
b <- bsims_all(s)
o <- b$replicate(n, cl=ncl)

fn <- paste0("n-", n, "_D-", D, "_phi-", phi, "_tau-", tau,
    "_M-", as.integer(MOVE), "_E-", as.integer(DISTERR), "_C-", as.integer(DOUBLECOUNT), "_J-", as.integer(JOINTDIST))

f <- file.path(DIR, paste0("sqpad-bsims_", fn, ".qs"))
make_dir(f)
qs::qsave(o, f)

message("Finished simulations in ", round((proc.time()[3] - t0[3])/60, 1), " min with ",
    paste(c("n", "D", "phi", "tau"), c(n, D, phi, tau), sep = "=", collapse=" "),
    " [", paste0(c("M", "E", "C", "J"), as.integer(c(MOVE, DISTERR, DOUBLECOUNT, JOINTDIST)), collapse="+"), "]")
message("Results saved into file `", f, "`")
message("DONE")

quit(save = "no", status = 0)
