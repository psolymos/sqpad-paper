#!/usr/bin/env Rscript

# Rscript --vanilla 05-field-data.R --B 2 --id 2

# SPP <- c("CHSP", "DEJU", "HETH", "MAWA", "MYWA", "OVEN", "RCKI", "REVI", "SWTH", "TEWA", "WIWR", "WTSP")

# Rscript --vanilla 04-field-data.R --id 1 # CHSP
# Rscript --vanilla 04-field-data.R --id 2 # DEJU
# Rscript --vanilla 04-field-data.R --id 3 # HETH
# Rscript --vanilla 04-field-data.R --id 4 # MAWA
# Rscript --vanilla 04-field-data.R --id 5 # MYWA
# Rscript --vanilla 04-field-data.R --id 6 # OVEN
# Rscript --vanilla 04-field-data.R --id 7 # RCKI
# Rscript --vanilla 04-field-data.R --id 8 # REVI
# Rscript --vanilla 04-field-data.R --id 9 # SWTH
# Rscript --vanilla 04-field-data.R --id 10 # TEWA
# Rscript --vanilla 04-field-data.R --id 11 # WIWR
# Rscript --vanilla 04-field-data.R --id 12 # WTSP


library(rconfig)
library(detect)
library(paired) # remotes::install_github("borealbirds/paired")
library(mefa4)
source("00-functions-sqpad.R")
source("01-functions-helpers.R")

CONFIG <- rconfig::rconfig()
str(CONFIG)

DIR <- rconfig::value(CONFIG$dir, "_tmp/paired_mc")
spp_id <- rconfig::value(CONFIG$id, 1) # 1-12
B <- rconfig::value(CONFIG$B, 200)


sqpad_est <- function(spp, x, id, type = "full", det = "joint", K = NULL, hessian = FALSE, ...) {
    x$UniqueID <- as.character(x$UniqueID)
    xx <- do.call(rbind, lapply(1:length(id), function(i) {
        z <- x[x$UniqueID==id[i],]
        z$UniqueID <- paste0(z$UniqueID, "_", i)
        z
    }))
    xx$UniqueID <- as.factor(xx$UniqueID)
    xt <- Xtab(Count ~ UniqueID + Interval + DISTANCE, xx[xx$SPECIES == spp,])
    d <- nonDuplicated(xx, UniqueID, TRUE)
    n <- nrow(d)

    y11 <- xt[["0-49 m"]][,c("0-3 min")]
    y12 <- rowSums(xt[["0-49 m"]][,c("0-3 min", "3-5 min")])
    y13 <- rowSums(xt[["0-49 m"]][,c("0-3 min", "3-5 min", "5-10 min")])
    y21 <- y11 + xt[["50-100 m"]][,c("0-3 min")]
    y22 <- y12 + rowSums(xt[["50-100 m"]][,c("0-3 min", "3-5 min")])
    y23 <- y13 + rowSums(xt[["50-100 m"]][,c("0-3 min", "3-5 min", "5-10 min")])
    y31 <- y11 + y21 + xt[[">100 m"]][,c("0-3 min")]
    y32 <- y12 + y22 + rowSums(xt[[">100 m"]][,c("0-3 min", "3-5 min")])
    y33 <- y13 + y23 + rowSums(xt[[">100 m"]][,c("0-3 min", "3-5 min", "5-10 min")])

    # dis <- rep(c(0.5, 1, 4), each = 3 * n)
    # dur <- c(rep(c(3, 5, 10), each = n), rep(c(3, 5, 10), each = n), rep(c(3, 5, 10), each = n))
    # pid <- rep(1:n, 3*3)
    # y <- c(y11, y12, y13, y21, y22, y23, y31, y32, y33)

    dis <- rep(c(0.5, 1), each = 3 * n)
    dur <- c(rep(c(3, 5, 10), each = n), rep(c(3, 5, 10), each = n))
    pid <- rep(1:n, 2*3)
    y <- c(y11, y12, y13, y21, y22, y23)

    X <- model.matrix(~scale(Latitude), d)[pid,]
    Z <- model.matrix(~TSSR+I(TSSR^2), d)[pid,]

    mod <- sqpad.fit(Y=y, dis=dis, dur=dur, X=X, Z=Z, type = type, det = det, K = K, hessian = hessian, ...)
    o <- unname(mean(exp(X %*% coef(mod)[1:ncol(X)])))
    attr(o, "ll") <- mod$loglik
    o
}

qpad_est <- function(spp, x, id, truncate = FALSE) {
    d <- nonDuplicated(x, UniqueID, TRUE)
    n <- nrow(d)
    X <- model.matrix(~scale(Latitude), d)
    Z <- model.matrix(~TSSR+I(TSSR^2), d)
    rn <- rownames(d)
    if (missing(id)) {
        rn <- id
    } else {
        rn <- id
    }

    # distance
    # yydis <- as.matrix(Xtab(Count ~ UniqueID + DISTANCE, x[x$SPECIES == spp,]))[rn,c("0-49 m", "50-100 m", ">100 m")]
    yydis <- as.matrix(Xtab(Count ~ UniqueID + DISTANCE, x[x$SPECIES == spp,]))[rn,c("0-49 m", "50-100 m")]
    # dddis <- matrix(c(0.5, 1, Inf), nrow(yydis), 3, byrow=TRUE)
    dddis <- matrix(c(0.5, 1), nrow(yydis), 2, byrow=TRUE)
    mdis <- cmulti.fit(Y=yydis, D=dddis, type = "dis")

    # time
    yydur <- as.matrix(Xtab(Count ~ UniqueID + Interval,
        x[x$SPECIES == spp & x$DISTANCE %in% c("0-49 m", "50-100 m"),]))[rn,c("0-3 min", "3-5 min", "5-10 min")]
    dddur <- matrix(c(3, 5, 10), nrow(yydur), 3, byrow=TRUE)
    mdur <- cmulti.fit(Y=yydur, D=dddur, X=Z, type = "rem")

    # QPAD
    Y <- rowSums(yydis)
    TAU <- exp(coef(mdis))
    PHI <- exp(Z %*% coef(mdur))
    # Corr <- TAU^2 * pi * (1 - exp(-PHI * 10))
    Area <- rep(1^2 * pi, n) # 100m radius
    p <- 1 - exp(-PHI * 10)
    r <- 1
    q <- TAU^2 * (1 - exp(-r^2/TAU^2))/r^2
    Off <- log(Area * p * q)
    if (truncate) {
        mod <- glm(ifelse(Y>0,1,0) ~ lat, data.frame(lat=X[rn,2]), offset = Off, family=binomial("cloglog"))
    } else {
        mod <- glm(Y ~ lat, data.frame(lat=X[rn,2]), offset = Off, family=poisson)
    }
    unname(mean(exp(X %*% coef(mod))))
}

naive_est <- function(spp, x, id, truncate = FALSE) {
    d <- nonDuplicated(x, UniqueID, TRUE)
    n <- nrow(d)
    X <- model.matrix(~scale(Latitude), d)
    rn <- rownames(d)

    # distance
    # yydis <- as.matrix(Xtab(Count ~ UniqueID + DISTANCE, x[x$SPECIES == spp,]))[rn,c("0-49 m", "50-100 m", ">100 m")]
    yydis <- as.matrix(Xtab(Count ~ UniqueID + DISTANCE, x[x$SPECIES == spp,]))[rn,c("0-49 m", "50-100 m")]
    # dddis <- matrix(c(0.5, 1, Inf), nrow(yydis), 3, byrow=TRUE)
    Y <- rowSums(yydis)
    Area <- rep(1^2 * pi, n) # 100m radius

    if (missing(id)) {
        rn <- id
    } else {
        rn <- id
    }

    if (truncate) {
        mod <- glm(ifelse(count>0,1,0) ~ lat, data.frame(count = Y[rn], lat=X[rn,2]), offset=log(Area), family=binomial("cloglog"))
    } else {
        mod <- glm(count ~ lat, data.frame(count = Y[rn], lat=X[rn,2]), offset=log(Area), family=poisson)
    }

    unname(mean(exp(X %*% coef(mod))))
}

x <- paired::paired
x1 <- droplevels(x[x$SurveyType == "HUM" & x$Visit == 1 & !is.na(x$Latitude),])
x2 <- droplevels(x[x$SurveyType == "HUM" & x$Visit == 2 & !is.na(x$Latitude),])
# x3 <- droplevels(x[x$SurveyType == "HUM" & x$Visit == 3 & !is.na(x$Latitude),])
ID <- levels(x1$UniqueID)
x <- x1

SPP <- c("CHSP", "DEJU", "HETH", "MAWA", "MYWA", "OVEN", "RCKI", "REVI", "SWTH", "TEWA", "WIWR", "WTSP")
methods = c(
    "SVAaj", "SVAapq",
    "SVAfj", "SVAfpq", 
    "SVAfj05", "SVAfjmc",
    "SVOaj", "SVOapq",
    "SVOfj", "SVOfpq", 
    "SVOfj05", "SVOfjmc",
    "QPAD", "QPADt",
    "NVE", "NVEt")
RES <- NULL

spp <- SPP[spp_id]
for (b in 1:B) {
    if (b == 1) {
        id <- ID
    } else {
        id <- sample(ID, length(ID), replace=TRUE)
    }
    res <- data.frame(iter = b, species = spp, method = methods, estimate = NA_real_, loglik = NA_real_)
    for (met in methods) {
        message("species: ", spp, " - iter: ", b, "/", B, " - method: ", met)
        m <- try(switch(met,
            "SVAfj"  = sqpad_est(spp, x, id, type = "full",   det = "joint"), 
            "SVAfj05"= sqpad_est(spp, x, id, type = "full",   det = "joint", distcorr = 0.5), 
            "SVAfjmc"= sqpad_est(spp, x, id, type = "full",   det = "joint", montecarlo = TRUE), 
            "SVAaj"  = sqpad_est(spp, x, id, type = "approx", det = "joint"), 
            "SVAfpq" = sqpad_est(spp, x, id, type = "full",   det = "pq"), 
            "SVAapq" = sqpad_est(spp, x, id, type = "approx", det = "pq"),
            "SVOfj"  = sqpad_est(spp, x, id, type = "full",   det = "joint", K = 1), 
            "SVOfj05"  = sqpad_est(spp, x, id, type = "full",   det = "joint", K = 1, distcorr = 0.5), 
            "SVOfjmc"  = sqpad_est(spp, x, id, type = "full",   det = "joint", K = 1, montecarlo = TRUE), 
            "SVOaj"  = sqpad_est(spp, x, id, type = "approx", det = "joint", K = 1), 
            "SVOfpq" = sqpad_est(spp, x, id, type = "full",   det = "pq",    K = 1), 
            "SVOapq" = sqpad_est(spp, x, id, type = "approx", det = "pq",    K = 1),
            "QPAD"   = qpad_est(spp, x, id, truncate = FALSE), 
            "QPADt"  = qpad_est(spp, x, id, truncate = TRUE),
            "NVE"    = naive_est(spp, x, id, truncate = FALSE), 
            "NVEt"   = naive_est(spp, x, id, truncate = TRUE)))
        if (!inherits(m, "try-error")) {
            res$estimate[res$method == met] <- m
            if (!is.null(attr(m, "ll")))
                res$loglik[res$method == met] <- attr(m, "ll")
        }
    }
    RES <- rbind(RES, res)
}

f <- file.path(DIR, paste0("estimates_field-data_", spp, ".rds"))
make_dir(f)
saveRDS(RES, f)



EST <- sapply(SPP, function(spp) {
    d <- nonDuplicated(x, UniqueID, TRUE)
    n <- nrow(d)
    X <- model.matrix(~scale(Latitude), d)
    Z <- model.matrix(~TSSR+I(TSSR^2), d)
    rn <- rownames(d)

    # distance
    yydis <- as.matrix(Xtab(Count ~ UniqueID + DISTANCE, x[x$SPECIES == spp,]))[rn,c("0-49 m", "50-100 m")]
    dddis <- matrix(c(0.5, 1), nrow(yydis), 2, byrow=TRUE)
    mdis <- cmulti.fit(Y=yydis, D=dddis, type = "dis")

    # time
    yydur <- as.matrix(Xtab(Count ~ UniqueID + Interval,
        x[x$SPECIES == spp & x$DISTANCE %in% c("0-49 m", "50-100 m"),]))[rn,c("0-3 min", "3-5 min", "5-10 min")]
    dddur <- matrix(c(3, 5, 10), nrow(yydur), 3, byrow=TRUE)
    mdur <- cmulti.fit(Y=yydur, D=dddur, X=NULL, type = "rem")

    # QPAD
    Y <- rowSums(yydis)
    TAU <- exp(coef(mdis))
    PHI <- exp(coef(mdur))
    # Corr <- TAU^2 * pi * (1 - exp(-PHI * 10))
    Area <- rep(1^2 * pi, n) # 100m radius
    p <- 1 - exp(-PHI * 10)
    r <- 1
    q <- TAU^2 * (1 - exp(-r^2/TAU^2))/r^2
    Off <- log(Area * p * q)
    mod <- glm(Y ~ 1, offset = Off, family=poisson)
    c(phi = PHI, tau = TAU, D = unname(exp(coef(mod))))
})
plot(t(EST[-2,]), type="n")
axis(1)
axis(2)
text(t(EST[-2,]), labels=SPP)
# OVEN Ovenbird Seiurus aurocapilla: high phi high D
# DEJU Dark-eyed Junco Junco hyemalis: low phi low D
# WIWR Winter Wren Troglodytes hiemalis: high phi low D
# CHSP Chipping Sparrow Spizella passerina: low phi high D