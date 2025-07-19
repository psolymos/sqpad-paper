# SQPAD helper functions for simulations

p_fun <- function (t, phi) {
    1 - exp(-t * phi)
}

q_fun <- function (r, tau) {
    tau^2 * (1 - exp(-r^2/tau^2))/r^2
}

make_dir <- function(file) {
    dir <- dirname(file)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    invisible(NULL)
}

est_fun <- function(n, all_sims,
    replace = FALSE,
    verbose = FALSE) {

    dur_vals <-  c(2, 3, 4, 5, 6, 7, 8, 9, 10)
    dur_probs <- c(1, 1, 1, 1, 1, 1, 1, 1, 1)
    nbinst <- 3

    dis_vals <- c(40, 60, 80, 100, 120, 140, 160, 180, 200)/100
    dis_probs <- c(1,  1,  1,  1,   1,   1,   1,   1,   1)
    nbinsr <- 3

    if (verbose) message("STARTING")
    # take a random sample w/o replacement
    ind <- sample(length(all_sims), n, replace=replace)
    olist <- all_sims[ind]

    # empty data to store results in
    SV <- REM <- DIS <- DREM <- DDIS <- NULL

    # data: take one site at a time
    for (j in seq_len(n)) {

        o <- olist[[j]]

        ti <- sample(dur_vals, 1, replace = FALSE, prob = dur_probs)
        ri <- sample(dis_vals, 1, replace = FALSE, prob = dis_probs)
        tint <- seq(0, ti, ti/nbinst)[-1]
        rint <- seq(0, ri, ri/nbinsr)[-1]

        d1 <- o$detections$removal
        d1 <- d1[!duplicated(d1$j),] # perceived individual ID is unique
        d1$derr <- d1$d + d1$error
        d1$rint_qpad <- cut(d1$derr, c(0, rint), include.lowest=TRUE)
        d1$tint_qpad <- cut(d1$t, c(0, tint), include.lowest=TRUE)

        d2 <- d1[d1$t < ti & d1$derr <= ri,] # duration and (perceived) distance respected
        xt1 <- as.matrix(Xtab(~ rint_qpad + tint_qpad, d1[!is.na(d1$rint_qpad) & !is.na(d1$tint_qpad),,drop=FALSE]))

        d_ve <- data.frame(
            Y = nrow(d2),
            dur = ti,
            dis = ri,
            dij = NA)
        d_ve$dij[1] <- list(d2$derr)
        d_rem <- colSums(xt1)
        d_dis <- rowSums(xt1)

        SV <- rbind(SV, d_ve)
        REM <- rbind(REM, d_rem)
        DREM <- rbind(DREM, tint)
        DIS <- rbind(DIS, d_dis)
        DDIS <- rbind(DDIS, rint)
    }

    # estimation
    x <- SV
    x$Area <- x$dis^2 * pi # truncation for MV etc.

    Try <- function(..., prefix = "", verbose = FALSE) {
        if (verbose)
            message("--- ", prefix)
        x <- try(...)
        if (inherits(x, "try-error")) {
            o <- rep(NA_real_, 4)
        } else {
            o <- c(exp(x$coef), x$loglik)
        }
        names(o) <- paste(prefix, c("D.hat", "phi.hat", "tau.hat", "LogLik"), sep = "_")
        o
    }

    res <- numeric(0)

    Y <- x$Y
    dur <- x$dur
    dis <- x$dis

    r0 <- Try(sqpad.fit(Y = Y, dis = dis, dur = dur, dislist = x$dij, type = "conv", det = "joint", K = NULL, hessian = FALSE, Nmax = 20),
        prefix = "CONV", verbose = verbose)
    res <- c(res, r0)

    r1 <- Try(sqpad.fit(Y = Y, dis = dis, dur = dur, type = "approx", det = "joint", K = NULL, hessian = FALSE),
        prefix = "SVAaj", verbose = verbose)
    res <- c(res, r1)
    r2 <- Try(sqpad.fit(Y = Y, dis = dis, dur = dur, type = "full", det = "joint", K = NULL, hessian = FALSE),
        prefix = "SVAfj", verbose = verbose)
    res <- c(res, r2)

    r5 <- Try(sqpad.fit(Y = Y, dis = dis, dur = dur, type = "approx", det = "joint", K = 1, hessian = FALSE),
        prefix = "SVOaj", verbose = verbose)
    res <- c(res, r5)
    r6 <- Try(sqpad.fit(Y = Y, dis = dis, dur = dur, type = "full", det = "joint", K = 1, hessian = FALSE),
        prefix = "SVOfj", verbose = verbose)
    res <- c(res, r6)

    # tru different approaches to p
    r2_05 <- Try(sqpad.fit(Y = Y, dis = dis, dur = dur, type = "full", det = "joint", K = NULL, hessian = FALSE, distcorr = 0.5),
        prefix = "SVAfj05", verbose = verbose)
    res <- c(res, r2_05)
    r2_mc <- Try(sqpad.fit(Y = Y, dis = dis, dur = dur, type = "full", det = "joint", K = NULL, hessian = FALSE, montecarlo = TRUE),
        prefix = "SVAfjmc", verbose = verbose)
    res <- c(res, r2_mc)

    r6_05 <- Try(sqpad.fit(Y = Y, dis = dis, dur = dur, type = "full", det = "joint", K = 1, hessian = FALSE, distcorr = 0.5),
        prefix = "SVOfj05", verbose = verbose)
    res <- c(res, r6_05)
    r6_mc <- Try(sqpad.fit(Y = Y, dis = dis, dur = dur, type = "full", det = "joint", K = 1, hessian = FALSE, montecarlo = TRUE),
        prefix = "SVOfjmc", verbose = verbose)
    res <- c(res, r6_mc)


    if (verbose) message("--- QPAD")
    fit_rem <- try(cmulti.fit(
        Y = REM,
        D = DREM,
        X = NULL,
        type = "rem"))
    fit_dis <- try(cmulti.fit(
        Y = DIS,
        D = DDIS,
        X = NULL,
        type = "dis"))

    if (inherits(fit_rem, "try-error") || inherits(fit_dis, "try-error")) {
        est_QPAD <- rep(NA_real_, 3)
        est_QPADt <- rep(NA_real_, 1)
        est_PAD <- rep(NA_real_, 2)
        est_QAD <- rep(NA_real_, 2)
    } else {
        p <- p_fun(ti, exp(coef(fit_rem)))
        q <- q_fun(ri, exp(coef(fit_dis)))
        A <- ri^2*pi

        off1 <- rep(log(A*p*q), nrow(x))
        eqp <- glm(Y ~ 1 + offset(off1), data=x, family=poisson)
        est_QPAD <- exp(c(coef(eqp), coef(fit_rem), coef(fit_dis)))

        eqp_t <- glm(ifelse(Y > 0, 1, 0) ~ 1 + offset(off1), data=x, family=binomial(link = "cloglog"))
        est_QPADt <- exp(coef(eqp_t))

        off2 <- rep(log(A*p), nrow(x))
        ep <- glm(Y ~ 1 + offset(off2), data=x, family=poisson)
        est_PAD <- exp(c(coef(ep), coef(fit_rem)))

        off3 <- rep(log(A*q), nrow(x))
        eq <- glm(Y ~ 1 + offset(off3), data=x, family=poisson)
        est_QAD <- exp(c(coef(eq), coef(fit_dis)))
    }
    names(est_QPAD) <- paste0("QPAD_", c("D.hat", "phi.hat", "tau.hat"))
    names(est_QPADt) <- paste0("QPADt_", c("D.hat"))
    names(est_PAD) <- paste0("PAD_", c("D.hat", "phi.hat"))
    names(est_QAD) <- paste0("QAD_", c("D.hat", "tau.hat"))
    res <- c(res, est_QPAD, est_QPADt, est_PAD, est_QAD)

    if (verbose) message("--- Naive")
    naive <- glm(Y ~ 1 + offset(log(Area)), data=x, family=poisson)
    est_NAIVE <- exp(coef(naive))
    names(est_NAIVE) <- paste0("NVE_", c("D.hat"))

    naive_t <- glm(ifelse(Y > 0, 1, 0) ~ 1 + offset(log(Area)), data=x, family=binomial(link = "cloglog"))
    est_NAIVEt <- exp(coef(naive_t))
    names(est_NAIVEt) <- paste0("NVEt_", c("D.hat"))
    res <- c(res, est_NAIVE, est_NAIVEt)

    if (verbose) message("--- Assembly")

    e <- c(
        nobs=n,
        maxdur=max(dur_vals), maxdis=max(dis_vals),
        res)
    if (verbose) message("DONE")
    e
}

est_fun_QPAD <- function(n, all_sims,
    dur_max = 10,
    dis_max = 2,
    dur_nbins = 10,
    dis_nbins = 10,
    verbose = FALSE
) {

    tint <- seq(0, dur_max, dur_max/dur_nbins)[-1]
    rint <- seq(0, dis_max, dis_max/dis_nbins)[-1]
    ti <- max(tint)
    ri <- max(rint)

    if (verbose) message("STARTING")
    # take a random sample w/o replacement
    ind <- sample(length(all_sims), n)
    olist <- all_sims[ind]

    # empty data to store results in
    SV <- REM <- DIS <- DREM <- DDIS <- NULL

    # data: take one site at a time
    for (j in seq_len(n)) {

        o <- olist[[j]]

        d1 <- o$detections$removal
        d1 <- d1[!duplicated(d1$j),] # perceived individual ID is unique
        d1$derr <- d1$d + d1$error
        d1$rint_qpad <- cut(d1$derr, c(0, rint), include.lowest=TRUE)
        d1$tint_qpad <- cut(d1$t, c(0, tint), include.lowest=TRUE)

        d2 <- d1[d1$t < ti & d1$derr <= ri,] # duration and (perceived) distance respected
        xt1 <- as.matrix(Xtab(~ rint_qpad + tint_qpad, d1[!is.na(d1$rint_qpad) & !is.na(d1$tint_qpad),,drop=FALSE]))

        d_ve <- data.frame(
            Y = nrow(d2),
            dur = ti,
            dis = ri,
            dij = NA)
        d_ve$dij[1] <- list(d2$derr)
        d_rem <- colSums(xt1)
        d_dis <- rowSums(xt1)

        SV <- rbind(SV, d_ve)
        REM <- rbind(REM, d_rem)
        DREM <- rbind(DREM, tint)
        DIS <- rbind(DIS, d_dis)
        DDIS <- rbind(DDIS, rint)
    }

    # estimation
    x <- SV
    x$Area <- x$dis^2 * pi # truncation for MV etc.

    Try <- function(..., prefix = "", verbose = FALSE) {
        if (verbose)
            message("--- ", prefix)
        x <- try(...)
        if (inherits(x, "try-error")) {
            o <- rep(NA_real_, 4)
        } else {
            o <- c(exp(x$coef), x$loglik)
        }
        names(o) <- paste(prefix, c("D.hat", "phi.hat", "tau.hat", "LogLik"), sep = "_")
        o
    }

    res <- numeric(0)

    if (verbose) message("--- QPAD")
    fit_rem <- try(cmulti.fit(
        Y = REM,
        D = DREM,
        X = NULL,
        type = "rem"))
    fit_dis <- try(cmulti.fit(
        Y = DIS,
        D = DDIS,
        X = NULL,
        type = "dis"))

    if (inherits(fit_rem, "try-error") || inherits(fit_dis, "try-error")) {
        est_QPAD <- rep(NA_real_, 3)
        est_QPADt <- rep(NA_real_, 1)
        est_PAD <- rep(NA_real_, 2)
        est_QAD <- rep(NA_real_, 2)
    } else {
        p <- p_fun(ti, exp(coef(fit_rem)))
        q <- q_fun(ri, exp(coef(fit_dis)))
        A <- ri^2*pi

        off1 <- rep(log(A*p*q), nrow(x))
        eqp <- glm(Y ~ 1 + offset(off1), data=x, family=poisson)
        est_QPAD <- exp(c(coef(eqp), coef(fit_rem), coef(fit_dis)))

        eqp_t <- glm(ifelse(Y > 0, 1, 0) ~ 1 + offset(off1), data=x, family=binomial(link = "cloglog"))
        est_QPADt <- exp(coef(eqp_t))

        off2 <- rep(log(A*p), nrow(x))
        ep <- glm(Y ~ 1 + offset(off2), data=x, family=poisson)
        est_PAD <- exp(c(coef(ep), coef(fit_rem)))

        off3 <- rep(log(A*q), nrow(x))
        eq <- glm(Y ~ 1 + offset(off3), data=x, family=poisson)
        est_QAD <- exp(c(coef(eq), coef(fit_dis)))
    }
    names(est_QPAD) <- paste0("QPAD_", c("D.hat", "phi.hat", "tau.hat"))
    names(est_QPADt) <- paste0("QPADt_", c("D.hat"))
    names(est_PAD) <- paste0("PAD_", c("D.hat", "phi.hat"))
    names(est_QAD) <- paste0("QAD_", c("D.hat", "tau.hat"))
    res <- c(res, est_QPAD, est_QPADt, est_PAD, est_QAD)

    if (verbose) message("--- Naive")
    naive <- glm(Y ~ 1 + offset(log(Area)), data=x, family=poisson)
    est_NAIVE <- exp(coef(naive))
    names(est_NAIVE) <- paste0("NVE_", c("D.hat"))

    if (verbose) message("--- Assembly")

    e <- c(
        nobs=n,
        maxdur=ti, maxdis=ri,
        res)
    if (verbose) message("DONE")
    e
}
