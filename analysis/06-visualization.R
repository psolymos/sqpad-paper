# make plots for the paper - revision

library(ggplot2)
library(dplyr)
library(arrow)
library(qpdf)

outlier.size = 1
base_size = 12

## --------- bSims simulations -------------

xb0 <- as.data.frame(read_parquet("estimates_bsims_mc.parquet"))

xb200 <- xb0[xb0$nobs == 200,]
xb1000 <- xb0[xb0$nobs == 1000,]

xb <- xb200
methods <- colnames(xb)[grep("_D\\.hat$", colnames(xb))]
xb2 <- NULL
for (i in methods) {
    xb2 <- rbind(xb2,
        data.frame(
            xb[,c("D", "phi", "tau", "M", "E", "C", "J", "nobs", "maxdur", "maxdis", "run")],
            Estimate = xb[,i],
            Method = i,
            LogLik = if (gsub("_D\\.hat$", "_LogLik", i) %in% colnames(xb)) xb[,gsub("_D\\.hat$", "_LogLik", i)] else NA_real_
        )
    )
}
xb2$Estimate[xb2$Estimate > 10*xb2$D] <- NA
xb2$RelBias <- (xb2$Estimate - xb2$D) / xb2$D
xb2$Method <- factor(xb2$Method, methods)
levels(xb2$Method) <- gsub("_D\\.hat$", "", levels(xb2$Method))
xb2$setup <- paste0("M", xb2$M, "-E", xb2$E, "-C", xb2$C, "-J", xb2$J)
xb2$Density <- paste("Density =", xb2$D)
xb2$CueRate <- paste("Cue Rate =", xb2$phi)

table(xb2$Method, is.na(xb2$RelBias))

xbD <- droplevels(xb2[xb2$Method %in% c("CONV", "SVAfj",  "SVOfj", "QPAD", "PAD", "QAD", "NVE"),])
xbD$Method <- factor(as.character(xbD$Method), c("CONV", "SVOfj",   "SVAfj",   "QPAD", "PAD", "QAD", "NVE"))
levels(xbD$Method) <-                          c("Conv", "SQPAD-O", "SQPAD-C", "QPAD", "PAD", "QAD", "Naïve")

mnD <- xbD |> group_by(setup, Method, CueRate, Density) |>
    summarize(RelBias = mean(RelBias, na.rm=T)) |> as.data.frame()

plD1 <- xbD[xbD$setup == "M0-E0-C0-J1",] |> 
    ggplot(aes(x = Method, y = RelBias)) + 
    geom_boxplot(outlier.shape = 1, outlier.size = outlier.size) +
    # geom_boxplot() +
    geom_abline(slope = 0, intercept = 0, lty = 2) +
    geom_point(aes(x = Method, y = RelBias), mnD[mnD$setup == "M0-E0-C0-J1",], shape=4, size=2.5) +
    # coord_flip() +
    # ylim(-1, NA) +
    coord_cartesian(ylim = c(-1, 5)) +
    theme_bw(base_size = base_size) +
    facet_grid(rows=vars(CueRate), cols=vars(Density)) +
    labs(x = "Estimator", y = "Relative Bias (Density)")

plD2 <- xbD[xbD$setup == "M1-E1-C1-J1",] |> 
    ggplot(aes(x = Method, y = RelBias)) + 
    geom_boxplot(outlier.shape = 1, outlier.size = outlier.size) +
    # geom_boxplot(outlier.shape = NA) +
    geom_abline(slope = 0, intercept = 0, lty = 2) +
    geom_point(aes(x = Method, y = RelBias), mnD[mnD$setup == "M1-E1-C1-J1",], shape=4, size=2.5) +
    # ylim(-1, NA) +
    coord_cartesian(ylim = c(-1, 5)) +
    theme_bw(base_size = base_size) +
    # coord_flip() +
    facet_grid(rows=vars(CueRate), cols=vars(Density)) +
    labs(x = "Estimator", y = "Relative Bias (Density)")

plD3 <- xbD[xbD$setup == "M1-E1-C0-J1",] |> 
    ggplot(aes(x = Method, y = RelBias)) + 
    geom_boxplot(outlier.shape = 1, outlier.size = outlier.size) +
    # geom_boxplot(outlier.shape = NA) +
    geom_abline(slope = 0, intercept = 0, lty = 2) +
    geom_point(aes(x = Method, y = RelBias), mnD[mnD$setup == "M1-E1-C1-J1",], shape=4, size=2.5) +
    # ylim(-1, NA) +
    coord_cartesian(ylim = c(-1, 5)) +
    theme_bw(base_size = base_size) +
    # coord_flip() +
    facet_grid(rows=vars(CueRate), cols=vars(Density)) +
    labs(x = "Estimator", y = "Relative Bias (Density)")

plD4 <- xbD[xbD$setup == "M1-E0-C1-J1",] |> 
    ggplot(aes(x = Method, y = RelBias)) + 
    geom_boxplot(outlier.shape = 1, outlier.size = outlier.size) +
    # geom_boxplot(outlier.shape = NA) +
    geom_abline(slope = 0, intercept = 0, lty = 2) +
    geom_point(aes(x = Method, y = RelBias), mnD[mnD$setup == "M1-E1-C1-J1",], shape=4, size=2.5) +
    # ylim(-1, NA) +
    coord_cartesian(ylim = c(-1, 5)) +
    theme_bw(base_size = base_size) +
    # coord_flip() +
    facet_grid(rows=vars(CueRate), cols=vars(Density)) +
    labs(x = "Estimator", y = "Relative Bias (Density)")

# Cue rate

methods2 <- colnames(xb)[grep("_phi\\.hat$", colnames(xb))]
xb2 <- NULL
for (i in methods2) {
    xb2 <- rbind(xb2,
        data.frame(
            xb[,c("D", "phi", "tau", "M", "E", "C", "J", "nobs", "maxdur", "maxdis", "run")],
            Estimate = xb[,i],
            Method = i
        )
    )
}
xb2$Estimate[xb2$Estimate > 10*xb2$phi] <- NA
xb2$RelBias <- (xb2$Estimate - xb2$phi) / xb2$phi
xb2$Method <- factor(xb2$Method, methods2)
levels(xb2$Method) <- gsub("_phi\\.hat$", "", levels(xb2$Method))
xb2$setup <- paste0("M", xb2$M, "-E", xb2$E, "-C", xb2$C, "-J", xb2$J)
xb2$Density <- paste("Density =", xb2$D)
xb2$CueRate <- paste("Cue Rate =", xb2$phi)

table(xb2$Method, is.na(xb2$RelBias))

xbPhi <- droplevels(xb2[xb2$Method %in% c("CONV", "SVAfj",  "SVOfj", "QPAD", "PAD", "NVE"),])
xbPhi$Method <- factor(as.character(xbPhi$Method), c("CONV", "SVOfj",   "SVAfj",   "QPAD", "PAD", "NVE"))
levels(xbPhi$Method) <-                            c("Conv", "SQPAD-O", "SQPAD-C", "QPAD", "PAD", "Naïve")

mnPhi <- xbPhi |> group_by(setup, Method, CueRate, Density) |>
    summarize(RelBias = mean(RelBias, na.rm=T)) |> as.data.frame()

plPhi1 <- xbPhi[xbPhi$setup == "M0-E0-C0-J1",] |> 
    ggplot(aes(x = Method, y = RelBias)) + 
    geom_boxplot(outlier.shape = 1, outlier.size = outlier.size) +
    # geom_boxplot() +
    geom_abline(slope = 0, intercept = 0, lty = 2) +
    geom_point(aes(x = Method, y = RelBias), mnPhi[mnPhi$setup == "M0-E0-C0-J1",], shape=4, size=2.5) +
    # coord_flip() +
    # ylim(-1, NA) +
    coord_cartesian(ylim = c(-1, 2)) +
    theme_bw(base_size = base_size) +
    facet_grid(rows=vars(CueRate), cols=vars(Density)) +
    labs(x = "Estimator", y = "Relative Bias (Cue Rate)")

plPhi2 <- xbPhi[xbPhi$setup == "M1-E1-C1-J1",] |> 
    ggplot(aes(x = Method, y = RelBias)) + 
    geom_boxplot(outlier.shape = 1, outlier.size = outlier.size) +
    # geom_boxplot(outlier.shape = NA) +
    geom_abline(slope = 0, intercept = 0, lty = 2) +
    geom_point(aes(x = Method, y = RelBias), mnPhi[mnPhi$setup == "M1-E1-C1-J1",], shape=4, size=2.5) +
    # ylim(-1, NA) +
    coord_cartesian(ylim = c(-1, 2)) +
    theme_bw(base_size = base_size) +
    # coord_flip() +
    facet_grid(rows=vars(CueRate), cols=vars(Density)) +
    labs(x = "Estimator", y = "Relative Bias (Cue Rate)")

plPhi3 <- xbPhi[xbPhi$setup == "M1-E1-C0-J1",] |> 
    ggplot(aes(x = Method, y = RelBias)) + 
    geom_boxplot(outlier.shape = 1, outlier.size = outlier.size) +
    # geom_boxplot(outlier.shape = NA) +
    geom_abline(slope = 0, intercept = 0, lty = 2) +
    geom_point(aes(x = Method, y = RelBias), mnPhi[mnPhi$setup == "M1-E1-C0-J1",], shape=4, size=2.5) +
    # ylim(-1, NA) +
    coord_cartesian(ylim = c(-1, 2)) +
    theme_bw(base_size = base_size) +
    # coord_flip() +
    facet_grid(rows=vars(CueRate), cols=vars(Density)) +
    labs(x = "Estimator", y = "Relative Bias (Cue Rate)")

plPhi4 <- xbPhi[xbPhi$setup == "M1-E0-C1-J1",] |> 
    ggplot(aes(x = Method, y = RelBias)) + 
    geom_boxplot(outlier.shape = 1, outlier.size = outlier.size) +
    # geom_boxplot(outlier.shape = NA) +
    geom_abline(slope = 0, intercept = 0, lty = 2) +
    geom_point(aes(x = Method, y = RelBias), mnPhi[mnPhi$setup == "M1-E0-C1-J1",], shape=4, size=2.5) +
    # ylim(-1, NA) +
    coord_cartesian(ylim = c(-1, 2)) +
    theme_bw(base_size = base_size) +
    # coord_flip() +
    facet_grid(rows=vars(CueRate), cols=vars(Density)) +
    labs(x = "Estimator", y = "Relative Bias (Cue Rate)")

# Distance parameter

methods2 <- colnames(xb)[grep("_tau\\.hat$", colnames(xb))]
xb2 <- NULL
for (i in methods2) {
    xb2 <- rbind(xb2,
        data.frame(
            xb[,c("D", "phi", "tau", "M", "E", "C", "J", "nobs", "maxdur", "maxdis", "run")],
            Estimate = xb[,i],
            Method = i
        )
    )
}
xb2$Estimate[xb2$Estimate > 10*xb2$tau] <- NA
xb2$RelBias <- (xb2$Estimate - xb2$tau) / xb2$tau
xb2$Method <- factor(xb2$Method, methods2)
levels(xb2$Method) <- gsub("_tau\\.hat$", "", levels(xb2$Method))
xb2$setup <- paste0("M", xb2$M, "-E", xb2$E, "-C", xb2$C, "-J", xb2$J)
xb2$Density <- paste("Density =", xb2$D)
xb2$CueRate <- paste("Cue Rate =", xb2$phi)

table(xb2$Method, is.na(xb2$RelBias))

xbTau <- droplevels(xb2[xb2$Method %in% c("CONV", "SVAfj",  "SVOfj", "QPAD", "QAD", "NVE"),])
xbTau$Method <- factor(as.character(xbTau$Method), c("CONV", "SVOfj",   "SVAfj",   "QPAD", "QAD", "NVE"))
levels(xbTau$Method) <-                            c("Conv", "SQPAD-O", "SQPAD-C", "QPAD", "QAD", "Naïve")

mnTau <- xbTau |> group_by(setup, Method, CueRate, Density) |>
    summarize(RelBias = mean(RelBias, na.rm=T)) |> as.data.frame()

plTau1 <- xbTau[xbTau$setup == "M0-E0-C0-J1",] |> 
    ggplot(aes(x = Method, y = RelBias)) + 
    geom_boxplot(outlier.shape = 1, outlier.size = outlier.size) +
    # geom_boxplot() +
    geom_abline(slope = 0, intercept = 0, lty = 2) +
    geom_point(aes(x = Method, y = RelBias), mnTau[mnTau$setup == "M0-E0-C0-J1",], shape=4, size=2.5) +
    # coord_flip() +
    # ylim(-1, NA) +
    coord_cartesian(ylim = c(-1, 1)) +
    theme_bw(base_size = base_size) +
    facet_grid(rows=vars(CueRate), cols=vars(Density)) +
    labs(x = "Estimator", y = "Relative Bias (Distance Parameter)")

plTau2 <- xbTau[xbTau$setup == "M1-E1-C1-J1",] |> 
    ggplot(aes(x = Method, y = RelBias)) + 
    geom_boxplot(outlier.shape = 1, outlier.size = outlier.size) +
    # geom_boxplot(outlier.shape = NA) +
    geom_abline(slope = 0, intercept = 0, lty = 2) +
    geom_point(aes(x = Method, y = RelBias), mnTau[mnTau$setup == "M1-E1-C1-J1",], shape=4, size=2.5) +
    # ylim(-1, NA) +
    coord_cartesian(ylim = c(-1, 1)) +
    theme_bw(base_size = base_size) +
    # coord_flip() +
    facet_grid(rows=vars(CueRate), cols=vars(Density)) +
    labs(x = "Estimator", y = "Relative Bias", title = "(B) Distance Parameter")

plTau3 <- xbTau[xbTau$setup == "M1-E1-C0-J1",] |> 
    ggplot(aes(x = Method, y = RelBias)) + 
    geom_boxplot(outlier.shape = 1, outlier.size = outlier.size) +
    # geom_boxplot(outlier.shape = NA) +
    geom_abline(slope = 0, intercept = 0, lty = 2) +
    geom_point(aes(x = Method, y = RelBias), mnTau[mnTau$setup == "M1-E1-C0-J1",], shape=4, size=2.5) +
    # ylim(-1, NA) +
    coord_cartesian(ylim = c(-1, 1)) +
    theme_bw(base_size = base_size) +
    # coord_flip() +
    facet_grid(rows=vars(CueRate), cols=vars(Density)) +
    labs(x = "Estimator", y = "Relative Bias", title = "(B) Distance Parameter")

plTau4 <- xbTau[xbTau$setup == "M1-E0-C1-J1",] |> 
    ggplot(aes(x = Method, y = RelBias)) + 
    geom_boxplot(outlier.shape = 1, outlier.size = outlier.size) +
    # geom_boxplot(outlier.shape = NA) +
    geom_abline(slope = 0, intercept = 0, lty = 2) +
    geom_point(aes(x = Method, y = RelBias), mnTau[mnTau$setup == "M1-E0-C1-J1",], shape=4, size=2.5) +
    # ylim(-1, NA) +
    coord_cartesian(ylim = c(-1, 1)) +
    theme_bw(base_size = base_size) +
    # coord_flip() +
    facet_grid(rows=vars(CueRate), cols=vars(Density)) +
    labs(x = "Estimator", y = "Relative Bias", title = "(B) Distance Parameter")

ggsave("fig2-bsims-0001.pdf", plD1, width=10, height=8)
# ggsave("fig3-bsims-1111.pdf", plD2, width=10, height=8)
ggsave("fig3-bsims-1101.pdf", plD3, width=10, height=8)
ggsave("fig4-bsims-1011.pdf", plD4, width=10, height=8)

ggsave("fig2-bsims-0001.png", plD1, width=10, height=8)
ggsave("fig3-bsims-1101.png", plD3, width=10, height=8)
ggsave("fig4-bsims-1011.png", plD4, width=10, height=8)

ggsave("fig2-bsims-0001phi.pdf", plPhi1, width=10, height=8)
ggsave("fig3-bsims-1111phi.pdf", plPhi2, width=10, height=8)

ggsave("fig2-bsims-0001tau.pdf", plTau1, width=10, height=8)
ggsave("fig3-bsims-1111tau.pdf", plTau2, width=10, height=8)

qpdf::pdf_combine(
    input = c(
        "fig2-bsims-0001.pdf", 
        "fig2-bsims-0001phi.pdf", 
        "fig2-bsims-0001tau.pdf", 
        "fig3-bsims-1111.pdf", 
        "fig3-bsims-1111phi.pdf", 
        "fig3-bsims-1111tau.pdf"),
    output = "figs-bsims-D-phi-tau.pdf")

pdf("figs-bsims-all-setups-D-phi-tau.pdf", width=10, height=8, onefile = TRUE)
for (i in 1:4) {
    setup <- c("M0-E0-C0-J1", "M1-E0-C1-J1", "M1-E1-C0-J1", "M1-E1-C1-J1")[i]
    setup_text <- c("No bias", "Double counting", "Distance estimation error", "Both biases")[i]

    plD <- xbD[xbD$setup == setup,] |> 
        ggplot(aes(x = Method, y = RelBias)) + 
        geom_boxplot(outlier.shape = 1, outlier.size = outlier.size) +
        # geom_boxplot() +
        geom_abline(slope = 0, intercept = 0, lty = 2) +
        geom_point(aes(x = Method, y = RelBias), mnD[mnD$setup == setup,]) +
        # coord_flip() +
        # ylim(-1, NA) +
        coord_cartesian(ylim = c(-1, 5)) +
        theme_bw(base_size = base_size) +
        facet_grid(rows=vars(CueRate), cols=vars(Density)) +
        labs(x = "Estimator", y = "Relative Bias (Density)", title = setup_text)
    plPhi <- xbPhi[xbPhi$setup == setup,] |> 
        ggplot(aes(x = Method, y = RelBias)) + 
        geom_boxplot(outlier.shape = 1, outlier.size = outlier.size) +
        # geom_boxplot() +
        geom_abline(slope = 0, intercept = 0, lty = 2) +
        geom_point(aes(x = Method, y = RelBias), mnPhi[mnPhi$setup == setup,]) +
        # coord_flip() +
        # ylim(-1, NA) +
        coord_cartesian(ylim = c(-1, 2)) +
        theme_bw(base_size = base_size) +
        facet_grid(rows=vars(CueRate), cols=vars(Density)) +
        labs(x = "Estimator", y = "Relative Bias (Cue Rate)", title = setup_text)
    plTau <- xbTau[xbTau$setup == setup,] |> 
        ggplot(aes(x = Method, y = RelBias)) + 
        geom_boxplot(outlier.shape = 1, outlier.size = outlier.size) +
        # geom_boxplot() +
        geom_abline(slope = 0, intercept = 0, lty = 2) +
        geom_point(aes(x = Method, y = RelBias), mnTau[mnTau$setup == setup,]) +
        # coord_flip() +
        # ylim(-1, NA) +
        coord_cartesian(ylim = c(-1, 1)) +
        theme_bw(base_size = base_size) +
        facet_grid(rows=vars(CueRate), cols=vars(Density)) +
        labs(x = "Estimator", y = "Relative Bias (Distance Parameter)", title = setup_text)
    print(plD)
    print(plPhi)
    print(plTau)

}
dev.off()

## --------- field data -------------

xf <- as.data.frame(read_parquet("estimates_field_data_mc.parquet"))

# 4 species
# OVEN Ovenbird Seiurus aurocapilla: high phi high D
# DEJU Dark-eyed Junco Junco hyemalis: low phi low D

# WIWR Winter Wren Troglodytes hiemalis: high phi low D
# SWTH Swainson's Thrush Catharus ustulatus: high phi low D

# CHSP Chipping Sparrow Spizella passerina: low phi high D

# make 1 plot (2x2)

z <- xf[xf$species %in% c("OVEN", "DEJU", "SWTH", "CHSP") &
    xf$method %in% c("NVE", "QPAD", "SVAfj", "SVOfj"),]

zq <- z |> group_by(method, species) |>
    summarize(
        mean = mean(estimate, na.rm=TRUE),
        median = median(estimate, na.rm=TRUE),
        lower = quantile(estimate, 0.025, na.rm=TRUE),
        upper = quantile(estimate, 0.975, na.rm=TRUE)
    ) |> as.data.frame()
zq$method <- factor(zq$method, rev(c("NVE",   "QPAD", "SVAfj",   "SVOfj")))
levels(zq$method) <- rev(c(          "Naïve", "QPAD", "SQPAD-C", "SQPAD-O"))
zq$shape <- 19
zq$shape[zq$method %in% c("Naïve", "QPAD")] <- 21

plField <- zq |> ggplot(aes(x = method, y = median, ymin = lower, ymax = upper)) +
    geom_point(size = 2, shape = zq$shape) +
    geom_linerange() +
    coord_flip() +
    theme_bw(base_size = base_size) +
    ylim(0, NA) +
    facet_wrap(vars(species), scales = "free_x") +
    labs(x = "Estimator", y = "Density [males/ha]")

ggsave("fig6-species.pdf", plField)
ggsave("fig6-species.png", plField)


# convolution plot
png("fig7-convolution.png", width=2.5*1200, height=2.5*500, res=300)

op <- par(mfrow=c(1,2), mar=c(2,2,2,2))

st <- 20

d <- seq(0, 200, st)
plot(0, 0, asp=1, type="p", pch="+", ann=F, axes=F, xlim=c(-1,1)*max(d), ylim=c(-1,1)*max(d))
for (i in d)
    bSims:::.draw_ellipse(0, 0, a = i, b = i, lty = 1, lwd = if (i==max(d)) 1 else 0.5)

par(op[2])

tau <- 1
d2 <- seq(0, 200, 0.1)
p2 <- exp(-(d2/100)^2/tau^2)
pt <- seq(0, max(d2), st)
xx <- pt[-1] - st/2
yy <- exp(-(xx/100)^2/tau^2)
# md <- aggregate(p2, list(cut(d2, pt, include.lowest=T)), mean)
plot(d2, p2, type="l", xlab="Distance [m]", ylab="Probability", axes=F)
axis(1)
axis(2)
plot(stepfun(pt[-length(pt)], c(1,yy)), add=T, do.points=F)
lines(pt, c(1, yy), type="h", lwd=0.5)
abline(h=0,lwd=0.5)

par(op)
dev.off()
