## ----"knitr-options", echo = FALSE--------------------------------------------
library(knitr)
opts_chunk$set(size = "small",
               fig.height = 4,
               fig.align = "center",
               cache = FALSE,
               message = FALSE,
               warning = FALSE)


## ----"data-loading"-----------------------------------------------------------
library(ReplicationSuccess)
data("RProjects")
str(RProjects, width = 80, strict.width = "cut")

## computing zo, zr, c
RProjects$zo <- with(RProjects, fiso/se_fiso)
RProjects$zr <- with(RProjects, fisr/se_fisr)
RProjects$c <- with(RProjects, se_fiso^2/se_fisr^2)

## ----"plot-projects", fig.height = 6------------------------------------------
## plots of effect estimates

par(mfrow = c(2, 2), las = 1, pty = "s",
    mar = c(4, 2.5, 2.5, 1))
for (p in unique(RProjects$project)) {
    data_project <- subset(RProjects, project == p)
    significant <- ifelse(data_project$pr1 < 0.025, "#DF536BCC", "#00000080")
    plot(rr ~ ro, data = data_project, ylim = c(-0.5, 1), col = significant,
         xlim = c(-0.5, 1), main = p, xlab = expression(italic(r)[o]),
         cex = 0.7, pch = 19, ylab = expression(italic(r)[r]))
    legend("topleft", legend = "replication significant", cex = 0.8, pch = 20,
           col = "#DF536BCC", bty = "n")
    abline(h = 0, lty = 2)
    abline(a = 0, b = 1, col = "grey")
}

## ----forestPlots, fig.height = 5----------------------------------------------
data("protzko2020")

## forestplots of effect estimates
parOld <- par(mar = c(5, 8, 4, 2), mfrow = c(2, 2))
experiments <- unique(protzko2020$experiment)[1:4]
for (ex in experiments) {
  ## compute CIs
  dat <- subset(protzko2020, experiment == ex)
  za <- qnorm(p = 0.975)
  plotDF <- data.frame(lower = dat$smd - za*dat$se,
                       est = dat$smd,
                       upper = dat$smd + za*dat$se)
colpalette <- c("#000000", "#1B9E77", "#D95F02")
cols <- colpalette[dat$type]
yseq <- seq(1, nrow(dat))

## forestplot
plot(x = plotDF$est, y = yseq, xlim = c(-0.15, 0.8),
     ylim = c(0.8*min(yseq), 1.05*max(yseq)), type = "n",
     yaxt = "n", xlab = "Effect estimate (SMD)", ylab = "")
abline(v = 0, col = "#0000004D")
arrows(x0 = plotDF$lower, x1 = plotDF$upper, y0 = yseq, angle = 90,
       code = 3, length = 0.05, col = cols)
points(y = yseq, x = plotDF$est, pch = 20, lwd = 2, col = cols)
axis(side = 2, at = yseq, las = 1, labels = dat$type, cex.axis = 0.85)
title(main = ex)
}

## -----------------------------------------------------------------------------
for (p in unique(RProjects$project)) {
    data_project <- subset(RProjects, project == p)
    significant_O <- data_project$po1 < 0.025
    significant_R <- data_project$pr1 < 0.025
    success <- significant_O & significant_R
    cat(paste0(p, ": \n"))
    cat(paste0(round(mean(significant_O)*100, 1), "% original studies significant (",
               sum(significant_O), "/", length(significant_O), ")\n"))
    cat(paste0(round(mean(significant_R)*100, 1), "% replications significant (",
               sum(significant_R), "/", length(significant_R), ")\n"))
    cat(paste0(round(mean(success)*100, 1),
               "% both studies significant in the same direction (",
               sum(success), "/", length(success), ")\n \n"))
}

## -----------------------------------------------------------------------------
sampleSizeSignificance(zo = 2.5, power = 0.8, level = 0.05, designPrior = "conditional")
sampleSizeSignificance(zo = 2.5, power = 0.8, level = 0.05, designPrior = "predictive")
sampleSizeSignificance(zo = 2.5, power = 0.8, level = 0.05, designPrior = "conditional",
                       shrinkage = 0.25)

## ----"plot-powerSignificance", echo = FALSE, fig.height = 4-------------------
po <- seq(0.0001, 0.05, 0.0001)/2

## plot power
plot(po, powerSignificance(zo = p2z(po, alternative = "one.sided"),
                           designPrior = "conditional")*100,
     type = "l", ylim = c(0, 100), lwd = 1.5, ylab = "Power (%)",
     xlab = expression(italic(p)[o]), las = 1, yaxt = "n")
axis(side = 2, at = seq(0, 100, 25), las = 1)
axis(side = 3, at = seq(0.0, 0.025, by = 0.005),
     labels = c("", round(p2z(p = seq(0.005, 0.025, by = 0.005),
                              alternative = "one.sided"), 2)))
mtext(text = expression(paste(italic(z)[o])), side = 3, line = 2)
abline(h = seq(0, 100, 25), lty = 1, col = adjustcolor("lightgrey", alpha.f = 0.3))
# abline(h = 50, lty = 1, lwd = 1.5, col = adjustcolor("lightgrey", alpha.f = 0.4))
# abline(h = 50, col = "#333333B3", lty = 3)
lines(po, powerSignificance(zo = p2z(po, alternative = "one.sided"),
                            designPrior = "predictive")*100, lwd = 2, lty = 2)
legend("topright", legend = c("conditional", "predictive"),
       title = "Design prior", lty = c(1, 2), lwd = 1.5, bty = "n")

## ----"plot-sampleSizeSignificance", echo = FALSE, fig.height = 4--------------
## plot sample size
plot(po, sampleSizeSignificance(zo = p2z(po, alternative = "one.sided"),
                                designPrior = "conditional", power = 0.8),
     type = "l", ylim = c(0.5, 8), log = "y", lwd = 1.5,
     ylab = expression(paste("Relative sample size ", n[r]/n[o])),
     xlab = expression(italic(p)[o]), las = 1, yaxt = "n")
axis(side = 3, at = seq(0.0, 0.025, by = 0.005),
     labels = c("", round(p2z(p = seq(0.005, 0.025, by = 0.005),
                              alternative = "one.sided"), 2)))
axis(side = 2, las = 1, at = c(0.5, 1, 2, 4, 8, 16, 32),
     labels = c("1/2", "1", "2", "4", "8", "16", "32"))
mtext(text = expression(paste(italic(z)[o])), side = 3, line = 2)
abline(h = c(0.5, 1, 2, 4, 8), lty = 1, col = adjustcolor("lightgrey", alpha.f = 0.3))
# abline(h = 1, col = "#333333B3", lty = 3)
abline(h = 1, lty = 1, lwd = 1.5, col = adjustcolor("lightgrey", alpha.f = 0.4))
lines(po, sampleSizeSignificance(zo = p2z(po, alternative = "one.sided"),
                                 designPrior = "predictive", power = 0.8),
      lwd = 2, lty = 2)
legend("topleft", legend = c("conditional", "predictive"),
       title = "Design prior", lty = c(1, 2), lwd = 1.5, bty = "n")

## ----"plot-predictionInterval", fig.height = 6--------------------------------
## compute prediction intervals for replication projects
par(mfrow = c(2, 2), las = 1, mai = rep(0.65, 4))
for (p in unique(RProjects$project)) {
    data_project <- subset(RProjects, project == p)
    pQ <- Qtest(thetao = data_project$fiso,
                thetar = data_project$fisr,
                seo = data_project$se_fiso,
                ser = data_project$se_fisr)
    PI <- predictionInterval(thetao = data_project$fiso,
                             seo = data_project$se_fiso,
                             ser = data_project$se_fisr)
    ## transforming back to correlation scale
    PI <- tanh(PI)
    incompatible <- pQ < 0.05 ## incompatible at 5% level
    color <- ifelse(incompatible == FALSE, "#00000099", "#DF536BCC")
    study <- seq(1, nrow(data_project))
    plot(data_project$rr, study, col = color, pch = 20, cex = 0.5,
         xlim = c(-0.5, 1), xlab = expression(italic(r)[r]), ylab = "Study",
         main = paste0(p, ": ", round(mean(incompatible)*100, 0), "% incompatible"))
    arrows(PI$lower, study, PI$upper, study, length = 0.02,
           angle = 90, code = 3, col = color)
    abline(v = 0, lty = 3)
}

## ----example-morewedge, echo = FALSE------------------------------------------
## data from study
morewedge <- subset(RProjects, study == "Morewedge et al. (2010), Science")

## Functions
ssv <- function(zo, so, a) {
    tau2 <- so^2/(zo^2/qnorm(p = a/2, lower.tail = FALSE)^2 - 1)
    return(tau2)
}

## Parameters and computations
optionsOld <- options(scipen = 5)
theta_o <- morewedge$fiso
theta_r <- morewedge$fisr
so <- morewedge$se_fiso
sr <- morewedge$se_fisr
c <- so^2/sr^2
zo <- theta_o/so
zr <- theta_r/sr
alpha <- 0.05
za <- qnorm(p = alpha/2, lower.tail = FALSE)
ps <- signif(pSceptical(zo = zo, zr = zr, c = c, alternative = "two.sided", type = "nominal"), 2)
ps1 <- signif(pSceptical(zo = zo, zr = zr, c = c,  alternative = "one.sided", type = "nominal"), 2)
ps1r <- signif(pSceptical(zo = zo, zr = zr, c = c,  alternative = "one.sided", type = "golden"), 2)
tau <- sqrt(ssv(zo, so, alpha)) # sufficiently sceptical prior sd at alpha = 0.05
s2_p <- 1/(1/so^2 + 1/tau^2) # posterior variance
mu_p <- s2_p*theta_o/so^2 # posterior mean

## ----"plot-pSceptical", echo = FALSE, fig.height = 4--------------------------

## Plot
df <- data.frame(estimate = factor(c("Original Study",
                                     # "Posterior",
                                     # "Sceptical Prior",
                                     "Replication Study"),
                                   levels = c("Original Study",
                                              # "Posterior",
                                              # "Sceptical Prior",
                                              "Replication Study")),
                 x = c(1,  2),
                 col = c(1, 4),
                 theta = c(theta_o,
                           theta_r),
                 lower = c(theta_o - za*so,
                           theta_r - za*sr),
                 upper = c(theta_o + za*so,
                           theta_r + za*sr),
                 p = c(signif(2*pnorm(q = zo, lower.tail = FALSE), 1),
                       signif(2*pnorm(q = zr, lower.tail = FALSE), 2)))

ylims <- c(-0.1, 1)
xlims <- c(0.5, 2.5)
cex <- 1
plot(x = df$x, y = df$theta, pch = " ",
     xlim = xlims, ylim = ylims,
     xlab = "", ylab = "Effect Size", xaxt = "n", las = 1)
axis(side = 1, at =  df$x, labels =  df$estimate,
     cex.axis = 0.9, mgp = c(3, 2, 0))
abline(h = seq(-2, 2, 0.25), lty = 1, col = adjustcolor("lightgrey", alpha.f = 0.4))
abline(h = 0, lty = 2)

# ## one-sided CIs
# arrows(x0 = df$x[-3], x1 = df$x[-3], y0 = df$theta[-3] + 0.05, y1 = df$upper[-3],
#        col = df$col[-3], code = 2, angle = 90, length = 0, lwd = 1.5)
# points(x = df$x[-3], y = df$upper[-3], pch = 17, cex = 1.3*cex, col = df$col[-3])

## two-sided CI
arrows(x0 = df$x, x1 = df$x, y0 = df$theta + 0.05, y1 = df$upper,
       col = df$col, code = 2, angle = 90, length = 0.1, lwd = 1.5)
arrows(x0 = df$x, x1 = df$x, y0 = df$theta - 0.05, y1 = df$lower,
       col = df$col, code = 2, angle = 90, length = 0.1, lwd = 1.5)

## arrows and text
points(x = df$x, y = df$theta, pch = 20, cex = cex*1.5, col = df$col)

text(x = df$x[1] - 0.2, y = df$theta[1],
       labels = bquote(paste(hat(theta)[o], " = ", .(round(df$theta[1], 2)))))


text(x = df$x[2] + 0.2, y = df$theta[2],
       labels = bquote(paste(hat(theta)[r], " = ", .(round(df$theta[2], 2)))))

text(x = df$x[1], y = df$upper[1] + 0.1,
       labels = bquote(paste(p[o], " = ", .(round(morewedge$po1, 4)))))


text(x = df$x[2], df$upper[2] + 0.1,
       labels = bquote(paste(p[r], " = ", .(round(morewedge$pr1, 4)))))


text(x = 1.5, y = 0.3,
     labels = bquote(paste(p[S], " = ", .(ps1))))




# arrows(x0 = 2.1, x1 = 2.9, y0 = 0.9, y1 = 0.9, length = 0.05, lwd = 1.2, col = "#000000B2")
# text(x = 2.5, y = 0.98, labels = "reverse-Bayes", cex = 0.75, col = "#000000B2")
# arrows(x0 = 3.15, x1 = 4.5, y0 = 0.9, y1 = 0.9, length = 0.05, lwd = 1.2, col = "#000000B2")
# text(x = 3.8, y = 0.98, labels = "prior-data conflict assessment", cex = 0.75, col = "#000000B2")
# points(x = 2, y = 0, pch = 1, cex = 2)
# arrows(x0 = 1.85, x1 = 1.98, y0 = -0.2, y1 = -0.03, length = 0.05, col = "#000000B2")
# text(x = 1.85, y = -0.25, labels = "fixed at zero", col = "#000000B2", cex = 0.75)
box()


## ----morewedge, echo = TRUE---------------------------------------------------
morewedge <- subset(RProjects, study == "Morewedge et al. (2010), Science")

## ----echo = FALSE-------------------------------------------------------------
ps_gold <- signif(pSceptical(zo = zo, zr = zr, c = c,  alternative = "one.sided", type = "golden"), 2)
ps_contr <- signif(pSceptical(zo = zo, zr = zr, c = c,  alternative = "one.sided", type = "controlled"), 2)

## ----echo = TRUE--------------------------------------------------------------
print(pS_nominal <- pSceptical(zo = morewedge$zo, zr = morewedge$zr,
                          c = morewedge$c,  alternative = "one.sided",
                          type = "nominal"))
print(pS_golden <- pSceptical(zo = morewedge$zo, zr = morewedge$zr,
                         c = morewedge$c,  alternative = "one.sided",
                         type = "golden"))
print(pS_controlled <- pSceptical(zo = morewedge$zo, zr = morewedge$zr,
                             c = morewedge$c,  alternative = "one.sided",
                             type = "controlled"))

## ----"plot-pSceptical-projects", fig.height = 4-------------------------------
## computing one-sided golden and controlled sceptical p-value for replication projects
RProjects$psG <- with(RProjects,
                     pSceptical(zo = zo, zr = zr, c = c,
                                alternative = "one.sided", type = "golden"))
RProjects$psC <- with(RProjects,
                     pSceptical(zo = zo, zr = zr, c = c,
                                alternative = "one.sided", type = "controlled"))

## ----echo = FALSE-------------------------------------------------------------

for (p in unique(RProjects$project)) {
    data_project <- subset(RProjects, project == p)
    cat(paste0(p, ": \n"))
    success_sceptG <- (data_project$psG < 0.025)
    cat(paste0(round(mean(success_sceptG)*100, 2),
               "% smaller than 0.025 (one-sided golden sceptical p-value) \n"))
    success_sceptC <- (data_project$psC < 0.025)
    cat(paste0(round(mean(success_sceptC)*100, 2),
               "% smaller than 0.025 (one-sided controlled sceptical p-value) \n"))
    success_tradit <- (data_project$po1 < 0.025) & (data_project$pr1 < 0.025)
    cat(paste0(round(mean(success_tradit)*100, 2),
               "% smaller than 0.025 (both one-sided traditional p-values) \n"))
    if(sum(success_sceptG != success_tradit) > 0){
        discrep <- data_project[(success_sceptG != success_tradit),
                                c("ro", "rr", "c", "po1", "pr1", "psG")]
        ## print effect estimates, 1sided p-values, and c of discrepant studies
        cat("Discrepant studies (golden vs significance): \n")
        print(signif(discrep, 2), row.names = FALSE)
    }

        if(sum(success_sceptC != success_tradit) > 0){
        discrep <- data_project[(success_sceptC != success_tradit),
                                c("ro", "rr", "c", "po1", "pr1", "psC")]
        ## print effect estimates, 1sided p-values, and c of discrepant studies
        cat("Discrepant studies (controlled vs significance): \n")
        print(signif(discrep, 2), row.names = FALSE)
        }

        if(sum(success_sceptC != success_sceptG) > 0){
        discrep <- data_project[(success_sceptC != success_sceptG),
                                c("ro", "rr", "c", "po1", "pr1", "psG", "psC")]
        ## print effect estimates, 1sided p-values, and c of discrepant studies
        cat("Discrepant studies (golden vs controlled): \n")
        print(signif(discrep, 2), row.names = FALSE)
  }
  cat("\n \n")
}

## ----golden-vs-controlled, fig.height = 6-------------------------------------
par(mfrow = c(2, 2), las = 1, pty = "s",
    mar = c(4, 2.5, 2.5, 1))
  myaxis <- c(10e-5, 0.001, 0.01, 0.1, 1)
for (p in unique(RProjects$project)) {
  data_project <- subset(RProjects, project == p)
  plot(psG ~ psC, data = data_project, ylim = c(10e-5, 1),
       xlim = c(10e-5, 1), main = p, xlab = expression(paste("controlled ", p[S])),
       ylab = expression(paste("golden ", p[S])),
       pch = 19,
       col = "#00000099",
       cex = 1.3,
       axes = F,
       log = "xy")
  abline(h = 0, lty = 2)
  abline(a = 0, b = 1, col = "grey")
  abline(h = 0.025, lty = 2, col = "grey")
  abline(v = 0.025, lty = 2, col = "grey")

axis(1, at = myaxis, labels = myaxis)
axis(2,at = myaxis, labels = myaxis)
box()
}


## ----"threshold-p-sceptical2", echo = FALSE-----------------------------------
## computing nominal, controlled, and golden replication success levels 
## for one-sided uncalibrated sceptical p-value
thresh_gol <- levelSceptical(level = 0.025, alternative = "one.sided",
                            type = "golden")
thresh_contr <- levelSceptical(level = 0.025, alternative = "one.sided",
                              type = "controlled", c = c(1, 2))
thresh_nom <- levelSceptical(level = 0.025, alternative = "one.sided",
                            type = "nominal")

## ----"threshold-p-sceptical"--------------------------------------------------
## computing nominal, golden and controlled replication success levels
## for one-sided uncalibrated sceptical p-value

print(rs_level_nom <- levelSceptical(level = 0.025, alternative = "one.sided",
                              type = "nominal"))
print(rs_level_gol <- levelSceptical(level = 0.025, alternative = "one.sided",
                              type = "golden"))
print(rs_level_contr <- levelSceptical(level = 0.025, alternative = "one.sided",
                                type = "controlled", c = c(1, 2)))


## -----------------------------------------------------------------------------
sampleSizeReplicationSuccess(zo = 2.5, power = 0.8, level = 0.025,
                             alternative = "one.sided",
                             designPrior = "conditional",
                             type = c("golden", "controlled"))

sampleSizeReplicationSuccess(zo = 2.5, power = 0.8, level = 0.025,
                             alternative = "one.sided",
                             designPrior = "predictive",
                             type = c("golden", "controlled"))

## ----"plot-powerReplicationSuccess", echo = FALSE, fig.height = 4-------------
## plot power
po <- seq(0.001, 0.05, length.out = 100)/2
plot(po, powerReplicationSuccess(zo = p2z(po, alternative  = "one.sided"),
                                 designPrior = "conditional",
                                 level = 0.025,
                                 alternative = "one.sided",
                                 type = "golden")*100,
     type = "l", ylim = c(0, 100), lwd = 1.5, ylab = "Power (%)", las = 1,
     xlab = expression(italic(p)[o]), yaxt = "n")
abline(h = seq(0, 100, 25), lty = 1, col = adjustcolor("lightgrey", alpha.f = 0.3))
# abline(h = 50, lty = 1, lwd = 1.5, col = adjustcolor("lightgrey", alpha.f = 0.4))
axis(side = 2, at = seq(0, 100, 25), las = 1)
axis(side = 3, at = seq(0.0, 0.025, by = 0.005),
     labels = c("", round(p2z(p = seq(0.005, 0.025, by = 0.005),
                              alternative = "one.sided"), 2)))
mtext(text = expression(paste( italic(z)[o])), side = 3, line = 2)
lines(po, powerReplicationSuccess(zo = p2z(po, alternative = "one.sided"),
                                  designPrior = "predictive",
                                  level = 0.025,
                                  alternative = "one.sided",
                                  type = "golden")*100,
      lwd = 2, lty = 2)
lines(po, powerReplicationSuccess(zo = p2z(po, alternative  = "one.sided"),
                                 designPrior = "conditional",
                                 level = 0.025,
                                 alternative = "one.sided",
                                 type = "controlled")*100,
      col = "red", lwd = 2)
lines(po, powerReplicationSuccess(zo = p2z(po, alternative  = "one.sided"),
                                 designPrior = "predictive",
                                 level = 0.025,
                                 alternative = "one.sided",
                                 type = "controlled")*100,
      col = "red", lty = 2, lwd = 2)

legend("topright", legend = c("conditional", "predictive"),
       title = "Design prior", lty = c(1, 2), lwd = 1.5, bty = "n")
# abline(h = 50, lty = 3)

legend("bottomleft",
       legend = c(expression(paste("golden ", p[S])),
                  expression(paste("controlled ", p[S]))),
       col = c("black", "red"),
       lty = 1, lwd = 2,
       bty = "n")

## ----"plot-sampleSizeReplicationSuccess", echo = FALSE, fig.height = 4--------
# po <- seq(0.0001, 0.05, 0.0001)
po <- seq(0.001, 0.05, length.out = 100)/2

## plot sample size
plot(po,
     sampleSizeReplicationSuccess(zo = p2z(po, alternative = "one.sided"),
                                  power = 0.8,
                                  level = 0.025,
                                  designPrior = "conditional",
                                  alternative = "one.sided",
                                  type = "golden"),
     type = "l", log = "y", lwd = 1.5, las = 1, ylim = c(0.5, 20),
     ylab = expression(paste("Relative sample size ", n[r]/n[o])),
     xlab = expression(italic(p)[o]), yaxt = "n")
axis(side = 2, las = 1, at = c(0.5, 1, 2, 4, 8, 16, 32),
     labels = c("1/2", "1", "2", "4", "8", "16", "32"))
# abline(h = 1, lty = 3)
abline(h = c(0.5, 1, 2, 4, 8, 16, 32), lty = 1, col = adjustcolor("lightgrey", alpha.f = 0.3))
abline(h = 1, lty = 1, lwd = 1.5, col = adjustcolor("lightgrey", alpha.f = 0.4))
axis(side = 3, at = seq(0.0, 0.025, by = 0.005),
     labels = c("", round(p2z(p = seq(0.005, 0.025, by = 0.005),
                              alternative = "one.sided"), 2)))
mtext(text = expression(paste(italic(z)[o])), side = 3, line = 2)
suppressWarnings({
    lines(po, sampleSizeReplicationSuccess(zo = p2z(po, alternative = "one.sided"),
                                           power = 0.8, level = 0.025,
                                           designPrior = "predictive",
                                           alternative = "one.sided",
                                           type = "golden"),
          lwd = 2, lty = 2)

  lines(po,
     sampleSizeReplicationSuccess(zo = p2z(po, alternative = "one.sided"),
                                  power = 0.8,
                                  level = 0.025,
                                  designPrior = "conditional",
                                  alternative = "one.sided",
                                  type = "controlled"),
     col = "red", lwd = 2)

    lines(po,
     sampleSizeReplicationSuccess(zo = p2z(po, alternative = "one.sided"),
                                  power = 0.8,
                                  level = 0.025,
                                  designPrior = "predictive",
                                  alternative = "one.sided",
                                  type = "controlled"),
     col = "red", lwd = 2, lty = 2)
})
legend("topleft", legend = c("conditional", "predictive"),
       title = "Design prior", lty = c(1, 2), lwd = 1.5, bty = "n")

legend("bottomright",
       legend = c(expression(paste("golden ", p[S])),
                  expression(paste("controlled ", p[S]))),
       col = c("black", "red"),
       lty = 1, lwd = 2,
       bty = "n")

## ----echo = FALSE-------------------------------------------------------------
d <- morewedge$fisr/morewedge$fiso
dminrs <- effectSizeReplicationSuccess(zo = morewedge$zo, c = morewedge$c)
dminsign <- effectSizeSignificance(zo = morewedge$zo, c = morewedge$c)

## ----echo = FALSE, fig.height = 3.5, eval = FALSE-----------------------------
#  data("SSRP")
#  intpow_cond <- with(SSRP, powerSignificanceInterim(zo = fiso/se_fiso, zi = fisi/se_fisi, c = nr/no,
#                                                     f = ni/nr, designPrior = "conditional"))
#  intpow_pred <- with(SSRP, powerSignificanceInterim(zo = fiso/se_fiso, zi = fisi/se_fisi, c = nr/no,
#                                                     f = ni/nr, designPrior = "predictive"))
#  plot(intpow_cond*100, intpow_pred*100,
#       xlab = "Conditional power (in %)",
#       ylab = "Predictive power (in %)",
#       pch = 20,
#       cex = 1.2,
#       xlim = c(80, 100),
#       ylim = c(0, 100),
#       yaxt = "n", las = 1,
#       col = "#00000099")
#  axis(side = 2, at = seq(0, 100, 25), las = 1)
#  abline(a = 0, b = 1, col = "grey")

## ----"shrinkage", echo = FALSE, fig.height = 3.5------------------------------
zo <- seq(0, 4, 0.01)
s <- pmax(1 - 1/zo^2, 0)
shrinkage <- (1 - s)*100
plot(zo, shrinkage, type = "l", ylim = c(0, 100), las = 1,
     xlab = expression(paste(italic(z)[o])), ylab = "Shrinkage (%)",
     yaxt = "n")
axis(side = 2, at = seq(0, 100, 25), las = 1)
axis(side = 3, at = seq(0, 4, by = 1),
     labels = c(signif(z2p(seq(0, 3, by = 1), alternative = "one.sided"), 2),
                signif(z2p(4, alternative = "one.sided"), 1)))
abline(h = seq(0, 100, 25), lty = 1, col = adjustcolor("lightgrey", alpha.f = 0.3))
mtext(text = expression(italic(p)[o]), side = 3, line = 2)

## ----include=FALSE------------------------------------------------------------
options(optionsOld)
par(parOld)

