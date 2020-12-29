
# ANOVA main type I ------------------------------------------------------------

#' ANOVA type I for main effects
#'
#' @inheritParams univariate
#'
#' @keywords internal
anova_main_I <- function(x, es_aov, digits, CI, lower.tail) {
  M <- x$M.mu
  F <- x$F.mu
  d <- M - F
  nM <- x$m
  nF <- x$f
  N <- nM + nF
  M.sd <- x$M.sdev
  F.sd <- x$F.sdev
  p <- 1
  r <- NROW(M)
  o.p <- rep(1, p)
  o.r <- rep(1, r)

  Xm <- nM %*% t(o.p) * M
  Xf <- nF %*% t(o.p) * F
  N <- sum(nM) + sum(nF)
  Gmean <- t(as.numeric(t(o.r) %*% (Xm + Xf) / N))
  Mc <- M - o.r %*% Gmean
  Fc <- F - o.r %*% Gmean
  J <- matrix(1, nrow = r, ncol = r)
  D <- M - F
  w <- (nM * nF) / (nM + nF)
  weighted.D <- as.numeric(t(D) %*% w)

  SS.b <- as.numeric(t(Mc) %*% diag(nM) %*% Mc + t(Fc) %*% diag(nF) %*% Fc)
  SS.i <- as.numeric(t(D) %*% (w %*% t(o.p) * D) - weighted.D %o% weighted.D / sum(w))
  SS.e <- as.numeric(sum((nM - 1) * M.sd^2 + (nF - 1) * F.sd^2))
  SS.samp <- as.numeric(t(Xm) %*% M + t(Xf) %*% F - t(Xm) %*% J %*% Xm / sum(nM) - t(Xf) %*% J %*% Xf / sum(nF) - SS.i)
  SS.sex <- as.numeric(SS.b - SS.samp - SS.i)

  SS.e <- SS.e + SS.i

  SS <- c(SS.sex, SS.samp, SS.e)
  Df <- c(1, r - 1, sum(nM) + sum(nF) - r - 1)
  MSQ <- SS / Df
  SS <- round(SS, 1)
  F <- MSQ[-3] / MSQ[3]
  MSQ <- round(MSQ, 1)

  matrix(data = NA, nrow = 2, ncol = 3, byrow = T) -> eff_all
  rep(NA, 2) -> upper
  rep(NA, 2) -> lower
  eta <- rep(NA, 2)
  omega <- rep(NA, 2)
  cohen_f <- rep(NA, 2)
  prob <- rep(NA, 2)
  for (i in 1:2) {
    prob[i] <- pf(F[i], Df[i], Df[3],lower.tail = lower.tail)
    eta[i] <- SS[i] / (SS[i] + SS[3])
    cohen_f[i] <- (eta[i]) / (1 - eta[i])
    omega[i] <- (Df[i] * (MSQ[i] - MSQ[3])) / ((Df[i] * MSQ[i]) + (N - Df[i]) * MSQ[3])
  }
  F <- c(round(F, digits), NA)
  prob <- round(prob, digits)
  prob[3] <- NA
  eff <- switch(es_aov, f = cohen_f, eta = eta, omega = omega, none = eta)
  for (i in 1:2) {
    eff_all[i, ] <- as.matrix(eff_CI(f = F[i], CI = CI, eff = eff[i], df1 = Df[i], df2 = Df[3]))
    lower [i] <- eff_all[i, 2]
    upper [i] <- eff_all[i, 3]
  }
  term <- c("Sex", "Pop", "Residuals")
  sto <- cbind_fill(
    term = term, df = Df, sumsq = round(SS, digits), meansq = round(MSQ, digits),
    statistic = F, p.value = prob, round(eff, digits), lower = round(
      lower,
      digits
    ), upper = round(upper, digits)
  )
  names(sto)[7] <- es_aov
  sto <- add_sig(x = sto)
  names(sto)[9] <- paste0(names(sto)[9], ".", es_aov)
  names(sto)[10] <- paste0(names(sto)[10], ".", es_aov)
  if (es_aov == "none") {
    sto <- sto[1:7]
  }

  SS.sex <- as.numeric(t(d) %*% w %*% t(w) %*% d / as.numeric(t(o.r) %*% w))
  SS.samp <- SS.b - SS.sex - SS.i

  SS <- c(SS.samp, SS.sex, SS.e)
  Df[1] <- r - 1
  Df[2] <- 1
  MSQ <- SS / Df
  SS <- round(SS, 1)
  F <- MSQ[-3] / MSQ[3]
  MSQ <- round(MSQ, 1)

  matrix(data = NA, nrow = 2, ncol = 3, byrow = T) -> eff_all
  rep(NA, 2) -> upper
  rep(NA, 2) -> lower
  eta <- rep(NA, 2)
  omega <- rep(NA, 2)
  cohen_f <- rep(NA, 2)
  prob <- rep(NA, 2)
  for (i in 1:2) {
    prob[i] <- pf(F[i], Df[i], Df[3],lower.tail = lower.tail)
    eta[i] <- SS[i] / (SS[i] + SS[3])
    cohen_f[i] <- (eta[i]) / (1 - eta[i])
    omega[i] <- (Df[i] * (MSQ[i] - MSQ[3])) / ((Df[i] * MSQ[i]) + (N - Df[i]) * MSQ[3])
  }
  F <- c(round(F, digits), NA)
  prob <- round(prob, digits)
  prob[3] <- NA
  eff <- switch(es_aov, f = cohen_f, eta = eta, omega = omega, none = eta)


  for (i in 1:2) {
    eff_all[i, ] <- as.matrix(eff_CI(f = F[i], CI = CI, eff = eff[i], df1 = Df[i], df2 = Df[3]))
    lower [i] <- eff_all[i, 2]
    upper [i] <- eff_all[i, 3]
  }


  term <- c("Pop", "Sex", "Residuals")
  sto2 <- cbind_fill(
    term = term, df = Df, sumsq = round(SS, digits), meansq = round(MSQ, digits),
    statistic = F, p.value = prob, round(eff, digits), lower = round(
      lower,
      digits
    ), upper = round(upper, digits)
  )
  names(sto2)[7] <- es_aov
  sto2 <- add_sig(x = sto2)
  names(sto2)[9] <- paste0(names(sto2)[9], ".", es_aov)
  names(sto2)[10] <- paste0(names(sto2)[10], ".", es_aov)
  if (es_aov == "none") {
    sto2 <- sto2[1:7]
  }
  list(sto, sto2)
}


# ANOVA main type II ------------------------------------------------------

#' ANOVA type II for main effects
#'
#' @inheritParams univariate
#'
#' @keywords internal
#'
anova_main_II <- function(x, es_aov, digits, CI, lower.tail) {
  M <- x$M.mu
  F <- x$F.mu
  nM <- x$m
  nF <- x$f
  M.sd <- x$M.sdev
  F.sd <- x$F.sdev
  p <- 1
  r <- NROW(M)
  o.p <- rep(1, p)
  o.r <- rep(1, r)

  Xm <- nM %*% t(o.p) * M
  Xf <- nF %*% t(o.p) * F
  N <- sum(nM) + sum(nF)
  J <- matrix(1, nrow = r, ncol = r)
  D <- M - F
  w <- (nM * nF) / (nM + nF)
  weighted.D <- as.numeric(t(D) %*% w)

  SS.i <- as.numeric(t(D) %*% (w %*% t(o.p) * D) - weighted.D %o% weighted.D / sum(w))
  SS.e <- as.numeric(sum((nM - 1) * M.sd^2 + (nF - 1) * F.sd^2))
  SS.samp <- as.numeric(t(Xm) %*% M + t(Xf) %*% F - t(Xm) %*% J %*% Xm / sum(nM) - t(Xf) %*% J %*% Xf / sum(nF) - SS.i)
  SS.sex <- as.numeric(t(D) %*% w %*% t(w) %*% D / as.numeric(t(o.r) %*% w))
  SS.e <- SS.e + SS.i

  SS <- c(SS.sex, SS.samp, SS.e)
  Df <- c(1, r - 1, sum(nM) + sum(nF) - r - 1)
  MSQ <- SS / Df
  SS <- round(SS, 1)
  F <- MSQ[-3] / MSQ[3]
  MSQ <- round(MSQ, 1)

  matrix(data = NA, nrow = 2, ncol = 3, byrow = T) -> eff_all
  rep(NA, 2) -> upper
  rep(NA, 2) -> lower
  eta <- rep(NA, 2)
  omega <- rep(NA, 2)
  cohen_f <- rep(NA, 2)
  prob <- rep(NA, 2)
  for (i in 1:2) {
    prob[i] <- pf(F[i], Df[i], Df[3],lower.tail = lower.tail)
    eta[i] <- SS[i] / (SS[i] + SS[3])
    cohen_f[i] <- (eta[i]) / (1 - eta[i])
    omega[i] <- (Df[i] * (MSQ[i] - MSQ[3])) / ((Df[i] * MSQ[i]) + (N - Df[i]) * MSQ[3])
  }
  F <- c(round(F, digits), NA)
  prob <- round(prob, digits)
  prob[3] <- NA
  eff <- switch(es_aov, f = cohen_f, eta = eta, omega = omega, none = eta)
  for (i in 1:2) {
    eff_all[i, ] <- as.matrix(eff_CI(f = F[i], CI = CI, eff = eff[i], df1 = Df[i], df2 = Df[3]))
    lower [i] <- eff_all[i, 2]
    upper [i] <- eff_all[i, 3]
  }
  term <- c("Sex", "Pop", "Residuals")
  sto <- cbind_fill(
    term = term, df = Df, sumsq = round(SS, digits), meansq = round(MSQ, digits),
    statistic = F, p.value = prob, round(eff, digits), lower = round(
      lower,
      digits
    ), upper = round(upper, digits)
  )
  names(sto)[7] <- es_aov
  sto <- add_sig(x = sto)
  names(sto)[9] <- paste0(names(sto)[9], ".", es_aov)
  names(sto)[10] <- paste0(names(sto)[10], ".", es_aov)
  if (es_aov == "none") {
    sto <- sto[1:7]
  }
  sto
}


# ANOVA type I ----------------------------------------------------------

#' ANOVA type I
#'
#' @inheritParams univariate
#'
#' @keywords internal
#'
anova_I <- function(x, es_aov, digits, CI, lower.tail) {
  M <- x$M.mu
  F <- x$F.mu
  d <- M - F
  nM <- x$m
  nF <- x$f
  M.sd <- x$M.sdev
  F.sd <- x$F.sdev
  p <- 1
  r <- NROW(M)
  o.p <- rep(1, p)
  o.r <- rep(1, r)

  Xm <- nM %*% t(o.p) * M
  Xf <- nF %*% t(o.p) * F
  N <- sum(nM) + sum(nF)
  Gmean <- as.numeric(t(o.r) %*% (Xm + Xf) / N)
  Mc <- M - rep(Gmean, r)
  Fc <- F - rep(Gmean, r)
  J <- matrix(1, nrow = r, ncol = r)
  D <- M - F
  w <- (nM * nF) / (nM + nF)
  weighted.D <- as.numeric(t(D) %*% w)

  SS.b <- as.numeric(t(Mc) %*% diag(nM) %*% Mc + t(Fc) %*% diag(nF) %*% Fc)
  SS.i <- as.numeric(t(D) %*% (w %*% t(o.p) * D) - weighted.D %o% weighted.D / sum(w))
  SS.e <- as.numeric(sum((nM - 1) * M.sd^2 + (nF - 1) * F.sd^2))
  SS.samp <- as.numeric(t(Xm) %*% M + t(Xf) %*% F - t(Xm) %*% J %*% Xm / sum(nM) - t(Xf) %*% J %*% Xf / sum(nF) - SS.i)
  SS.sex <- as.numeric(SS.b - SS.samp - SS.i)

  SS <- c(SS.sex, SS.samp, SS.i, SS.e)
  Df <- c(1, r - 1, r - 1, sum(nF) + sum(nM) - 2 * r)
  MSQ <- SS / Df
  SS <- round(SS, 1)
  F <- MSQ[-4] / MSQ[4]
  MSQ <- round(MSQ, 1)

  matrix(data = NA, nrow = 3, ncol = 3, byrow = T) -> eff_all
  rep(NA, 3) -> upper
  rep(NA, 3) -> lower
  eta <- rep(NA, 3)
  omega <- rep(NA, 3)
  cohen_f <- rep(NA, 3)
  prob <- rep(NA, 3)
  for (i in 1:3) {
    prob[i] <-  pf(F[i], Df[i], Df[4],lower.tail = lower.tail)
    eta[i] <- SS[i] / (SS[i] + SS[4])
    cohen_f[i] <- (eta[i]) / (1 - eta[i])
    omega[i] <- (Df[i] * (MSQ[i] - MSQ[4])) / ((Df[i] * MSQ[i]) + (N - Df[i]) * MSQ[4])
  }
  F <- c(round(F, digits), NA)
  prob <- round(prob, digits)
  prob[4] <- NA

  eff <- switch(es_aov, f = cohen_f, eta = eta, omega = omega, none = eta)


  for (i in 1:3) {
    eff_all[i, ] <- as.matrix(eff_CI(f = F[i], CI = CI, eff = eff[i], df1 = Df[i], df2 = Df[4]))
    lower [i] <- eff_all[i, 2]
    upper [i] <- eff_all[i, 3]
  }

  term <- c("Sex", "Pop", "Sex*Pop", "Residuals")
  sto <- cbind_fill(
    term = term, df = Df, sumsq = round(SS, digits), meansq = round(MSQ, digits),
    statistic = F, p.value = prob, round(eff, digits), lower = round(
      lower,
      digits
    ), upper = round(upper, digits)
  )
  names(sto)[7] <- es_aov
  sto <- add_sig(x = sto)
  names(sto)[9] <- paste0(names(sto)[9], ".", es_aov)
  names(sto)[10] <- paste0(names(sto)[10], ".", es_aov)
  if (es_aov == "none") {
    sto <- sto[1:7]
  }

  SS.sex <- as.numeric(t(d) %*% w %*% t(w) %*% d / as.numeric(t(o.r) %*% w))
  SS.samp <- SS.b - SS.sex - SS.i

  SS <- c(SS.samp, SS.sex, SS.i, SS.e)
  Df[1] <- r - 1
  Df[2] <- 1

  MSQ <- SS / Df
  SS <- round(SS, 1)
  F <- MSQ[-4] / MSQ[4]
  MSQ <- round(MSQ, 1)

  matrix(data = NA, nrow = 3, ncol = 3, byrow = T) -> eff_all
  rep(NA, 3) -> upper
  rep(NA, 3) -> lower
  eta <- rep(NA, 3)
  omega <- rep(NA, 3)
  cohen_f <- rep(NA, 3)
  prob <- rep(NA, 3)
  for (i in 1:3) {
    prob[i] <-  pf(F[i], Df[i], Df[4],lower.tail = lower.tail)
    eta[i] <- SS[i] / (SS[i] + SS[4])
    cohen_f[i] <- (eta[i]) / (1 - eta[i])
    omega[i] <- (Df[i] * (MSQ[i] - MSQ[4])) / ((Df[i] * MSQ[i]) + (N - Df[i]) * MSQ[4])
  }
  F <- c(round(F, digits), NA)
  prob <- round(prob, digits)
  prob[4] <- NA

  eff <- switch(es_aov, f = cohen_f, eta = eta, omega = omega, none = eta)


  for (i in 1:3) {
    eff_all[i, ] <- as.matrix(eff_CI(f = F[i], CI = CI, eff = eff[i], df1 = Df[i], df2 = Df[4]))
    lower [i] <- eff_all[i, 2]
    upper [i] <- eff_all[i, 3]
  }

  term <- c("Pop", "Sex", "Sex*Pop", "Residuals")
  sto2 <- cbind_fill(
    term = term, df = Df, sumsq = round(SS, digits), meansq = round(MSQ, digits),
    statistic = F, p.value = prob, round(eff, digits), lower = round(
      lower,
      digits
    ), upper = round(upper, digits)
  )
  names(sto2)[7] <- es_aov
  sto2 <- add_sig(x = sto2)
  names(sto2)[9] <- paste0(names(sto2)[9], ".", es_aov)
  names(sto2)[10] <- paste0(names(sto2)[10], ".", es_aov)
  if (es_aov == "none") {
    sto2 <- sto2[1:7]
  }

  list(sto, sto2)
}


# ANOVA type II -----------------------------------------------------------

#' ANOVA type II
#'
#' @inheritParams univariate
#'
#' @keywords internal
#'
anova_II <- function(x, es_aov, digits, CI, lower.tail) {
  M <- x$M.mu
  F <- x$F.mu
  d <- M - F
  nM <- x$m
  nF <- x$f
  M.sd <- x$M.sdev
  F.sd <- x$F.sdev
  p <- 1
  r <- NROW(M)
  o.p <- rep(1, p)
  o.r <- rep(1, r)

  Xm <- nM %*% t(o.p) * M
  Xf <- nF %*% t(o.p) * F
  N <- sum(nM) + sum(nF)
  Gmean <- as.numeric(t(o.r) %*% (Xm + Xf) / N)
  Mc <- M - rep(Gmean, r)
  Fc <- F - rep(Gmean, r)
  J <- matrix(1, nrow = r, ncol = r)
  D <- M - F
  w <- (nM * nF) / (nM + nF)
  weighted.D <- as.numeric(t(D) %*% w)

  SS.i <- as.numeric(t(D) %*% (w %*% t(o.p) * D) - weighted.D %o% weighted.D / sum(w))
  SS.e <- as.numeric(sum((nM - 1) * M.sd^2 + (nF - 1) * F.sd^2))
  SS.samp <- as.numeric(t(Xm) %*% M + t(Xf) %*% F - t(Xm) %*% J %*% Xm / sum(nM) - t(Xf) %*% J %*% Xf / sum(nF) - SS.i)
  SS.sex <- as.numeric(t(D) %*% w %*% t(w) %*% D / as.numeric(t(o.r) %*% w))


  SS <- c(SS.sex, SS.samp, SS.i, SS.e)
  Df <- c(1, r - 1, r - 1, sum(nF) + sum(nM) - 2 * r)
  MSQ <- SS / Df
  SS <- round(SS, 1)
  F <- MSQ[-4] / MSQ[4]
  MSQ <- round(MSQ, 1)

  matrix(data = NA, nrow = 3, ncol = 3, byrow = T) -> eff_all
  rep(NA, 3) -> upper
  rep(NA, 3) -> lower
  eta <- rep(NA, 3)
  omega <- rep(NA, 3)
  cohen_f <- rep(NA, 3)
  prob <- rep(NA, 3)
  for (i in 1:3) {
    prob[i] <-  pf(F[i], Df[i], Df[4],lower.tail = lower.tail)
    eta[i] <- SS[i] / (SS[i] + SS[4])
    cohen_f[i] <- (eta[i]) / (1 - eta[i])
    omega[i] <- (Df[i] * (MSQ[i] - MSQ[4])) / ((Df[i] * MSQ[i]) + (N - Df[i]) * MSQ[4])
  }
  F <- c(round(F, digits), NA)
  prob <- round(prob, digits)
  prob[4] <- NA

  eff <- switch(es_aov, f = cohen_f, eta = eta, omega = omega, none = eta)


  for (i in 1:3) {
    eff_all[i, ] <- as.matrix(eff_CI(f = F[i], CI = CI, eff = eff[i], df1 = Df[i], df2 = Df[4]))
    lower [i] <- eff_all[i, 2]
    upper [i] <- eff_all[i, 3]
  }

  term <- c("Sex", "Pop", "Sex*Pop", "Residuals")
  sto <- cbind_fill(
    term = term, df = Df, sumsq = round(SS, digits), meansq = round(MSQ, digits),
    statistic = F, p.value = prob, round(eff, digits), lower = round(
      lower,
      digits
    ), upper = round(upper, digits)
  )
  names(sto)[7] <- es_aov
  sto <- add_sig(x = sto)
  names(sto)[9] <- paste0(names(sto)[9], ".", es_aov)
  names(sto)[10] <- paste0(names(sto)[10], ".", es_aov)
  if (es_aov == "none") {
    sto <- sto[1:7]
  }
  sto
}


# ANOVA type III ----------------------------------------------------------

#' ANOVA type III
#'
#' @inheritParams univariate
#'
#' @keywords internal
#'
anova_III <- function(x, es_aov, digits, CI, lower.tail) {
  M <- x$M.mu
  F <- x$F.mu
  nM <- x$m
  nF <- x$f
  N.M <- sum(nM)
  N.F <- sum(nF)
  N <- N.M + N.F
  M.sd <- x$M.sdev
  F.sd <- x$F.sdev
  p <- 1
  r <- NROW(M)
  o.p <- rep(1, p)
  o.r <- rep(1, r)

  SS.e <- as.numeric(sum((nM - 1) * M.sd^2 + (nF - 1) * F.sd^2))

  Pop.Block <- t(contr.sum(r)) %*% diag(nM + nF) %*% contr.sum(r)

  contrast <- contr.sum(2) %x% contr.sum(r)
  first <- t(contrast) %*% diag(c(nF, nM))
  first <- first[, 1:r] + first[, (r + 1):(2 * r)]
  second <- contr.sum(r)
  Interact.Block <- first %*% second

  Xm <- nM %*% t(o.p) * M
  Xf <- nF %*% t(o.p) * F
  Male.sum <- apply(Xm, 2, sum)
  Female.sum <- apply(Xf, 2, sum)
  Pop.sums <- Xm + Xf
  Pops <- t(contr.sum(r)) %*% Pop.sums
  XY.interact <- matrix(0, nrow = r - 1, ncol = p)
  XX.1 <- vector()
  XX.2 <- vector()
  for (i in 1:(r - 1)) {
    XY.interact[i, ] <- Xf[i, ] - Xf[r, ] + Xm[r, ] - Xm[i, ]
    XX.1[i] <- nF[i] - nF[r] + nM[r] - nM[i]
    XX.2[i] <- (nF[i] + nM[i]) - (nF[r] + nM[r])
  }
  XY <- rbind(Female.sum + Male.sum, Female.sum - Male.sum, Pops, XY.interact)

  XX <- matrix(NA, ncol = 2 * r, nrow = 2 * r)
  XX[1:2, 1:2] <- diag(rep(N, 2))
  XX[1, 2] <- N.F - N.M
  XX[2, 1] <- XX[1, 2]
  XX[3:(2 * r), 1] <- c(XX.2, XX.1)
  XX[1, 3:(2 * r)] <- c(XX.2, XX.1)
  XX[3:(2 * r), 2] <- c(XX.1, XX.2)
  XX[2, 3:(2 * r)] <- c(XX.1, XX.2)
  XX[3:(1 + r), 3:(1 + r)] <- Pop.Block
  XX[(2 + r):(2 * r), (2 + r):(2 * r)] <- Pop.Block
  XX[(r + 2):(2 * r), 3:(1 + r)] <- Interact.Block
  XX[3:(1 + r), (r + 2):(2 * r)] <- Interact.Block

  inv.XX <- solve(XX)
  B <- inv.XX %*% XY

  A <- t(c(0, 1, rep(0, 2 * r - 2)))
  AB <- A %*% B
  inv.contr <- A %*% inv.XX %*% t(A)
  SS.sex <- as.numeric(t(AB) %*% solve(inv.contr) %*% AB)

  A <- cbind(matrix(0, nrow = r - 1, ncol = 2), diag(r - 1), matrix(0, nrow = r - 1, ncol = r - 1))
  AB <- A %*% B
  inv.contr <- A %*% inv.XX %*% t(A)
  SS.samp <- as.numeric(t(AB) %*% solve(inv.contr) %*% AB)


  A <- cbind(matrix(0, nrow = r - 1, ncol = 2), matrix(0, nrow = r - 1, ncol = r - 1), diag(r - 1))
  AB <- A %*% B
  inv.contr <- A %*% inv.XX %*% t(A)
  SS.i <- as.numeric(t(AB) %*% solve(inv.contr) %*% AB)

  SS <- c(SS.sex, SS.samp, SS.i, SS.e)
  Df <- c(1, r - 1, r - 1, sum(nF) + sum(nM) - 2 * r)
  MSQ <- SS / Df
  SS <- round(SS, 1)
  F <- MSQ[-4] / MSQ[4]
  MSQ <- round(MSQ, 1)

  matrix(data = NA, nrow = 3, ncol = 3, byrow = T) -> eff_all
  rep(NA, 3) -> upper
  rep(NA, 3) -> lower
  eta <- rep(NA, 3)
  omega <- rep(NA, 3)
  cohen_f <- rep(NA, 3)
  prob <- rep(NA, 3)
  for (i in 1:3) {
    prob[i] <-  pf(F[i], Df[i], Df[4],lower.tail = lower.tail)
    eta[i] <- SS[i] / (SS[i] + SS[4])
    cohen_f[i] <- (eta[i]) / (1 - eta[i])
    omega[i] <- (Df[i] * (MSQ[i] - MSQ[4])) / ((Df[i] * MSQ[i]) + (N - Df[i]) * MSQ[4])
  }
  F <- c(round(F, digits), NA)
  prob <- round(prob, digits)
  prob[4] <- NA

  eff <- switch(es_aov, f = cohen_f, eta = eta, omega = omega, none = eta)


  for (i in 1:3) {
    eff_all[i, ] <- as.matrix(eff_CI(f = F[i], CI = CI, eff = eff[i], df1 = Df[i], df2 = Df[4]))
    lower [i] <- eff_all[i, 2]
    upper [i] <- eff_all[i, 3]
  }

  term <- c("Sex", "Pop", "Sex*Pop", "Residuals")
  sto <- cbind_fill(
    term = term, df = Df, sumsq = round(SS, digits), meansq = round(MSQ, digits),
    statistic = F, p.value = prob, round(eff, digits), lower = round(
      lower,
      digits
    ), upper = round(upper, digits)
  )
  names(sto)[7] <- es_aov
  sto <- add_sig(x = sto)
  names(sto)[9] <- paste0(names(sto)[9], ".", es_aov)
  names(sto)[10] <- paste0(names(sto)[10], ".", es_aov)
  if (es_aov == "none") {
    sto <- sto[1:7]
  }
  sto
}
