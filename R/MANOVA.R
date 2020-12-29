# manova_main_I -----------------------------------------------------------

#' MANOVA type I for main effects
#'
#' @inheritParams multivariate
#'
#' @keywords internal
manova_main_I <- function(x,
                          es_manova,
                          test,
                          digits, CI, lower.tail) {
  R <- x$R.res
  M <- x$M.mu
  F <- x$F.mu
  d <- M - F
  nM <- x$m
  nF <- x$f
  M.sd <- x$M.sdev
  F.sd <- x$F.sdev
  p <- NROW(R)
  r <- NROW(M)
  o.p <- rep(1, p)
  o.r <- rep(1, r)

  Xm <- nM %*% t(o.p) * M
  Xf <- nF %*% t(o.p) * F
  N <- sum(nM) + sum(nF)
  # Must transpose so that R properly treats as a row vector
  Gmean <- t(as.numeric(t(o.r) %*% (Xm + Xf) / N))
  Mc <- M - o.r %*% Gmean
  Fc <- F - o.r %*% Gmean
  J <- matrix(1, nrow = r, ncol = r)
  D <- M - F
  w <- (nM * nF) / (nM + nF)
  weighted.D <- as.numeric(t(D) %*% w)
  T <- diag(sqrt(apply((nM - 1) %o% o.p * M.sd^2 + (nF - 1) %o% o.p * F.sd^2, 2, sum)))

  SSCPb <- t(Mc) %*% diag(nM) %*% Mc + t(Fc) %*% diag(nF) %*% Fc
  SSCPi <- t(D) %*% (w %*% t(o.p) * D) - weighted.D %o% weighted.D / sum(w)
  SSCPe <- T %*% R %*% T
  SSCPsamp <- t(Xm) %*% M + t(Xf) %*% F - t(Xm) %*% J %*% Xm / sum(nM) - t(Xf) %*% J %*% Xf / sum(nF) - SSCPi
  SSCPsex <- SSCPb - SSCPsamp - SSCPi

  vh <- c(1, r - 1)
  ve <- rep(N - 1 - r, 2)

  H <- array(dim = c(p, p, 2))

  H[, , 1] <- SSCPsex
  H[, , 2] <- SSCPsamp

  E <- SSCPe + SSCPi

  vh <- c(1, r - 1)
  ve <- rep(N - 1 - r, 2)

  Stats <- rep(0, 2)
  F.stat <- rep(0, 2)
  df1 <- rep(0, 2)
  df2 <- rep(0, 2)
  prob <- rep(0, 2)
  exact <- rep("", 2)
  eta <- rep(0, 2)
  upper <- rep(0, 2)
  lower <- rep(0, 2)

  for (i in 1:2) {
    if (test == "W") sto <- Wilks(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "R") sto <- Roy(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "P") sto <- Pillai(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "HL") sto <- HL(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    Stats[i] <- sto$Stats
    F.stat[i] <- sto$F.stat
    df1[i] <- sto$df1
    df2[i] <- sto$df2
    prob[i] <- sto$prob
    exact[i] <- sto$exact
    eta [i] <- sto$eta
    lower[i] <- sto$lower.eta
    upper[i] <- sto$upper.eta
  }
  if (test == "W") test.type <- "Wilks"
  if (test == "R") test.type <- "Roy"
  if (test == "P") test.type <- "Pillai"
  if (test == "HL") test.type <- "Hotelling-Lawley"
  term <- c(paste("Sex", exact[1], sep = ""), paste("Pop", exact[2], sep = ""))


  out <- cbind.data.frame(term, vh, round(Stats, digits), round(F.stat, digits), round(df1, 0), round(df2, 3),
    p.value = round(prob, digits), round(eta, digits),
    round(lower, digits), round(upper, digits)
  )

  out <- add_sig(out)


  colnames(out) <- c(
    "term", "df", test.type, "approx.f", "num.df", "den.df", "p.value", "signif",

    "eta", "lower.eta", "upper.eta"
  )



  SSCPsex <- t(d) %*% w %*% t(w) %*% d / as.numeric(t(o.r) %*% w)
  SSCPsamp <- SSCPb - SSCPsex - SSCPi

  H[, , 1] <- SSCPsamp
  H[, , 2] <- SSCPsex

  vh <- c(r - 1, 1)

  for (i in 1:2) {
    if (test == "W") sto <- Wilks(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "R") sto <- Roy(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "P") sto <- Pillai(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "HL") sto <- HL(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    Stats[i] <- sto$Stats
    F.stat[i] <- sto$F.stat
    df1[i] <- sto$df1
    df2[i] <- sto$df2
    prob[i] <- sto$prob
    exact[i] <- sto$exact
    eta [i] <- sto$eta
    lower[i] <- sto$lower.eta
    upper[i] <- sto$upper.eta
  }
  if (test == "W") test.type <- "Wilks"
  if (test == "R") test.type <- "Roy"
  if (test == "P") test.type <- "Pillai"
  if (test == "HL") test.type <- "Hotelling-Lawley"
  term <- c(paste("Pop", exact[1], sep = ""), paste("Sex", exact[2], sep = ""))
  out2 <- cbind.data.frame(term, vh, round(Stats, digits), round(F.stat, digits), round(df1, 0), round(df2, 3),
    p.value = round(prob, digits), round(eta, digits),
    round(lower, digits), round(upper, digits)
  )
  out2 <- add_sig(out2)

  colnames(out2) <- colnames(out)
  if (es_manova == "none") {
    out <- out[1:8]
    out2 <- out2[1:8]
  }
  list(out, out2)
}


# manova_main_II ----------------------------------------------------------

#' MANOVA type II for main effects
#'
#' @inheritParams multivariate
#'
#' @keywords internal
manova_main_II <- function(x,
                           es_manova,
                           test,
                           digits, CI, lower.tail) {
  R <- x$R.res
  M <- x$M.mu
  F <- x$F.mu
  nM <- x$m
  nF <- x$f
  M.sd <- x$M.sdev
  F.sd <- x$F.sdev
  p <- NROW(R)
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
  T <- diag(sqrt(apply((nM - 1) %o% o.p * M.sd^2 + (nF - 1) %o% o.p * F.sd^2, 2, sum)))

  SSCPi <- t(D) %*% (w %*% t(o.p) * D) - weighted.D %o% weighted.D / sum(w)
  SSCPe <- T %*% R %*% T
  SSCPsamp <- t(Xm) %*% M + t(Xf) %*% F - t(Xm) %*% J %*% Xm / sum(nM) - t(Xf) %*% J %*% Xf / sum(nF) - SSCPi
  SSCPsex <- t(D) %*% w %*% t(w) %*% D / as.numeric(t(o.r) %*% w)

  H <- array(dim = c(p, p, 2))

  H[, , 1] <- SSCPsex
  H[, , 2] <- SSCPsamp

  E <- SSCPe + SSCPi

  vh <- c(1, r - 1)
  ve <- rep(N - 1 - r, 2)

  Stats <- rep(0, 2)
  F.stat <- rep(0, 2)
  df1 <- rep(0, 2)
  df2 <- rep(0, 2)
  prob <- rep(0, 2)
  eta <- rep(0, 2)
  lower <- rep(0, 2)
  upper <- rep(0, 2)
  exact <- rep("", 2)

  for (i in 1:2) {
    if (test == "W") sto <- Wilks(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "R") sto <- Roy(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "P") sto <- Pillai(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "HL") sto <- HL(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    Stats[i] <- sto$Stats
    F.stat[i] <- sto$F.stat
    df1[i] <- sto$df1
    df2[i] <- sto$df2
    prob[i] <- sto$prob
    exact[i] <- sto$exact
    eta [i] <- sto$eta
    lower[i] <- sto$lower.eta
    upper[i] <- sto$upper.eta
  }


  if (test == "W") test.type <- "Wilks"
  if (test == "R") test.type <- "Roy"
  if (test == "P") test.type <- "Pillai"
  if (test == "HL") test.type <- "Hotelling-Lawley"
  term <- c(paste("Sex", exact[1], sep = ""), paste("Pop", exact[2], sep = ""))
  out <- cbind.data.frame(term, vh, round(Stats, digits), round(F.stat, digits), round(df1, 0), round(df2, 3),
    p.value = round(prob, digits), round(eta, digits),
    round(lower, digits), round(upper, digits)
  )
  out <- add_sig(out)


  colnames(out) <- c(
    "term", "df", test.type, "approx.f", "num.df", "den.df", "p.value", "signif",

    "eta", "lower.eta", "upper.eta"
  )
  if (es_manova == "none") {
    out <- out[1:8]
  }
  out
}


# manova_I ----------------------------------------------------------------

#' MANOVA type I
#'
#' @inheritParams multivariate
#'
#' @keywords internal
manova_I <- function(x,
                     es_manova,
                     test,
                     digits, CI, lower.tail) {
  R <- x$R.res
  M <- x$M.mu
  F <- x$F.mu
  d <- M - F
  nM <- x$m
  nF <- x$f
  M.sd <- x$M.sdev
  F.sd <- x$F.sdev
  p <- NROW(R)
  r <- NROW(M)
  o.p <- rep(1, p)
  o.r <- rep(1, r)

  Xm <- nM %*% t(o.p) * M
  Xf <- nF %*% t(o.p) * F
  N <- sum(nM) + sum(nF)
  # Must transpose so that R properly treats as a row vector
  Gmean <- t(as.numeric(t(o.r) %*% (Xm + Xf) / N))
  Mc <- M - o.r %*% Gmean
  Fc <- F - o.r %*% Gmean
  J <- matrix(1, nrow = r, ncol = r)
  D <- M - F
  w <- (nM * nF) / (nM + nF)
  weighted.D <- as.numeric(t(D) %*% w)
  T <- diag(sqrt(apply((nM - 1) %o% o.p * M.sd^2 + (nF - 1) %o% o.p * F.sd^2, 2, sum)))

  SSCPb <- t(Mc) %*% diag(nM) %*% Mc + t(Fc) %*% diag(nF) %*% Fc
  SSCPi <- t(D) %*% (w %*% t(o.p) * D) - weighted.D %o% weighted.D / sum(w)
  SSCPe <- T %*% R %*% T
  SSCPsamp <- t(Xm) %*% M + t(Xf) %*% F - t(Xm) %*% J %*% Xm / sum(nM) - t(Xf) %*% J %*% Xf / sum(nF) - SSCPi
  SSCPsex <- SSCPb - SSCPsamp - SSCPi

  vh <- c(1, r - 1, r - 1)
  ve <- rep(N - 2 * r, 3)

  H <- array(dim = c(p, p, 3))

  H[, , 1] <- SSCPsex
  H[, , 2] <- SSCPsamp
  H[, , 3] <- SSCPi

  E <- SSCPe

  vh <- c(1, r - 1, r - 1)
  ve <- rep(N - 2 * r, 3)

  Stats <- rep(0, 3)
  F.stat <- rep(0, 3)
  df1 <- rep(0, 3)
  df2 <- rep(0, 3)
  prob <- rep(0, 3)
  exact <- rep("", 3)
  eta <- rep(0, 3)
  upper <- rep(0, 3)
  lower <- rep(0, 3)

  for (i in 1:3) {
    if (test == "W") sto <- Wilks(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "R") sto <- Roy(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "P") sto <- Pillai(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "HL") sto <- HL(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    Stats[i] <- sto$Stats
    F.stat[i] <- sto$F.stat
    df1[i] <- sto$df1
    df2[i] <- sto$df2
    prob[i] <- sto$prob
    exact[i] <- sto$exact
    eta [i] <- sto$eta
    lower[i] <- sto$lower.eta
    upper[i] <- sto$upper.eta
  }
  if (test == "W") test.type <- "Wilks"
  if (test == "R") test.type <- "Roy"
  if (test == "P") test.type <- "Pillai"
  if (test == "HL") test.type <- "Hotelling-Lawley"

  term <- c(paste("Sex", exact[1], sep = ""), paste("Pop", exact[2], sep = ""), paste("Sex*Pop", exact[3], sep = ""))
  out <- cbind.data.frame(term, vh, round(Stats, digits), round(F.stat, digits), round(df1, 0), round(df2, 3),
    p.value = round(prob, digits), round(eta, digits),
    round(lower, digits), round(upper, digits)
  )
  out <- add_sig(out)


  colnames(out) <- c(
    "term", "df", test.type, "approx.f", "num.df", "den.df", "p.value", "signif",

    "eta", "lower.eta", "upper.eta"
  )


  SSCPsex <- t(d) %*% w %*% t(w) %*% d / as.numeric(t(o.r) %*% w)
  SSCPsamp <- SSCPb - SSCPsex - SSCPi

  vh <- c(r - 1, 1, r - 1)
  ve <- rep(N - 2 * r, 3)

  H[, , 1] <- SSCPsamp
  H[, , 2] <- SSCPsex

  vh <- c(r - 1, 1, r - 1)
  ve <- rep(N - 2 * r, 3)

  for (i in 1:3) {
    if (test == "W") sto <- Wilks(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "R") sto <- Roy(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "P") sto <- Pillai(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "HL") sto <- HL(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    Stats[i] <- sto$Stats
    F.stat[i] <- sto$F.stat
    df1[i] <- sto$df1
    df2[i] <- sto$df2
    prob[i] <- sto$prob
    exact[i] <- sto$exact
    eta [i] <- sto$eta
    lower[i] <- sto$lower.eta
    upper[i] <- sto$upper.eta
  }


  if (test == "W") test.type <- "Wilks"
  if (test == "R") test.type <- "Roy"
  if (test == "P") test.type <- "Pillai"
  if (test == "HL") test.type <- "Hotelling-Lawley"

  term <- c(paste("Pop", exact[1], sep = ""), paste("Sex", exact[2], sep = ""), paste("Sex*Pop", exact[3], sep = ""))
  out2 <- cbind.data.frame(term, vh, round(Stats, digits), round(F.stat, digits), round(df1, 0), round(df2, 3),
    p.value = round(prob, digits), round(eta, digits),
    round(lower, digits), round(upper, digits)
  )

  out2 <- add_sig(out2)


  colnames(out2) <- colnames(out)
  if (es_manova == "none") {
    out <- out[1:8]
    out2 <- out2[1:8]
  }
  list(out, out2)
}


# manova_II ---------------------------------------------------------------
#' MANOVA type II
#'
#' @inheritParams multivariate
#'
#' @keywords internal
manova_II <- function(x,
                      es_manova,
                      test,
                      digits, CI, lower.tail) {
  R <- x$R.res
  M <- x$M.mu
  F <- x$F.mu
  nM <- x$m
  nF <- x$f
  M.sd <- x$M.sdev
  F.sd <- x$F.sdev
  p <- NROW(R)
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
  T <- diag(sqrt(apply((nM - 1) %o% o.p * M.sd^2 + (nF - 1) %o% o.p * F.sd^2, 2, sum)))

  SSCPi <- t(D) %*% (w %*% t(o.p) * D) - weighted.D %o% weighted.D / sum(w)
  SSCPe <- T %*% R %*% T
  SSCPsamp <- t(Xm) %*% M + t(Xf) %*% F - t(Xm) %*% J %*% Xm / sum(nM) - t(Xf) %*% J %*% Xf / sum(nF) - SSCPi
  SSCPsex <- t(D) %*% w %*% t(w) %*% D / as.numeric(t(o.r) %*% w)

  H <- array(dim = c(p, p, 3))

  H[, , 1] <- SSCPsex
  H[, , 2] <- SSCPsamp
  H[, , 3] <- SSCPi

  E <- SSCPe

  vh <- c(1, r - 1, r - 1)
  ve <- rep(N - 2 * r, 3)

  Stats <- rep(0, 3)
  F.stat <- rep(0, 3)
  df1 <- rep(0, 3)
  df2 <- rep(0, 3)
  prob <- rep(0, 3)
  exact <- rep("", 3)
  eta <- rep(0, 3)
  lower <- rep(0, 3)
  upper <- rep(0, 3)

  for (i in 1:3) {
    if (test == "W") sto <- Wilks(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "R") sto <- Roy(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "P") sto <- Pillai(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "HL") sto <- HL(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    Stats[i] <- sto$Stats
    F.stat[i] <- sto$F.stat
    df1[i] <- sto$df1
    df2[i] <- sto$df2
    prob[i] <- sto$prob
    exact[i] <- sto$exact
    eta [i] <- sto$eta
    lower[i] <- sto$lower.eta
    upper[i] <- sto$upper.eta
  }


  if (test == "W") test.type <- "Wilks"
  if (test == "R") test.type <- "Roy"
  if (test == "P") test.type <- "Pillai"
  if (test == "HL") test.type <- "Hotelling-Lawley"
  term <- c(paste("Sex", exact[1], sep = ""), paste("Pop", exact[2], sep = ""), paste("Sex*Pop", exact[3], sep = ""))
  out <- cbind.data.frame(term, vh, round(Stats, digits), round(F.stat, digits), round(df1, 0), round(df2, 3),
    p.value = round(prob, digits), round(eta, digits),
    round(lower, digits), round(upper, digits)
  )

  out <- add_sig(out)


  colnames(out) <- c(
    "term", "df", test.type, "approx.f", "num.df", "den.df", "p.value", "signif",

    "eta", "lower.eta", "upper.eta"
  )
  if (es_manova == "none") {
    out <- out[1:8]
  }
  out
}


# manova_III --------------------------------------------------------------
#' MANOVA type III
#'
#' @inheritParams multivariate
#'
#' @keywords internal
manova_III <- function(x,
                       es_manova,
                       test,
                       digits, CI, lower.tail) {
  R <- x$R.res
  M <- x$M.mu
  F <- x$F.mu
  nM <- x$m
  nF <- x$f
  N.M <- sum(nM)
  N.F <- sum(nF)
  N <- N.M + N.F
  M.sd <- x$M.sdev
  F.sd <- x$F.sdev
  p <- NROW(R)
  r <- NROW(M)
  o.p <- rep(1, p)
  o.r <- rep(1, r)

  T <- diag(sqrt(apply((nM - 1) %o% o.p * M.sd^2 + (nF - 1) %o% o.p * F.sd^2, 2, sum)))
  SSCPe <- T %*% R %*% T

  Pop.Block <- t(stats::contr.sum(r)) %*% diag(nM + nF) %*% stats::contr.sum(r)

  contrast <- stats::contr.sum(2) %x% stats::contr.sum(r)
  first <- t(contrast) %*% diag(c(nF, nM))
  first <- first[, 1:r] + first[, (r + 1):(2 * r)]
  second <- stats::contr.sum(r)
  Interact.Block <- first %*% second

  Xm <- nM %*% t(o.p) * M
  Xf <- nF %*% t(o.p) * F
  Male.sum <- apply(Xm, 2, sum)
  Female.sum <- apply(Xf, 2, sum)
  Pop.sums <- Xm + Xf
  Pops <- t(stats::contr.sum(r)) %*% Pop.sums
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
  SSCPsex <- t(AB) %*% solve(inv.contr) %*% AB

  A <- cbind(matrix(0, nrow = r - 1, ncol = 2), diag(r - 1), matrix(0,
    nrow =
      r - 1, ncol = r - 1
  ))
  AB <- A %*% B
  inv.contr <- A %*% inv.XX %*% t(A)
  SSCPsamp <- t(AB) %*% solve(inv.contr) %*% AB


  A <- cbind(
    matrix(0, nrow = r - 1, ncol = 2), matrix(0, nrow = r - 1, ncol = r - 1),
    diag(r - 1)
  )
  AB <- A %*% B
  inv.contr <- A %*% inv.XX %*% t(A)
  SSCPi <- t(AB) %*% solve(inv.contr) %*% AB

  H <- array(dim = c(p, p, 3))

  H[, , 1] <- SSCPsex
  H[, , 2] <- SSCPsamp
  H[, , 3] <- SSCPi

  E <- SSCPe

  vh <- c(1, r - 1, r - 1)
  ve <- rep(N - 2 * r, 3)

  Stats <- rep(0, 3)
  F.stat <- rep(0, 3)
  df1 <- rep(0, 3)
  df2 <- rep(0, 3)
  prob <- rep(0, 3)
  exact <- rep("", 3)
  eta <- rep(0, 3)
  lower <- rep(0, 3)
  upper <- rep(0, 3)

  for (i in 1:3) {
    if (test == "W") sto <- Wilks(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "R") sto <- Roy(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "P") sto <- Pillai(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    if (test == "HL") sto <- HL(p, H[, , i], E, vh[i], ve[i], CI, lower.tail)
    Stats[i] <- sto$Stats
    F.stat[i] <- sto$F.stat
    df1[i] <- sto$df1
    df2[i] <- sto$df2
    prob[i] <- sto$prob
    exact[i] <- sto$exact
    eta [i] <- sto$eta
    lower[i] <- sto$lower.eta
    upper[i] <- sto$upper.eta
  }


  if (test == "W") test.type <- "Wilks"
  if (test == "R") test.type <- "Roy"
  if (test == "P") test.type <- "Pillai"
  if (test == "HL") test.type <- "Hotelling-Lawley"
  term <- c(paste("Sex", exact[1], sep = ""), paste("Pop", exact[2], sep = ""), paste("Sex*Pop", exact[3], sep = ""))

  out <- cbind.data.frame(term, vh, round(Stats, digits), round(F.stat, digits), round(df1, 0), round(df2, 3),
    p.value = round(prob, digits), round(eta, digits),
    round(lower, digits), round(upper, digits)
  )

  out <- add_sig(out)
  colnames(out) <- c(
    "term", "df", test.type, "approx.f", "num.df", "den.df", "p.value", "signif",

    "eta", "lower.eta", "upper.eta"
  )
  if (es_manova == "none") {
    out <- out[1:8]
  }
  out
}
