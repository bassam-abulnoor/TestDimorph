
# Wilks lambda ------------------------------------------------------------
#' Wilks
#'
#' @param p Number of dependent variables (measurements)
#' @param H Hypothesis SSCP matrix
#' @param E Error SSCP matrix
#' @param v.h Hypothesis degrees of freedom
#' @param v.e Error degrees of freedom
#'
#' @keywords internal
Wilks <- function(p, H, E, v.h, v.e, CI, lower.tail) {
  r <- min(p, v.h)
  Lamb <- det(E) / det(E + H)
  t <- v.e + v.h
  q <- v.h
  if (q == 1) {
    df1 <- p
    df2 <- t - p
    R <- (1 - Lamb) / Lamb * (df2 / df1)
    prob <- pf(R, df1, df2, lower.tail = lower.tail)
    eta <- 1 - Lamb^(1 / r)
    eff <- eff_CI(
      f = R,
      CI = CI,
      eff = eta,
      df1 = df1,
      df2 = df2
    ) %>% as.numeric()
    lower <- eff[2]
    upper <- eff[3]
    return(list(
      Stats = Lamb, F.stat = R, df1 = df1, df2 = df2, prob = prob, eta = eta,
      lower.eta = lower, upper.eta = upper, exact = "(E)"
    ))
  }
  if (q == 2) {
    sqrt.L <- sqrt(Lamb)
    R <- (1 - sqrt.L) / sqrt.L * ((v.e - p + 1) / p)
    df1 <- 2 * p
    df2 <- 2 * (t - p - 1)
    prob <- pf(R, df1, df2, lower.tail = lower.tail)
    eta <- 1 - Lamb^(1 / r)
    eff <- eff_CI(
      f = R,
      CI = CI,
      eff = eta,
      df1 = df1,
      df2 = df2
    ) %>% as.numeric()
    lower <- eff[2]
    upper <- eff[3]
    return(list(
      Stats = Lamb, F.stat = R, df1 = df1, df2 = df2, prob = prob, eta = eta,
      lower.eta = lower, upper.eta = upper, exact = "(E)"
    ))
  }

  if (p == 2) {
    sqrt.L <- sqrt(Lamb)
    R <- (1 - sqrt.L) / sqrt.L * ((t - q - 1) / q)
    df1 <- 2 * q
    df2 <- 2 * (t - q - 1)
    prob <- pf(R, df1, df2, lower.tail = lower.tail)
    eta <- 1 - Lamb^(1 / r)
    eff <- eff_CI(
      f = R,
      CI = CI,
      eff = eta,
      df1 = df1,
      df2 = df2
    ) %>% as.numeric()
    lower <- eff[2]
    upper <- eff[3]
    return(list(
      Stats = Lamb, F.stat = R, df1 = df1, df2 = df2, prob = prob, eta = eta,
      lower.eta = lower, upper.eta = upper, exact = "(E)"
    ))
  }
  s <- sqrt((p^2 * q^2 - 4) / (p^2 + q^2 - 5))
  m <- t - (p + q + 1) / 2
  L <- (p * q - 2) / 4
  df1 <- p * q
  df2 <- m * s - 2 * L

  Lamb.s <- Lamb^(1 / s)
  R <- (1 - Lamb.s) / Lamb.s * (df2 / df1)
  prob <- pf(R, df1, df2, lower.tail = lower.tail)
  eta <- 1 - Lamb^(1 / r)
  eff <- eff_CI(
    f = R,
    CI = CI,
    eff = eta,
    df1 = df1,
    df2 = df2
  ) %>% as.numeric()
  lower <- eff[2]
  upper <- eff[3]
  return(list(
    Stats = Lamb, F.stat = R, df1 = df1, df2 = df2, prob = prob, eta = eta,
    lower.eta = lower, upper.eta = upper, exact = ""
  ))
}


# Pillai ------------------------------------------------------------------

#' Pillai
#'
#' @inheritParams Wilks
#' @keywords internal
Pillai <- function(p, H, E, v.h, v.e, CI, lower.tail) {
  s <- min(p, v.h)
  H.inv.E <- H %*% solve(E)
  c <- Re(eigen(H.inv.E)$val)
  Vs <- 0
  for (i in 1:s) {
    Vs <- Vs + c[i] / (1 + c[i])
  }
  m <- (abs(p - v.h) - 1) / 2
  n <- (v.e - p - 1) / 2
  F <- (2 * n + s + 1) / (2 * m + s + 1) * Vs / (s - Vs)

  df1 <- s * (2 * m + s + 1)
  df2 <- s * (2 * n + s + 1)
  prob <- pf(F, df1, df2, lower.tail = lower.tail)
  eta <- Vs / s
  eff <- eff_CI(
    f = F,
    CI = CI,
    eff = eta,
    df1 = df1,
    df2 = df2
  ) %>% as.numeric()
  lower <- eff[2]
  upper <- eff[3]
  if (v.h == 1) {
    return(list(
      Stats = Vs, F.stat = F, df1 = df1, df2 = df2, prob = prob, eta = eta,
      lower.eta = lower, upper.eta = upper, exact = "(E)"
    ))
  }
  return(list(
    Stats = Vs, F.stat = F, df1 = df1, df2 = df2, prob = prob, eta = eta,
    lower.eta = lower, upper.eta = upper, exact = ""
  ))
}

# Hotelling-Lawley --------------------------------------------------------


#' Hotelling-Lawley
#'
#' @inheritParams Wilks
#' @keywords internal
HL <- function(p, H, E, v.h, v.e, CI, lower.tail) {
  df.h <- v.h
  df.e <- v.e
  r <- min(c(p, df.h))
  b <- max(c(p, df.h))
  m <- df.e - (p - df.h + 1) / 2

  H.inv.E <- H %*% solve(E)
  e <- Re(eigen(H.inv.E)$val[1:r])
  T <- sum(e)

  df1 <- b * r
  df2 <- r * (df.e - p - 1) + 2
  F <- (r * (df.e - p - 1) + 2) / (r^2 * b) * T
  prob <- pf(F, df1, df2, lower.tail = lower.tail)
  eta <- (T / r) / (T / r + 1)
  eff <- eff_CI(
    f = F,
    CI = CI,
    eff = eta,
    df1 = df1,
    df2 = df2
  ) %>% as.numeric()
  lower <- eff[2]
  upper <- eff[3]
  if (v.h == 1) {
    return(list(
      Stats = T, F.stat = F, df1 = df1, df2 = df2, prob = prob, eta = eta,
      lower.eta = lower, upper.eta = eta, exact = "(E)"
    ))
  }
  return(list(
    Stats = T, F.stat = F, df1 = df1, df2 = df2, prob = prob, eta = eta,
    lower.eta = eta, upper.eta = upper, exact = ""
  ))
}


# Roy's largest root ------------------------------------------------------
#' Roy's largest root
#'
#' @inheritParams Wilks
#' @keywords internal

Roy <- function(p, H, E, v.h, v.e, CI, lower.tail) {
  H.inv.E <- H %*% solve(E)
  R <- Re(eigen(H.inv.E)$val[1])
  v1 <- max(p, v.h)
  v2 <- v.e - v1 + v.h
  df1 <- v1
  df2 <- v2
  F <- df2 / df1 * R
  prob <- pf(F, df1, df2, lower.tail = lower.tail)
  eta <- R / (R + 1)
  eff <- eff_CI(
    f = F,
    CI = CI,
    eff = eta,
    df1 = df1,
    df2 = df2
  ) %>% as.numeric()
  lower <- eff[2]
  upper <- eff[3]
  if (v.h == 1) {
    return(list(
      Stats = R, F.stat = F, df1 = df1, df2 = df2, prob = prob, eta = eta,
      lower.eta = lower, upper.eta = upper, exact = "(E)"
    ))
  }
  return(list(
    Stats = R, F.stat = F, df1 = df1, df2 = df2, prob = prob, eta = eta,
    lower.eta = lower, upper.eta = upper, exact = ""
  ))
}
