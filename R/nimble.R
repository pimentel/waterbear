#' dirichlet multinomial likelihood
#'
#' used for MCMC. slightly modified version from the NIMBLE manual.
#'
#' @param x a vector of data
#' @param alpha a vector of concentration parameters
#' @param size the sum of x.
#' @param log if TRUE return the log probability. otherwise, just return probability.
ddirchmulti <- nimbleFunction(
  run = function(x = double(1), alpha = double(1), size = double(0),
    log = integer(0, default = 0)) {
    returnType(double(0))
    logProb <- lgamma(size + 1) - sum(lgamma(x + 1)) + lgamma(sum(alpha)) -
      sum(lgamma(alpha)) + sum(lgamma(alpha + x)) - lgamma(sum(alpha) +
        size)
    if(log) return(logProb)
    else return(exp(logProb))
  })

rdirchmulti <- nimbleFunction(
  run = function(n = integer(0), alpha = double(1), size = double(0)) {
    returnType(double(1))
    if(n != 1) print("rdirchmulti only allows n = 1; using n = 1.")
    p <- rdirch(1, alpha)
    return(rmulti(1, size = size, prob = p))
  })

#' convert a N-dimensional dirichlet to a N-1-dimensional set of normal cutoffs
#'
#' match the quantiles from the Dirichlet to a standard normal and return the corresponding
#' cutoffs to match the bins
dirichlet_to_normal_bins = nimbleFunction(
  run = function(alpha = double(1)) {
    returnType(double(1))
    a0 = sum(alpha)
    p = alpha / a0
    qp = numeric(length(p))
    qp[1] = p[1]
    for (i in 2:length(p)) {
      qp[i] = qp[i - 1] + p[i]
    }
    cutoff = numeric(length(qp) - 1)
    for (i in 1:(length(qp) - 1)) {
      cutoff[i] = probit(qp[i])
    }
    return(cutoff)
  })

#' convert standard normal cutoffs to get the corresponding mass in each bin
#'
#' given a set of cutoffs, convert these bins to a probability mass function.
cutoff_to_probability = nimbleFunction(
  run = function(cutoff = double(1), offset = double(0)) {
    returnType(double(1))
    p = numeric(length(cutoff) + 1)
    p[1] = phi(cutoff[1] - offset)
    total = p[1]
    for (i in 2:length(cutoff)) {
      p[i] = phi(cutoff[i] - offset) - phi(cutoff[i - 1] - offset)
      total = total + p[i]
    }
    p[length(cutoff) + 1] = 1 - total
    return(p)
  })

sample_specific_dispersion_model = nimbleCode({
  dispersion ~ dexp(1.0 / dispersion_prior_mean)
  for (n in 1:N) {
    bin_alpha[n, 1:N_bins] ~ ddirch(bin_alpha_prior[n, 1:N_bins])
    cutoffs[n, 1:N_cutoffs] <- dirichlet_to_normal_bins(bin_alpha[n, 1:N_bins])
    sample_dispersion[n] ~ dexp(1.0 / dispersion)
  }
  psi ~ dbeta(10, 10)
  sigma_gene ~ dgamma(1, 0.10)
  sigma_guide ~ dgamma(1, 0.10)
  for (gene in 1:N_genes) {
    gene_inclusion[gene] ~ dbern(psi)
    gene_shift[gene] ~ dnorm(0, sd = sigma_gene)
  }
  for (g in 1:N_nt) {
    for (n in 1:N) {
      x[n, nt_data_index[g], 1:N_bins] ~ ddirchmulti(
        sample_dispersion[n] * bin_alpha[n, 1:N_bins],
        x_total[n, nt_data_index[g]])
    }
  }
  for (g in 1:N_guides) {
    guide_shift[g] ~ dnorm(gene_shift[guide_to_gene[g]], sd = sigma_guide)
    total_shift[g] <- gene_inclusion[guide_to_gene[g]] * guide_shift[g]
    for (n in 1:N) {
      q[n, g, 1:N_bins] <- cutoff_to_probability(cutoffs[n, 1:N_cutoffs], total_shift[g])
      x[n, guide_data_index[g], 1:N_bins] ~ ddirchmulti(
        sample_dispersion[n] * q[n, g, 1:N_bins],
        x_total[n, guide_data_index[g]])
    }
  }
})

shared_dispersion_model = nimbleCode({
  dispersion ~ dexp(1.0 / dispersion_prior_mean)
  for (n in 1:N) {
    bin_alpha[n, 1:N_bins] ~ ddirch(bin_alpha_prior[n, 1:N_bins])
    cutoffs[n, 1:N_cutoffs] <- dirichlet_to_normal_bins(bin_alpha[n, 1:N_bins])
  }
  psi ~ dbeta(10, 10)
  sigma_gene ~ dgamma(1, 0.10)
  sigma_guide ~ dgamma(1, 0.10)
  for (gene in 1:N_genes) {
    gene_inclusion[gene] ~ dbern(psi)
    gene_shift[gene] ~ dnorm(0, sd = sigma_gene)
  }
  for (g in 1:N_nt) {
    for (n in 1:N) {
      x[n, nt_data_index[g], 1:N_bins] ~ ddirchmulti(
        dispersion * bin_alpha[n, 1:N_bins],
        x_total[n, nt_data_index[g]])
    }
  }
  for (g in 1:N_guides) {
    guide_shift[g] ~ dnorm(gene_shift[guide_to_gene[g]], sd = sigma_guide)
    total_shift[g] <- gene_inclusion[guide_to_gene[g]] * guide_shift[g]
    for (n in 1:N) {
      q[n, g, 1:N_bins] <- cutoff_to_probability(cutoffs[n, 1:N_cutoffs], total_shift[g])
      x[n, guide_data_index[g], 1:N_bins] ~ ddirchmulti(
        dispersion * q[n, g, 1:N_bins],
        x_total[n, guide_data_index[g]])
    }
  }
})
