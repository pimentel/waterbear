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

# TODO: implement rank base proposal

# mass function for adding a variable from bayesian sparse linear mixed models
# mixture * Categorical(0, max) + (1 - mixture) * TruncatedGeometric(p_geometric)
rqt = nimbleFunction(
  run = function(n = integer(0), max = integer(0), mixture = double(0), p_geometric = double(0)) {
    returnType(double(1))
    p = rep(p_geometric, max)
    q = 1 - p[1]
    s = p[1]
    for (i in 2:max) {
      p[i] = p[i - 1] * q
      s = p[i] + s
    }
    p = p / s
    for (i in 1:max) {
      p[i] = mixture / max + (1 - mixture) * p[i]
    }
    return(rcat(n, p) - 1)
  })

get_nth_zero = nimbleFunction(
  run = function(indicator = integer(1), order = integer(1), n = integer(0)) {
    returnType(integer(i))
    nz = 0
    for (i in 1:length(indicator)) {
      if (indicator[order[i]] == 0) {
        nz = nz + 1
        if (nz == n) {
          return(order[i])
        }
      }
    }
    # n is larger than the number of zeros
    return(-1)
  })

rank_RW <- nimbleFunction(

    contains = sampler_BASE,

    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        r <- control$rank
        o <- control$order
    },

    run = function() {
        # initial model logProb
        model_lp_initial <- getLogProb(model, calcNodes)
        # choose to add, remove, or swap
        which_move <- rcat(1, c(0.45, 0.45, 0.1))
        current_values = model[[target]]
        proposal = current_values
        if (which_move == 1) {
          # TODO: run rqt and get a variable to add
          n_zeros = length(which(current_values == 0))
          which_rank = rqt(1, n_zeros, 0.3, 2 / n_zeros) + 1
          i = get_nth_zero(current_values, o, which_rank)
          proposal[i] = 1
        } else if (which_move == 2) {
          # TODO: uniformly pick one to remove
          which_ones = which(current_value == 1)
          i = rcat(1, rep(1, length(which_ones)))
          proposal[i] = 0
        } else {
          # TODO: uniformly pick one to remove and one to add
          which_zeros = which(current_value == 0)
          i = rcat(1, rep(1, length(which_zeros)))
          proposal[i] = 1
          which_ones = which(current_value == 1)
          j = rcat(1, rep(1, length(which_ones)))
          proposal[j] = 0
        }
        # store proposal into model
        model[[target]] <<- proposal
        # proposal model logProb
        model_lp_proposed <- model$calculate(calcNodes)

        # log-Metropolis-Hastings ratio
        log_MH_ratio <- model_lp_proposed - model_lp_initial

        # Metropolis-Hastings step: determine whether or
        # not to accept the newly proposed value
        u <- runif(1, 0, 1)
        if(u < exp(log_MH_ratio)) jump <- TRUE
        else                      jump <- FALSE

        # keep the model and mvSaved objects consistent
        if(jump) copy(from = model, to = mvSaved, row = 1,
                         nodes = calcNodes, logProb = TRUE)
        else     copy(from = mvSaved, to = model, row = 1,
                         nodes = calcNodes, logProb = TRUE)
    },

    methods = list(   reset = function () {}   )
)
