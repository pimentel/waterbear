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
    # if (length(cutoff) == 1) {
    #   p[2] = 1.0 - p[1]
    #   return(p)
    # }
    total = p[1]
    if (length(cutoff) > 1) {
      for (i in 2:length(cutoff)) {
        p[i] = phi(cutoff[i] - offset) - phi(cutoff[i - 1] - offset)
        total = total + p[i]
      }
    }
    p[length(cutoff) + 1] = 1 - total
    return(p)
  })

sample_specific_dispersion_model_no_controls = nimbleCode({
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
  if (N_nt > 0) {
    for (g in 1:N_nt) {
      for (n in 1:N) {
        x[n, nt_data_index[g], 1:N_bins] ~ ddirchmulti(
          sample_dispersion[n] * bin_alpha[n, 1:N_bins],
          x_total[n, nt_data_index[g]])
      }
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
    returnType(integer(0))
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
    returnType(integer(0))
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
    print('get_nth_zero: didnt find enough zeros? error')
    return(-1)
  })

get_nth_value = nimbleFunction(
  run = function(indicator = integer(0), n = integer(0), value = integer(0)) {
    returnType(integer(0))
    nvalue = 0
    for (i in 1:length(indicator)) {
      if (indicator[i] == value) {
        nvalue = nvalue + 1
        if (nvalue == n) {
          return(i)
        }
      }
    }
    print('get_nth_value: didnt find enough zeros? error')
    return(-1)
        })


rank_RW <- nimbleFunction(
  name = 'rank_RW',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    ts_order <- control$order
    times_ran = 0
    times_1 = 0
    times_2 = 0
    times_3 = 0
    times_accepted = 0
  },
  run = function() {
    # initial model logProb
    # print('start')
    times_ran <<- times_ran + 1
    model_lp_initial <- getLogProb(model, calcNodes)
    # choose to add, remove, or swap
    which_move <- rcat(1, c(0.45, 0.45, 0.1))
    current_values = integer(length(model[[target]]), value = model[[target]])
    total_sum = sum(current_values)
    if (total_sum == length(current_values)) {
      # all ones
      which_move = 2
    } else if (total_sum == 0) {
      # all zeroes
      which_move = 1
    }
    proposal = current_values
    if (which_move == 1) {
      # print('move 1')
      # run rqt and get a variable to add
      n_zeros = length(which(current_values == 0))
      geometric_mean = 2.0 / n_zeros
      which_rank = rqt(1, n_zeros, 0.3, geometric_mean)
      which_rank = which_rank + 1
      # print('attempting to get the ')
      # print(which_rank)
      # print('zero')
      i = get_nth_zero(current_values, ts_order, which_rank)
      # print('about to assign adding')
      proposal[i] = 1
      # print('')
      # print(current_values)
      # print(proposal)
      times_1 <<- times_1 + 1
      # print('done adding')
    } else if (which_move == 2) {
      # print('move 2')
      # uniformly pick one to remove
      # print('which')
      which_ones = which(current_values == 1)
      # print('rcat')
      # print(length(which_ones))
      i = rcat(1, rep(1, length(which_ones)))
      # print('about to remove one')
      proposal[which_ones[i]] = 0
      # print('done removing')
      times_2 <<- times_2 + 1
    } else {
      # print('move 3')
      # uniformly pick one to remove and one to add
      which_zeros = which(current_values == 0)
      i = rcat(1, rep(1, length(which_zeros)))
      which_ones = which(current_values == 1)
      j = rcat(1, rep(1, length(which_ones)))
      # print('about to do a swap')
      proposal[which_zeros[i]] = 1
      proposal[which_ones[j]] = 0
      times_3 <<- times_3 + 1
      # print('done swapping')
    }
    # store proposal into model
    model[[target]] <<- proposal
    # proposal model logProb
    model_lp_proposed <- model$calculate(calcNodes)

        # log-Metropolis-Hastings ratio
        log_MH_ratio <- model_lp_proposed - model_lp_initial
    # print('')
    # print(model_lp_initial)
    # print(model_lp_proposed)
    # print(exp(log_MH_ratio))


        # Metropolis-Hastings step: determine whether or
        # not to accept the newly proposed value
        u <- runif(1, 0, 1)
        if(u < exp(log_MH_ratio)) {
          jump = TRUE
          times_accepted <<- times_accepted + 1
          # print('accepted!')
        }
        else{
          jump = FALSE
        }

        # keep the model and mvSaved objects consistent
        if(jump) copy(from = model, to = mvSaved, row = 1,
                         nodes = calcNodes, logProb = TRUE)
        else     copy(from = mvSaved, to = model, row = 1,
                         nodes = calcNodes, logProb = TRUE)
        # print('stop')
    },

    methods = list(
      reset = function () {
        times_1 <<- 0
        times_2 <<- 0
        times_3 <<- 0
        times_accepted <<- 0
        times_ran <<- 0
      },
      getAcceptanceHistory = function() {
        returnType(double(1))
        # nimbleList(
        #   variable = character(6, values =
        #     c('times_ran',
        #     'times_accepted',
        #     'times_1',
        #     'times_2',
        #     'times_3',
        #     'acceptance_rate')
        #     ),
          # result = double(6, values = c(
        acceptance_rate = times_accepted / times_ran
        return(c(
            times_ran,
            times_accepted,
            times_1,
            times_2,
            times_3,
            acceptance_rate
          ))

      })
)
