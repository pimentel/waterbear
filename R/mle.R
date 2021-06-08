dirichlet_multinomial_null_ll = function(counts, bin_sizes, phi, aggregate = TRUE) {
  # make sure that bin_sizes is a probability distribution
  if (all(apply(bin_sizes, 1, sum) != rep(1, nrow(bin_sizes)))) {
    print(bin_sizes)
    stop('bin_sizes != 1')
  }
  # stopifnot(all(apply(bin_sizes, 1, sum) == rep(1, nrow(bin_sizes))))
  a = phi * bin_sizes
  total_counts = apply(counts, c(1, 2), sum)

  # for each replicate, iterate through and compute the log probability
  ll = sapply(1:dim(counts)[1],
    function(n) {
      sum(sapply(1:dim(counts)[2],
          function(g) {
            ddirchmulti(counts[n, g, ], a[n, ], total_counts[n, g], TRUE)
          }))
    })
  if (aggregate) {
    ll = sum(ll)
  }
  ll
}


estimate_dispersion_from_controls = function(counts, bin_sizes) {
  fit = optim(c(1),
    function(x) {
      -dirichlet_multinomial_null_ll(counts, bin_sizes, x)
    }, method = 'Brent', lower = 0.0005, upper = 1000)
  fit$par
}

estimate_bin_sizes = function(counts, prior_bin_sizes, dispersion) {
  N = dim(counts)[1]
  prior_bin_cutoffs = t(apply(prior_bin_sizes, 1, dirichlet_to_normal_bins))
  for (n in 1:N) {
    current_fit = optim(
      prior_bin_cutoffs[n, ],
      function(x) {
        y = prior_bin_cutoffs
        y[n, ] = x
        z = t(apply(y, 1, cutoff_to_probability, offset = 0))
        -dirichlet_multinomial_null_ll(counts, z, dispersion)
      })
    prior_bin_cutoffs[n, ] = current_fit$par
  }
  t(apply(prior_bin_cutoffs, 1, cutoff_to_probability, offset = 0))
  # prior_bin_sizes
}


# assumes counts is a matrix with samples on the rows and columns representing bins
# you will get this if you do all_counts[, g, ]
ddmult_per_guide_ll = function(counts, cutoffs, mu, dispersion) {
  ll = sapply(1:dim(counts)[1],
    function(n) {
      cur_guide_counts = counts[n, ]
      a = cutoff_to_probability(cutoffs[n, ], mu)
      a = a * dispersion
      total = sum(cur_guide_counts)
      ddirchmulti(cur_guide_counts, a, total, TRUE)
    })
  sum(ll)
}

estimate_guide_mu = function(counts, bin_sizes, dispersion) {
  cutoffs = t(apply(bin_sizes, 1, dirichlet_to_normal_bins))
  # go guide by guide and estimate each mu
  sapply(1:dim(counts)[2],
    function(g) {
      fit = optim(c(0),
      function(x) {
        -ddmult_per_guide_ll(counts[, g, ], cutoffs, x, dispersion)
      }, method = 'Brent', lower = -10, upper = 10)
      fit$par
    })
}

estimate_gene_mu = function(guide_mu, guide_to_gene) {
  mu_by_gene = split(guide_mu, guide_to_gene)
  mu_gene = sapply(mu_by_gene, mean)
  mu_gene = mu_gene[as.character(sort(unique(guide_to_gene)))]

  sd_mu_gene = sapply(mu_by_gene, sd)
  w = mu_gene / sd_mu_gene
  data.frame(
    gene_id = as.integer(names(mu_gene)),
    mu = mu_gene,
    sd_mu = sd_mu_gene,
    wald = w,
  stringsAsFactors = FALSE)
}

lrt = function(counts, bin_sizes, phi, mu, guide_to_gene) {
  cutoffs = t(apply(bin_sizes, 1, dirichlet_to_normal_bins))
  ll_0 = sapply(1:dim(counts)[2],
    function(g) {
      ddmult_per_guide_ll(counts[, g, ], cutoffs, 0, phi)
    })

  ll_1 = sapply(1:dim(counts)[2],
    function(g) {
      ddmult_per_guide_ll(counts[, g, ], cutoffs, mu[g], phi)
    })

  ll_index = split(1:length(ll_0), guide_to_gene)
  gene_lrt = lapply(ll_index,
    function(i) {
      gene_0 = sum(ll_0[i])
      gene_1 = sum(ll_1[i])
      data.frame(ll_0 = gene_0, ll_1 = gene_1, lambda = -2 * (gene_0 - gene_1))
    })
  gene_lrt = dplyr::bind_rows(gene_lrt)
  gene_lrt = dplyr::mutate(gene_lrt, mapping = names(ll_index))
  gene_lrt
}
