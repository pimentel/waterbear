#' convert a raw counts table to a three-dimensional array
#'
#' takes a two-dimensional data frame with guides on the rows and samples + bins on the
#' columns and converts it to a three-dimensional array [sample, guide, bins].
#'
#' @param raw_counts the data frame of the raw counts
#' @param sample_mapping a data frame with the following column names:
#'  - c_name: the name of the column as it is represented in raw_counts
#'  - sample: a unique identifier for the sample (replicate)
#'  - bin: the bin that this column refers to. must be consistent with `ordering`.
#' @param ordering an ordering corresponding to the bin column in `sample_mapping`.
#'  e.g.: c('low', 'low_mid', 'mid_high', 'high').
#' @return a three-dimensional array corresponding to [sample, guide, bins]
#' @export
wb_counts_to_array = function(raw_counts, sample_mapping, ordering) {
  sample_names = unique(sample_mapping$sample)
  sample_mapping = dplyr::mutate(sample_mapping, bin = factor(bin, levels = ordering))

  counts_by_sample = lapply(sample_names,
    function(s) {
      df = dplyr::filter(sample_mapping, sample == s)
      df = dplyr::arrange(df, bin)
      mat = as.matrix(raw_counts[, df$c_name])
      rownames(mat) = raw_counts$sgRNA
      colnames(mat) = df$bin
      mat
    })

  counts_array = array(NA,
    dim = c(
      length(counts_by_sample),
      nrow(counts_by_sample[[1]]),
      ncol(counts_by_sample[[1]])
      ))

  for (i in 1:length(counts_by_sample)) {
    counts_array[i, , ] = counts_by_sample[[i]]
  }

  dimnames(counts_array) = list(
    sample_names,
    rownames(counts_by_sample[[1]]),
    colnames(counts_by_sample[[1]]))

  counts_array
}

#' counts to waterbear object
#'
#' this function takes a count table and converts it to a water bear object
#'
#' @param counts_array an array organized in the following dimension:
#' @param gene_mapping a data frame mapping of guide names to gene names. requires the following column names: (1) guide, (2) gene.
#' @param control_guide_regex a regular expression used to find/match control guides. default is 'Non-'.
#' @param bin_size_prior the expected mass in each bin. If NULL, defaults to uniform (e.g. c(0.25, 0.25, 0.25, 0.25)).
#' @return a water bear object that inference can be performed on.
#' @export
wb_make_object = function(
  counts_array,
  gene_mapping,
  control_guide_regex = 'Non-',
  bin_size_prior = NULL
  ) {
  guide_names = data.frame(guide = dimnames(counts_array)[[2]])
  guide_names = dplyr::mutate(guide_names, i = 1:nrow(guide_names))
  N_bins = dim(counts_array)[3]

  if (is.null(bin_size_prior)) {
    bin_size_prior = rep(1 / N_bins, length.out = N_bins)
  } else {
    if (sum(bin_size_prior) != 1.0) {
      stop('bin_size_prior should sum to 1.0')
    }
  }

  if (is.na(control_guide_regex) || is.null(control_guide_regex) || control_guide_regex == '') {
    nt_guide_names = data.frame()
    control_guide_regex = 'ANYWHERE_BUT_HERE'
  } else {
    nt_guide_names = inner_join(guide_names,
      dplyr::filter(gene_mapping, grepl(control_guide_regex, gene)), by = 'guide')
    if (nrow(nt_guide_names) < 1) {
      stop('couldn\'t find any control guides. check the regular expression in `control_guide_regex` and make sure it finds genes in `gene_mapping`.')
  }
  }
  test_guide_names = dplyr::inner_join(guide_names,
    dplyr::filter(gene_mapping, !grepl(control_guide_regex, gene)), by = c('guide'))
  test_guide_names = dplyr::mutate(test_guide_names, mapping = as.integer(factor(gene)))
  test_guide_names = dplyr::arrange(test_guide_names, i)

  gg_data = list(x = counts_array)
  dispersion_init = 300
  gg_const = list(
    N = dim(gg_data$x)[1],
    N_guides = nrow(test_guide_names),
    x_total = apply(gg_data$x, c(1, 2), sum),
    N_genes = length(unique(test_guide_names$mapping)),
    N_bins = N_bins,
    N_nt = nrow(nt_guide_names),
    nt_index = (if(nrow(nt_guide_names) > 0) {1:nrow(nt_guide_names)} else {NULL}),
    nt_data_index = nt_guide_names$i,
    guide_index = 1:nrow(test_guide_names),
    guide_data_index = test_guide_names$i,
    guide_to_gene = test_guide_names$mapping,
    dispersion_prior_mean = dispersion_init
    )

  gg_const$bin_alpha_prior = matrix(
    rep(bin_size_prior, gg_const$N),
      nrow = gg_const$N, byrow = TRUE)
  gg_const$N_cutoffs = gg_const$N_bins - 1

  gg_init = list(
    psi = 0.2,
    sigma_gene = 5,
    sigma_guide = 5,
    dispersion = dispersion_init
    )
  gg_init$bin_alpha = gg_const$bin_alpha_prior
  gg_init$cutoffs = matrix(rep(dirichlet_to_normal_bins(gg_init$bin_alpha[1, ]), gg_const$N),
    nrow = gg_const$N, byrow = TRUE)
  q_init = array(0, dim = c(gg_const$N, gg_const$N_guides, N_bins))
  for (n in 1:dim(q_init)[1]) {
    for (g in 1:dim(q_init)[2]) {
      q_init[n, g, ] = bin_size_prior
    }
  }

  gg_init$q = q_init
  gg_init$gene_inclusion = rep(0, gg_const$N_genes)
  gg_init$gene_shift = rep(0, gg_const$N_genes)
  gg_init$guide_shift = rep(0, gg_const$N_guides)
  gg_init$total_shift = rep(0, gg_const$N_guides)

  list(data = gg_data, init = gg_init, const = gg_const, test_guide_names = test_guide_names)
}

#' @export
wb_em_start = function(wo, n_it = 10, n_random_guides = 1000) {
  message('estimating bin sizes')
  control_guides = wo$data$x[, wo$const$nt_data_index, , drop = FALSE]
  target_counts = wo$data$x[, wo$const$guide_data_index, , drop = FALSE]
  if (is.null(wo$const$nt_data_index)) {
    msg = paste0('no control guides, randomly sampling ', n_random_guides, '.')
    message(msg)
    control_guides = target_counts[, sample(1:dim(target_counts)[2], n_random_guides), ,
      drop = FALSE]
  }
  dispersion_hat = estimate_dispersion_from_controls(control_guides, wo$init$bin_alpha)
  # null_cutoffs = t(apply(control_guides, 1, dirichlet_to_normal_bins))
  bin_hat = estimate_bin_sizes(control_guides, wo$init$bin_alpha, dispersion_hat)

  for (i in 1:n_it) {
    message('.', appendLF = FALSE)
    dispersion_hat = estimate_dispersion_from_controls(control_guides, bin_hat)
    # null_cutoffs = t(apply(control_guides, 1, dirichlet_to_normal_bins))
    bin_hat = estimate_bin_sizes(control_guides, bin_hat, dispersion_hat)
  }
  message('')
  # browser()

  message('estimating guide and gene level effects')
  guide_mu  = estimate_guide_mu(target_counts, bin_hat, dispersion_hat)
  gene_mu = estimate_gene_mu(guide_mu, wo$const$guide_to_gene)
  message('ranking based on likelihood ratio')
  test_stat = lrt(target_counts, bin_hat, dispersion_hat, guide_mu, wo$const$guide_to_gene)
  test_stat = dplyr::mutate(test_stat, mapping = as.integer(mapping))

  # recode the mapping
  gm = dplyr::distinct(dplyr::select(wo$test_guide_names, gene, mapping))
  test_stat = dplyr::inner_join(test_stat, gm, by = 'mapping')
  test_stat = dplyr::arrange(test_stat, desc(lambda))

  updates = list(
    dispersion = dispersion_hat,
    bins = bin_hat,
    guide_mu = guide_mu,
    gene_mu = gene_mu,
    test_stat = test_stat,
    order = test_stat$mapping
  )
  wo$init$dispersion = updates$dispersion
  wo$init$bin_alpha = updates$bins
  cutoffs_hat = matrix(apply(updates$bins, 1, dirichlet_to_normal_bins),
    nrow = wo$const$N, byrow = TRUE)
  # cutoffs_hat = t(apply(updates$bins, 1, dirichlet_to_normal_bins))
  wo$init$cutoffs = cutoffs_hat
  n_sig = round(wo$init$psi * wo$const$N_genes)
  wo$init$gene_inclusion[1:length(wo$gene_inclusion)] = 0
  wo$init$gene_inclusion[updates$order[1:n_sig]] = 1
  wo$const$order = updates$order
  wo$updates = updates

  wo
}

#' @param bin_sizes either a vector the length of observed bins to be used for every sample,
#' or, a list the length of the number of samples with the estimated bin size for
#' the observed bins
#' @export
wb_estimate_unobserved = function(counts_array, estimated_fraction, bin_sizes,
  ordering) {

  if (is.numeric(bin_sizes)) {
    bin_sizes = lapply(1:dim(counts_array)[1], function(i) bin_sizes)
  }
  if (is.numeric(estimated_fraction)) {
    estimated_fraction = lapply(1:dim(counts_array)[1], function(i) estimated_fraction)
  }

  unobserved_by_sample = lapply(1:dim(counts_array)[1],
    function(i) {
      estimate_unobserved_per_sample(counts_array[i, , ],
        estimated_fraction[[i]],
        bin_sizes[[i]])
    })
  updated_counts = array(NA, dim = dim(counts_array) + c(0, 0, 1))
  updated_counts[
    1:dim(counts_array)[1],
    1:dim(counts_array)[2],
    1:dim(counts_array)[3]] = counts_array
  tmp = dimnames(counts_array)
  tmp[[3]] = c(tmp[[3]], 'unobserved')
  dimnames(updated_counts) = tmp
  for (i in 1:dim(updated_counts)[1]) {
    updated_counts[i, , 'unobserved'] = as.integer(pmax(unobserved_by_sample[[i]]$unobserved, 0))
  }
  updated_counts[, , ordering]
}

estimate_unobserved_per_sample = function(counts, gfpr, bin_sizes) {
  observed = rowSums(counts)
  n_hat = sum(observed) / sum(bin_sizes)
  gfpr = gfpr / sum(gfpr)
  unobserved = n_hat * gfpr - observed

  lower_bound = observed + 1
  upper_bound = round(2 * unobserved + lower_bound)

  list(
    unobserved = unobserved,
    n_hat = n_hat,
    lower_bound = lower_bound,
    upper_bound = upper_bound)
}
