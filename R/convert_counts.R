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

  nt_guide_names = inner_join(guide_names,
    dplyr::filter(gene_mapping, grepl(control_guide_regex, gene)), by = 'guide')
  if (nrow(nt_guide_names) < 1) {
    stop('couldn\'t find any control guides. check the regular expression in `control_guide_regex` and make sure it finds genes in `gene_mapping`.')
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
    N_nt = nrow(nt_guide_names),
    N_bins = N_bins,
    nt_index = 1:nrow(nt_guide_names),
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
