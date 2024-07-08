#' Get guide level summaries
#'
#' This functigon provids both guide and gene level summaries.
#' @param wo a fully initialized waterbear object.
#' @param wb_samples NIMBLE samples in waterbear format
#' @return a list with three components: (1) guides, (2) gene_inclusion estimates, and (3) gene_shift estimates
#' @export
wb_guide_summaries = function(wo, wb_samples) {
  sample_summary = wb_samples$summary$all.chains
  guide_to_guide_mapping = data.frame(guide_data_index = wo$const$guide_data_index,
    wb_guide_index = wo$const$guide_index)
  sample_summary = data.frame(variable = rownames(sample_summary), sample_summary)
  guides = dplyr::filter(sample_summary, grepl('guide_shift', variable))
  guides = dplyr::mutate(guides, wb_guide_index = as.integer(str_extract(variable, '(\\d)+')))
  guides = inner_join(guides, guide_to_guide_mapping, by = 'wb_guide_index')

  guides = inner_join(guides, wo$test_guide_names, by = c('guide_data_index' = 'i'))
  gene_inclusion = dplyr::filter(sample_summary, grepl('gene_inclusion', variable))
  gene_inclusion = dplyr::mutate(gene_inclusion,
    gene_mapping = as.integer(str_extract(variable, '(\\d)+')))
  gene_inclusion = inner_join(gene_inclusion,
    dplyr::distinct(dplyr::select(wo$test_guide_names, -c(i, guide))),
    by = c('gene_mapping' = 'mapping'))
  gene_shift = dplyr::filter(sample_summary, grepl('gene_shift', variable))
  gene_shift = dplyr::mutate(gene_shift,
    gene_mapping = as.integer(str_extract(variable, '(\\d)+')))
  gene_shift = inner_join(gene_shift,
    dplyr::distinct(dplyr::select(wo$test_guide_names, -c(i, guide))),
    by = c('gene_mapping' = 'mapping'))

  list(guides = guides, gene_inclusion = gene_inclusion, gene_shift = gene_shift)
}


#' @export
wb_gene_posterior_mass = function(wo, samples, alpha) {
  stopifnot(is.list(samples$samples))
  # assumes chains are organizing the same exact manner
  # this should happen assuming it is the exact same model
  gene_columns = grep('gene_shift', colnames(samples$samples$chain1))
  n_chains = length(samples$samples)
  n_samples = nrow(samples$samples$chain1)
  tmp_values = vector('numeric', n_chains * n_samples)
  one_chain = vector('numeric', n_samples)
  n_genes = sum(gene_columns)

  gene_summary = lapply(gene_columns,
    function(gene) {
      tmp = sapply(samples$samples,
        function(chain) {
          chain[, gene]
          # }, one_chain)
        }, simplify = FALSE)
      gene_name = colnames(samples$samples$chain1)[gene]
      tmp = unlist(tmp)
      qs = quantile(tmp, probs = c(0 + alpha/2, 1 - alpha/2))
      # gt0 = mean(tmp > .Machine$double.eps)
      # lt0 = mean(tmp < .Machine$double.eps)
      mu = mean(tmp)
      gt0 = mean(tmp >= 0)
      lt0 = mean(tmp <= 0)
      data.frame(nimble_id = gene_name, lower = qs[1], upper = qs[2],
        gt0 = gt0, lt0 = lt0, mu = mu, row.names = NULL)
    })

  gene_summary = dplyr::bind_rows(gene_summary)
  gene_summary = dplyr::mutate(gene_summary, mapping = as.integer(stringr::str_extract(nimble_id, '(\\d)+')))
  gene_summary = inner_join(gene_summary,
    dplyr::distinct(dplyr::select(wo$test_guide_names, -c(i, guide))),
    by = c('mapping' = 'mapping'))
  gene_summary = dplyr::mutate(gene_summary, lfsr = pmin(gt0, lt0))
  dplyr::select(gene_summary, -nimble_id, -mapping)
}


#' Get the posterior of several parameters
#'
#' After the sampling has been done using NIMBLE, you can extract several of the relevant parameters.
#'
#' @param wo a fully initialized waterbear object
#' @param wb_samples a sampled model
#' @param parameter the parameter that you wish to extract: 'guide', 'gene_pip', or 'gene'.
#' @return a matrix where the rows represent the parameter and of the columns represent a specific sample. Note that the samples will always be in the same order, so if you want to compute functions of the samples, you can simply do the appropriate element wise operations. An example might be getting the guide level effect while incorporating the posterior inclusion probability.
#' @export
wb_get_posterior = function(wo, wb_samples, parameter) {
  stopifnot(is.list(wb_samples$samples))
    n_chains = length(wb_samples$samples)
    n_samples = nrow(wb_samples$samples$chain1)
    tmp_values = vector('numeric', n_chains * n_samples)
    one_chain = vector('numeric', n_samples)
  sample_summary = wb_samples$summary$all.chains
  if (parameter == 'guide') {
    n_chains = length(wb_samples$samples)
    n_samples = nrow(wb_samples$samples$chain1)
    tmp_values = vector('numeric', n_chains * n_samples)
    one_chain = vector('numeric', n_samples)
    guide_columns = grep('guide_shift', colnames(wb_samples$samples$chain1),
      value = TRUE)
    n_guides = length(guide_columns)

    posterior = sapply(guide_columns,
      function(guide) {
        tmp = lapply(wb_samples$samples,
          function(chain) {
            chain[, guide]
          })
        unlist(tmp, use.names = FALSE)
      }, USE.NAMES = FALSE)
    posterior = t(posterior)

    tmp = data.frame(
      wb_guide_index = as.integer(str_extract(guide_columns, '(\\d)+')))
    guide_to_guide_mapping = data.frame(guide_data_index = wo$const$guide_data_index,
      wb_guide_index = wo$const$guide_index)
    tmp = dplyr::inner_join(tmp, guide_to_guide_mapping,
      by = 'wb_guide_index')
    tmp = inner_join(tmp,
      wo$test_guide_names, by = c('guide_data_index' = 'i'))
    rownames(posterior) = tmp$guide
  } else if (parameter == 'gene_pip') {
    columns = grep('gene_inclusion', colnames(wb_samples$samples$chain1), value = TRUE)

    posterior = sapply(columns,
      function(guide) {
        tmp = lapply(wb_samples$samples,
          function(chain) {
            chain[, guide]
          })
        unlist(tmp, use.names = FALSE)
      }, USE.NAMES = FALSE)
    posterior = t(posterior)

    tmp = dplyr::select(
      dplyr::arrange(wo$test_guide_names, mapping),
      gene, mapping)
    tmp = dplyr::distinct(tmp)


    rownames(posterior) = tmp$gene
  } else if (parameter == 'gene') {
    columns = grep('gene_shift', colnames(wb_samples$samples$chain1), value = TRUE)

    posterior = sapply(columns,
      function(guide) {
        tmp = lapply(wb_samples$samples,
          function(chain) {
            chain[, guide]
          })
        unlist(tmp, use.names = FALSE)
      }, USE.NAMES = FALSE)
    posterior = t(posterior)

    tmp = dplyr::select(
      dplyr::arrange(wo$test_guide_names, mapping),
      gene, mapping)
    tmp = dplyr::distinct(tmp)

    rownames(posterior) = tmp$gene
  }

  posterior
}
