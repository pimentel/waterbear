#' parallel waterbear inference
#'
#' run waterbear in parallel
#' TODO: add some extraction functions
#' TODO: add example
#'
#' @param wo waterbear data object from `counts_to_wb()`
#' @param cluster cluster object from `parallel::makeCluster()`. see example.
#' @param n_chains total number of MCMC chains. We recommend at least 4.
#' @param n_burnin the number of samples to burn in (discard) at the beginning of each chain. we recommend this be half the total number of MCMC samples.
#' @param n_samples the total number of MCMC samples generated from each chain. note, we will only keep the samples after the burn in period.
#' @param seed either a single number and consecutive numbers will be chosen for each chain or a unique length vector of seeds for each chain.
#' @param nimble_params additional parameters to be passed to nimble NIMBLE.
#' @return a set of MCMC samples.
#' @export
wb_run_parallel = function(
  wo,
  cluster,
  n_chains = 4,
  n_burnin = 20000,
  n_samples = 40000,
  seed = 42,
  thin = 1,
  summary = TRUE,
  nimble_params = list()
  ) {
  if (length(seed) > 1 && length(unique(seed)) != n_chains) {
    stop('number of unique seeds is not equal to the number of chains.')
  }
  if (n_samples - n_burnin <= 0) {
    stop('number of MCMC samples must be greater than number of burnin samples.')
  }

  if (length(seed) == 1) {
    seed = seed:(seed + n_chains - 1)
  }

  wb_lambda = function(s) {
    wb_run_sequential(wo, n_chains = 1, n_burnin = n_burnin,
    n_samples = n_samples,
    seed = s,
    thin = thin,
    summary = summary,
    nimble_params = nimble_params)
  }

  res = parLapply(
    cl = cluster,
    X = seed,
    fun = wb_lambda
  )
  res
}

#' single-core waterbear inference
#'
#' run waterbear sequentially
#' TODO: add some extraction functions
#' TODO: add example
#' @param wo waterbear data object from `counts_to_wb()`
#' @param n_chains total number of MCMC chains. We recommend at least 4.
#' @param n_burnin the number of samples to burn in (discard) at the beginning of each chain. we recommend this be half the total number of MCMC samples.
#' @param n_samples the total number of MCMC samples generated from each chain. note, we will only keep the samples after the burn in period.
#' @param seed either a single number and consecutive numbers will be chosen for each chain or a unique length vector of seeds for each chain.
#' @param nimble_params additional parameters to be passed to nimble NIMBLE.
#' @return a set of MCMC samples.
#' @export
wb_run_sequential = function(
  wo,
  n_chains = 4,
  n_burnin = 20000,
  n_samples = 40000,
  seed = 42,
  thin = 1,
  summary = TRUE,
  nimble_params = list()
  ) {
  if (length(seed) > 1 && length(unique(seed)) != n_chains) {
    stop('number of unique seeds is not equal to the number of chains.')
  }
  if (n_samples - n_burnin <= 0) {
    stop('number of MCMC samples must be greater than number of burnin samples.')
  }

  if (length(seed) == 1) {
    seed = seed:(seed + n_chains - 1)
  }
  n_model = nimbleModel(sample_specific_dispersion_model,
    data = wo$data, constants = wo$const, inits = wo$init)

  n_configuration = configureMCMC(n_model)

  n_configuration$addMonitors(
    c('gene_inclusion',
      'total_shift',
      'guide_shift',
      'gene_shift',
      'dispersion',
      'sample_dispersion',
      'psi'
      ))

  n_mcmc_build = buildMCMC(n_configuration)
  C_n_model = compileNimble(n_model)
  C_n_mcmc = compileNimble(n_mcmc_build, project = n_model)

  samples = runMCMC(
    C_n_mcmc, niter = n_samples, nburnin = n_burnin,
    nchains = n_chains,
    thin = thin,
    setSeed = seed,
    summary = summary)
  samples
}
