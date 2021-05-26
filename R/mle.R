dirichlet_multinomial_null_ll = function(counts, bin_sizes, phi) {
  # make sure that bin_sizes is a probability distribution
  stopifnot(all(apply(bin_sizes, 1, sum) == rep(1, nrow(bin_sizes))))
  a = phi * bin_sizes
  total_counts = apply(counts, c(1, 2), sum)

  # for each replicate, iterate through and compute the log probability
  sum(sapply(1:dim(counts)[1],
    function(n) {
      sum(sapply(1:dim(counts)[2],
        function(g) {
          ddirchmulti(counts[n, g, ], a[n, ], total_counts[n, g], TRUE)
        }))
    }))
}

estimate_dispersion_from_controls = function(counts, bin_sizes) {
  fit = optim(c(1),
    function(x) {
      -dirichlet_multinomial_null_ll(counts, bin_sizes, x)
    }, method = 'Brent', lower = 0.0005, upper = 1000)
  fit$par
}
