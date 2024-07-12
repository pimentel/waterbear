# Waterbear: A model for accurate quantification of CRISPR effects in pooled FACS screens

Waterbear is a statistical method for analyzing CRISPR FACS screens.
It uses a hierarchical model and MCMC sampling to infer the posterior and is implemented in R.

To get started:

```r
install.packages(c('NIMBLE', 'tidyverse'))
devtools::install_github('pimentel/waterbear')
```

Then, you can access the vignettes (currently only one) by typing:

```r
devtools::browseVignettes('waterbear')
```

## Questions?

Please post any questions or any bugs to GitHub: https://github.com/pimentel/waterbear/issues.

## Development

The current implementation depends on NIMBLE and thus, the interface leaves a bit to be desired.
There is currently a manual C++ implementaiton under development in the `gibbs` branch.
We expect it to be functional in the near future.
