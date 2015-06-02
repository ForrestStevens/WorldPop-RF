WorldPop-RF
===========

This is the home for the RandomForest-based production code for WorldPop's population maps.

Prediction of Gridded Population Density
--------------------------------------------------------

The source stored here represent the Python and R scripts backing the [WorldPop](http://worldpop.org) project's current mapping efforts to produce gridded population density estimates of people per ~100 m pixel. These are produced using a Random Forest (RF) model as described in Stevens, et al. ([2015](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0107042)).  The published article contains a description of the RF model and its covariates, their sources and any metadata collected for each covariate.  An estimated prediction weighting layer is used to dasymetrically redistribute the census counts and project counts to match estimated populations based on UN estimates for the final population maps provided by [WorldPop](http://worldpop.org).

###### Stevens, F. R., Gaughan, A. E., Linard, C., & Tatem, A. J. (2015). Disaggregating Census Data for Population Mapping Using Random Forests with Remotely-Sensed and Ancillary Data. [PLOS ONE, 10(2), e0107042. doi:10.1371/journal.pone.0107042](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0107042)

