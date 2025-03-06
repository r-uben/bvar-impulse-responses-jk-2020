For maximum likelihood estimation:
  Run main0_maxlik.m

For Bayesian estimation:
  1. Run main1_mh.m to run the Metropolis-Hastings chain.
  2. Optionally: run main2_normalize.m to normalize the draws in the chain.
     This is only needed if the chain visits multiple modes (visible regime switches in the traceplot of the draws)
  3. Run main3_plots.m to plot the results.

The data in 'data.csv' are taken from the dataset of Gurkaynak, Sack and Swanson (2005, GSS) updated until June 2019 by Gurkaynak, Karasoy-Can and Lee (2021, GKL) and available in their replication files at \url{http://www.bilkent.edu.tr/~refet/GKL_replication.zip}. To see how 'data.csv' is obtained see the script 'data/make_data.m'.
The file 'data.csv' is provided here for convenience, please cite GSS or GKL as the data source.

