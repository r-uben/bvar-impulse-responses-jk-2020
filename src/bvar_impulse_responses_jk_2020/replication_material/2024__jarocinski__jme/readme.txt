Replication package for:
Marek Jarocinski "Estimating the Fed's Unconventional Monetary Policy Shocks," Journal of Monetary Economics
2024-01-05

Folder structure:
daily/ - data for daily local projections: source data and a Stata script that merges them into a single dataset daily.csv
data/ - Matlab script that creates work/data.csv, which contains MP1, ONRUN2, ONRUN10, SP500 from the GSS dataset
results/ - results, Matlab scripts in the other folders save their results here
work/ - estimate model (1) by Maximum Likelihood and Bayesian methods, produce Figure 3
work1_plots/ - produce Figures 1, 2, 4 and Table 1
work2_lp/ - estimate local projections and produce Figure 5

To replicate the results in the paper:
1. Run work/main1_mh.m to estimate the shocks with both maximum likelihood and Bayesian methods. The script will save its output in the folder results/.
2. Run work/main3_plots.m to produce Figure 3 (and other results). The script will save the save its output in the folder results/.
3. Switch to work1_plots/ and run figure1.m, figure2.m, figure4.m, table1.m.
4. Switch to work2_lp/ and run main2.m to produce Figure 5.
Note: Steps 3 and 4 rely on the content of folder results/

These codes were tested using Matlab R2023b.