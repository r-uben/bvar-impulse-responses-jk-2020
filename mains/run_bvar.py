"""
Run Bayesian VAR analysis on macroeconomic data
"""
from bvar_impulse_responses_jk_2020.models.bayesian_var import BayesianVAR, VARPrior, GibbsSettings
from bvar_impulse_responses_jk_2020.aws_manager.bucket_manager import BucketManager
from mains.get_data import get_gdp_data, get_sp500, fed_funds
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time


def prepare_data():
    """Prepare and merge macroeconomic data for VAR analysis
    
    This function:
    1. Loads GDP data from S3 bucket
    2. Loads S&P 500 data from external source 
    3. Loads Federal Funds rate data
    4. Merges all data to monthly frequency using left joins
    5. Takes log of S&P 500 index
    
    Returns:
        pd.DataFrame: Combined dataset with GDP, S&P 500, and Fed Funds rate
                     at monthly frequency
    """
    # Load data from different sources
    bucket = BucketManager("macroeconomic-data")
    gdp_data = get_gdp_data(bucket)
    sp500_data = get_sp500()
    ff_data = fed_funds()
    
    # 2. Prepare the data
    print("Preparing data...")
    # Merge GDP with S&P 500 data, resampling S&P to month-end frequency
    data = pd.merge(
        gdp_data,
        sp500_data.resample('ME').last(),
        left_index=True,
        right_index=True,
        how='left'
    )
    
    # Merge result with Federal Funds rate data
    data = pd.merge(
        data,
        ff_data,
        left_index=True,
        right_index=True,
        how='left'
    )
    
    # Take natural log of S&P 500 for analysis
    data['log_sp500'] = np.log(data['sp500'])
    return data

def run_var_analysis():
    """Run Bayesian VAR analysis on macroeconomic data"""
    
    start_time = time.time()
    
    # 1. Get the data
    print("Loading data...")
    data = prepare_data()
    
    
    # Select variables matching the original paper's specification
    # Following main1.m: mnames = {'ff4_hf','sp500_hf'} for US baseline
    var_data = pd.DataFrame({
        # Monetary variables first (same as in main1.m)
        'ff4_hf': data['fed_funds'],
        'sp500_hf': data['log_sp500'],
        # Then other variables (same as in main1.m for us1)
        'gs1': data['fed_funds'],
        'logsp500': data['log_sp500'],
        'us_rgdp': data['gdp_real'],
        'us_gdpdef': data['gdp_nominal'] / data['gdp_real']
    })
    
    # Handle missing values
    var_data = var_data.ffill().bfill()
    
    # Print data summary
    print("\nData Summary:")
    print(var_data.describe())
    
    # 3. Set up the VAR with the paper's specifications
    print("\nSetting up BVAR model...")
    prior = VARPrior(
        lags=12,
        tightness=0.2,
        decay=1.0
    )
    
    gs_settings = GibbsSettings(
        n_draws=100,  # Reduced for testing
        burnin=100,   # Reduced for testing
        save_every=4,
        verbose=True
    )
    
    # Following main1.m: prior.Nm = length(mnames) = 2
    bvar = BayesianVAR(prior, gs_settings, n_monetary=2)
    
    # 4. Run the estimation
    print("\nEstimating BVAR model...")
    results = bvar.estimate(var_data)
    
    # 5. Compute impulse responses with error handling
    print("\nComputing impulse responses...")
    try:
        irfs = bvar.compute_impulse_responses(
            results['beta_draws'], 
            results['sigma_draws'],
            horizon=40
        )
        
        # 6. Print results
        print("\nVAR Estimation Results:")
        print(f"Sample period: {var_data.index[0]} to {var_data.index[-1]}")
        print(f"Variables: {', '.join(var_data.columns)}")
        print(f"Number of lags: {prior.lags}")
        print(f"Number of stored draws: {results['beta_draws'].shape[0]}")
        
        # 7. Plot impulse responses only if we have them
        if irfs is not None:
            plot_impulse_responses(irfs, var_data.columns)
        else:
            print("Warning: No impulse responses were generated")
            
    except Exception as e:
        print(f"Error computing impulse responses: {str(e)}")
        irfs = None
    
    elapsed_time = time.time() - start_time
    print(f"\nTotal execution time: {elapsed_time:.2f} seconds")
    
    return results, irfs

def plot_impulse_responses(irfs, variable_names, shock_idx=0, response_indices=None):
    """Plot impulse responses to a shock"""
    n_draws, N, _, horizon = irfs.shape
    
    if response_indices is None:
        response_indices = range(N)
    
    # Compute median and percentiles
    median_irfs = np.median(irfs, axis=0)
    lower_irfs = np.percentile(irfs, 16, axis=0)
    upper_irfs = np.percentile(irfs, 84, axis=0)
    
    # Create figure
    fig, axes = plt.subplots(len(response_indices), 1, figsize=(10, 2*len(response_indices)))
    if len(response_indices) == 1:
        axes = [axes]
    
    # Plot each response
    for i, resp_idx in enumerate(response_indices):
        ax = axes[i]
        
        # Plot median
        ax.plot(range(horizon), median_irfs[resp_idx, shock_idx, :], 'b-', linewidth=2)
        
        # Plot confidence bands
        ax.fill_between(
            range(horizon), 
            lower_irfs[resp_idx, shock_idx, :], 
            upper_irfs[resp_idx, shock_idx, :],
            color='b', alpha=0.2
        )
        
        # Add horizontal line at zero
        ax.axhline(y=0, color='k', linestyle='-', alpha=0.2)
        
        # Set title and labels
        ax.set_title(f"Response of {variable_names[resp_idx]} to {variable_names[shock_idx]} shock")
        ax.set_xlabel("Horizon")
        ax.set_ylabel("Response")
    
    plt.tight_layout()
    plt.savefig("impulse_responses.pdf")
    plt.show()

if __name__ == "__main__":
    results, irfs = run_var_analysis() 