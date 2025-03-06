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
import os


def load_shock_data(file_path="data/us_shocks.csv"):
    """Load high-frequency shock data from CSV file
    
    Args:
        file_path: Path to the shock data CSV file
        
    Returns:
        pd.DataFrame: Shock data with datetime index
    """
    df = pd.read_csv(file_path)
    
    # Convert year, month to datetime index
    df['date'] = pd.to_datetime(df[['year', 'month']].assign(day=1)) + pd.offsets.MonthEnd(0)
    df = df.set_index('date')
    
    # Filter to keep only the required shock columns
    shock_cols = ['ff4_hf', 'sp500_hf']
    return df[shock_cols]


def prepare_data(use_shock_data=True):
    """Prepare and merge macroeconomic data for VAR analysis"""
    # Load macroeconomic data
    bucket = BucketManager("macroeconomic-data")
    gdp_data = get_gdp_data(bucket)
    sp500_data = get_sp500()
    ff_data = fed_funds()
 
    if use_shock_data:
        # Load high-frequency shock data first
        print("Loading high-frequency shock data...")
        shock_data = load_shock_data()
        print("\nShock data shape:", shock_data.shape)
        print("Shock data date range:", shock_data.index[0], "to", shock_data.index[-1])
        print("Missing values in shock data:")
        print(shock_data.isnull().sum())

        # Restrict sample to periods where we have shock data
        start_date = shock_data.index[0]
        end_date = shock_data.index[-1]
        
        # Filter other data to match shock data sample
        gdp_data = gdp_data.loc[start_date:end_date]
        sp500_data = sp500_data.loc[start_date:end_date]
        ff_data = ff_data.loc[start_date:end_date]
    
    print("\nData shapes before merge:")
    print("GDP data:", gdp_data.shape)
    print("S&P500 data:", sp500_data.shape)
    print("Fed Funds data:", ff_data.shape)
    
    # Merge all data
    data = pd.merge(
        gdp_data,
        sp500_data.resample('ME').last(),
        left_index=True,
        right_index=True,
        how='inner'
    )
    
    print("\nAfter first merge:", data.shape)
    print("Missing values after first merge:")
    print(data.isnull().sum())
    
    data = pd.merge(
        data,
        ff_data,
        left_index=True,
        right_index=True,
        how='inner'
    )
    
    print("\nAfter second merge:", data.shape)
    print("Missing values after second merge:")
    print(data.isnull().sum())
    
    # Take natural log of S&P 500
    data['log_sp500'] = np.log(data['sp500'])
    
    if use_shock_data:
        # Merge with shock data
        data = pd.merge(
            shock_data,
            data,
            left_index=True,
            right_index=True,
            how='inner'
        )
        
        print("\nAfter merging with shock data:", data.shape)
        print("Missing values in final dataset:")
        print(data.isnull().sum())
        print("\nFirst few rows of final dataset:")
        print(data.head())
    
    return data


def run_var_analysis(use_shock_data=True):
    """Run Bayesian VAR analysis on macroeconomic data"""
    start_time = time.time()
    
    # Get the data
    print("Loading data...")
    data = prepare_data(use_shock_data)
    
    # Select and standardize variables
    var_data = pd.DataFrame({
        # Monetary shock variables first (high-frequency)
        'ff4_hf': data['ff4_hf'],
        'sp500_hf': data['sp500_hf'],
        # Then other variables
        'gs1': data['fed_funds'],
        'logsp500': data['log_sp500'],
        'us_rgdp': np.log(data['gdp_real']),  # Take log of GDP
        'us_gdpdef': np.log(data['gdp_nominal'] / data['gdp_real'])  # Take log of GDP deflator
    })
    
    # Standardize all variables
    var_data = (var_data - var_data.mean()) / var_data.std()
    
    # Verify no missing values
    if var_data.isnull().any().any():
        raise ValueError("Dataset contains missing values after preparation")
    
    print("\nData Summary:")
    print(var_data.describe())
    breakpoint()
    
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
            horizon=40,
            identification='sign',  # Use sign restrictions for identification
            sign_restrictions={
                'monetary_policy': {'ff4_hf': '+', 'sp500_hf': '-'},  # MP shock
                'cb_info': {'ff4_hf': '+', 'sp500_hf': '+'}           # CB info shock
            }
        )
        
        # 6. Print results
        print("\nVAR Estimation Results:")
        print(f"Sample period: {var_data.index[0]} to {var_data.index[-1]}")
        print(f"Variables: {', '.join(var_data.columns)}")
        print(f"Number of lags: {prior.lags}")
        print(f"Number of stored draws: {results['beta_draws'].shape[0]}")
        
        # 7. Plot impulse responses only if we have them
        if irfs is not None:
            # Define shock names based on identification
            shock_names = ['Monetary Policy', 'CB Information'] + [f'Shock {i+3}' for i in range(var_data.shape[1]-2)]
            plot_impulse_responses(irfs, var_data.columns, shock_names=shock_names)
        else:
            print("Warning: No impulse responses were generated")
            
    except Exception as e:
        print(f"Error computing impulse responses: {str(e)}")
        irfs = None
    
    elapsed_time = time.time() - start_time
    print(f"\nTotal execution time: {elapsed_time:.2f} seconds")
    
    return results, irfs

def plot_impulse_responses(irfs, variable_names, shock_names=None, shock_indices=[0, 1], response_indices=None):
    """Plot impulse responses to identified shocks
    
    Args:
        irfs: Impulse response function array [n_draws, n_vars, n_vars, horizon]
        variable_names: Names of variables in the VAR
        shock_names: Names of identified shocks (defaults to variable names)
        shock_indices: Indices of shocks to plot
        response_indices: Indices of response variables to plot (defaults to all)
    """
    n_draws, N, _, horizon = irfs.shape
    
    if shock_names is None:
        shock_names = variable_names
    
    if response_indices is None:
        response_indices = range(N)
    
    # Compute median and percentiles
    median_irfs = np.median(irfs, axis=0)
    lower_irfs = np.percentile(irfs, 16, axis=0)
    upper_irfs = np.percentile(irfs, 84, axis=0)
    
    # Create figure grid based on shock_indices and response_indices
    n_shocks = len(shock_indices)
    n_responses = len(response_indices)
    fig, axes = plt.subplots(n_responses, n_shocks, figsize=(4*n_shocks, 3*n_responses))
    
    # Handle single row or column case
    if n_responses == 1 and n_shocks == 1:
        axes = np.array([[axes]])
    elif n_responses == 1:
        axes = axes.reshape(1, -1)
    elif n_shocks == 1:
        axes = axes.reshape(-1, 1)
    
    # Plot each response to each shock
    for i, resp_idx in enumerate(response_indices):
        for j, shock_idx in enumerate(shock_indices):
            ax = axes[i, j]
            
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
            if i == 0:
                ax.set_title(f"{shock_names[shock_idx]} Shock", fontsize=12)
            if j == 0:
                ax.set_ylabel(f"{variable_names[resp_idx]}", fontsize=12)
            if i == n_responses - 1:
                ax.set_xlabel("Horizon (months)", fontsize=10)
    
    plt.tight_layout()
    plt.savefig("impulse_responses.pdf")
    plt.show()

if __name__ == "__main__":
    results, irfs = run_var_analysis(use_shock_data=True) 