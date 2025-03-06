from dataclasses import dataclass
from typing import List, Dict, Optional, Union
import numpy as np
import pandas as pd
from scipy import stats, linalg

@dataclass
class VARPrior:
    """Prior specification for Bayesian VAR
    
    Attributes:
        lags: Number of lags in the VAR
        tightness: Overall tightness of the Minnesota prior (lambda_1)
        decay: Lag decay parameter (lambda_2)
        minnesota_mean: Optional mean vector for Minnesota prior
    """
    lags: int = 12
    tightness: float = 0.2
    decay: float = 1.0
    minnesota_mean: Optional[np.ndarray] = None

@dataclass
class GibbsSettings:
    """Settings for Gibbs sampler
    
    Attributes:
        n_draws: Number of draws to keep
        burnin: Number of initial draws to discard
        save_every: Save every nth draw
        compute_marglik: Whether to compute marginal likelihood
    """
    n_draws: int = 4000
    burnin: int = 4000
    save_every: int = 4
    compute_marglik: bool = False

class BayesianVAR:
    """Implementation of Bayesian VAR with Minnesota prior and Gibbs sampling
    
    Based on Jarocinski & Karadi (2020) implementation in MATLAB
    """
    
    def __init__(
        self,
        prior: VARPrior,
        gs_settings: GibbsSettings
    ):
        self.prior = prior
        self.gs_settings = gs_settings
        self._validate_settings()
        
    def _validate_settings(self):
        """Validate prior and Gibbs sampler settings"""
        if self.prior.lags < 1:
            raise ValueError("Number of lags must be positive")
        if self.prior.tightness <= 0:
            raise ValueError("Prior tightness must be positive")
        if self.prior.decay <= 0:
            raise ValueError("Decay parameter must be positive")
        if self.gs_settings.n_draws < 1:
            raise ValueError("Number of draws must be positive")
        if self.gs_settings.burnin < 0:
            raise ValueError("Burnin must be non-negative")
        if self.gs_settings.save_every < 1:
            raise ValueError("save_every must be positive")
            
    def prepare_data(self, data: pd.DataFrame) -> Dict[str, np.ndarray]:
        """Transform pandas DataFrame into VAR-ready format
        
        Args:
            data: DataFrame with datetime index and variables as columns
            
        Returns:
            Dictionary containing:
                y: Endogenous variables matrix
                w: Exogenous variables matrix (if any)
                names: List of variable names
        """
        # Ensure data is sorted by date
        data = data.sort_index()
        
        # Convert to numpy array
        y = data.values
        
        # Create lags matrix
        T, N = y.shape
        y_lags = np.zeros((T - self.prior.lags, self.prior.lags * N))
        
        for p in range(self.prior.lags):
            y_lags[:, N*p:N*(p+1)] = y[self.prior.lags-p-1:T-p-1]
            
        return {
            'y': y[self.prior.lags:],  # T x N matrix
            'y_lags': y_lags,          # T x (N*p) matrix
            'names': list(data.columns)
        }
        
    def _setup_minnesota_prior(self, data: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """Set up Minnesota prior
        
        Args:
            data: Dictionary with prepared data
            
        Returns:
            Dictionary with prior matrices
        """
        T, N = data['y'].shape
        K = N * self.prior.lags
        
        # Compute residual variances for scaling
        sigma = np.zeros(N)
        for i in range(N):
            ar = np.zeros((T-1, 2))
            ar[:, 0] = 1  # constant
            ar[:, 1] = data['y'][:-1, i]
            beta = np.linalg.lstsq(ar, data['y'][1:, i], rcond=None)[0]
            sigma[i] = np.std(data['y'][1:, i] - ar @ beta)
            
        # Prior mean
        if self.prior.minnesota_mean is None:
            self.prior.minnesota_mean = np.zeros(N)
            
        # Prior variance matrix
        prior_var = np.zeros((K, N))
        
        for i in range(N):
            for j in range(N):
                for p in range(self.prior.lags):
                    if i == j:
                        prior_var[p*N + j, i] = (self.prior.tightness / ((p+1)**self.prior.decay))
                    else:
                        prior_var[p*N + j, i] = (self.prior.tightness * sigma[i]**2 / 
                                               (sigma[j]**2 * (p+1)**self.prior.decay))
                        
        return {
            'mean': self.prior.minnesota_mean,
            'var': prior_var,
            'sigma': sigma
        }
        
    def estimate(self, data: pd.DataFrame) -> Dict[str, np.ndarray]:
        """Main estimation routine
        
        Args:
            data: DataFrame with datetime index and variables as columns
            
        Returns:
            Dictionary containing estimation results
        """
        # Prepare data
        prepared_data = self.prepare_data(data)
        
        # Setup prior
        prior_matrices = self._setup_minnesota_prior(prepared_data)
        
        # Initialize storage for Gibbs sampling results
        n_save = self.gs_settings.n_draws // self.gs_settings.save_every
        T, N = prepared_data['y'].shape
        K = N * self.prior.lags
        
        beta_draws = np.zeros((n_save, K, N))
        sigma_draws = np.zeros((n_save, N, N))
        
        # TODO: Implement Gibbs sampler
        # This will be the next major component to implement
        
        return {
            'beta_draws': beta_draws,
            'sigma_draws': sigma_draws,
            'data': prepared_data,
            'prior': prior_matrices
        } 

    def compute_impulse_responses(self, beta_draws, sigma_draws, horizon=40, 
                                 identification='cholesky', sign_restrictions=None, 
                                 max_trials=1000):
        """Compute impulse responses based on posterior draws
        
        Args:
            beta_draws: Posterior draws of VAR coefficients
            sigma_draws: Posterior draws of VAR covariance matrices
            horizon: Horizon for impulse responses
            identification: Identification method ('cholesky' or 'sign')
            sign_restrictions: Dictionary with sign restrictions for identified shocks
            max_trials: Maximum number of rotation matrices to try (for sign restrictions)
            
        Returns:
            Array of impulse responses [n_draws, n_vars, n_vars, horizon]
        """
        n_draws = beta_draws.shape[0]
        N = sigma_draws.shape[1]  # Number of variables
        irfs = np.zeros((n_draws, N, N, horizon))
        
        # For each posterior draw
        for d in range(n_draws):
            # Extract coefficients for this draw
            beta = beta_draws[d].reshape(N, -1)  # [N, K]
            sigma = sigma_draws[d]  # [N, N]
            
            # Create companion form matrices
            companion = self._create_companion_matrix(beta, N, self.prior.lags)
            
            if identification == 'cholesky':
                # Cholesky identification
                impact = np.linalg.cholesky(sigma)
                
                # Compute impulse responses
                irfs[d] = self._compute_irf_from_impact(companion, impact, N, horizon)
                
            elif identification == 'sign':
                # Sign restriction identification
                if sign_restrictions is None:
                    raise ValueError("Sign restrictions must be provided for sign identification")
                
                # Try to find rotation matrices satisfying the restrictions
                for _ in range(max_trials):
                    # Start with Cholesky as the base
                    impact = np.linalg.cholesky(sigma)
                    
                    # Generate random rotation matrix for the first two shocks
                    theta = np.random.uniform(0, 2*np.pi)
                    rotation = np.eye(N)
                    rotation[0:2, 0:2] = np.array([
                        [np.cos(theta), -np.sin(theta)],
                        [np.sin(theta), np.cos(theta)]
                    ])
                    
                    # Apply rotation
                    candidate_impact = impact @ rotation
                    
                    # Compute impulse responses for checking restrictions
                    candidate_irf = self._compute_irf_from_impact(companion, candidate_impact, N, horizon)
                    
                    # Check if restrictions are satisfied
                    mp_satisfied = (
                        candidate_irf[0, 0, 0] > 0 and  # ff4_hf responds positively to MP shock
                        candidate_irf[1, 0, 0] < 0       # sp500_hf responds negatively to MP shock
                    )
                    
                    cbi_satisfied = (
                        candidate_irf[0, 1, 0] > 0 and  # ff4_hf responds positively to CB info shock
                        candidate_irf[1, 1, 0] > 0       # sp500_hf responds positively to CB info shock
                    )
                    
                    if mp_satisfied and cbi_satisfied:
                        # Save the IRFs that satisfy the restrictions
                        irfs[d] = candidate_irf
                        break
                else:
                    # If no rotation satisfies restrictions after max_trials
                    print(f"Warning: No rotation matrix found satisfying sign restrictions for draw {d}")
                    # Just use the last candidate, despite not satisfying restrictions
                    irfs[d] = candidate_irf
        
        return irfs

    def _create_companion_matrix(self, beta, N, lags):
        """Create companion form matrix from VAR coefficients
        
        Args:
            beta: VAR coefficients [N, N*lags]
            N: Number of variables
            lags: Number of lags
            
        Returns:
            Companion form matrix
        """
        # Initialize companion matrix
        companion = np.zeros((N * lags, N * lags))
        
        # Fill in the coefficient matrices
        for i in range(lags):
            companion[0:N, i*N:(i+1)*N] = beta[:, i*N:(i+1)*N]
        
        # Fill in the identity matrices
        for i in range(1, lags):
            companion[i*N:(i+1)*N, (i-1)*N:i*N] = np.eye(N)
        
        return companion

    def _compute_irf_from_impact(self, companion, impact, N, horizon):
        """Compute impulse responses from companion form and impact matrix
        
        Args:
            companion: Companion form matrix
            impact: Impact matrix (e.g., Cholesky decomposition of covariance)
            N: Number of variables
            horizon: Horizon for impulse responses
            
        Returns:
            Array of impulse responses [N, N, horizon]
        """
        # Initialize impulse responses
        irf = np.zeros((N, N, horizon))
        
        # Impact responses (first period)
        irf[:, :, 0] = impact
        
        # Propagate the impulses forward
        for h in range(1, horizon):
            # Compute shock response at horizon h
            temp = np.zeros((N * self.prior.lags, N))
            temp[0:N, :] = irf[:, :, h-1]
            
            # Multiply by companion matrix
            temp = companion @ temp
            
            # Store the responses
            irf[:, :, h] = temp[0:N, :]
        
        return irf 