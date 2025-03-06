"""
Bayesian Vector Autoregression (BVAR) implementation based on Jarocinski & Karadi (2020)
"""
from dataclasses import dataclass
from typing import Dict, List, Optional, Union
import numpy as np
import pandas as pd
from scipy import stats, linalg
from scipy import sparse

@dataclass
class VARPrior:
    """Prior specification for Bayesian VAR
    
    Attributes:
        lags: Number of lags in the VAR
        tightness: Overall tightness of the Minnesota prior (lambda_1)
        decay: Lag decay parameter (lambda_2)
        exog_std: Standard deviation for exogenous variables
        minnesota_mean: Optional mean vector for Minnesota prior
    """
    lags: int = 12
    tightness: float = 0.2
    decay: float = 1.0
    exog_std: float = 1e5
    minnesota_mean: Optional[np.ndarray] = None

@dataclass
class GibbsSettings:
    """Settings for Gibbs sampler
    
    Attributes:
        n_draws: Number of draws to keep
        burnin: Number of initial draws to discard
        save_every: Save every nth draw
        compute_marglik: Whether to compute marginal likelihood
        verbose: Print progress information
    """
    n_draws: int = 4000
    burnin: int = 4000
    save_every: int = 4
    compute_marglik: bool = False
    verbose: bool = True

class BayesianVAR:
    """Bayesian Vector Autoregression with Minnesota prior and Gibbs sampling
    
    Implementation based on Jarocinski & Karadi (2020) MATLAB code.
    """
    
    def __init__(
        self,
        prior: VARPrior,
        gs_settings: GibbsSettings,
        n_monetary: int
    ):
        """Initialize the BVAR model
        
        Args:
            prior: Prior specification
            gs_settings: Gibbs sampler settings
            n_monetary: Number of monetary variables
        """
        self.prior = prior
        self.gs_settings = gs_settings
        self.n_monetary = n_monetary
        self._validate_settings()
    
    def _validate_settings(self) -> None:
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
            Dictionary with prepared data
        """
        # Ensure data is sorted by date
        data = data.sort_index()
        
        # Convert to numpy array
        y = data.values
        
        # Get dimensions
        T_full, N = y.shape
        Nm = self.n_monetary
        Ny = N - Nm
        P = self.prior.lags
        
        # Extract initial lags
        Y0 = y[:P, :]  # Initial lags
        
        # Create lagged data matrix
        y_lags = np.zeros((T_full - P, P * N))
        for p in range(P):
            y_lags[:, p*N:(p+1)*N] = y[P-p-1:T_full-p-1, :]
        
        # No exogenous variables for now
        X = y_lags
        Y = y[P:, :]
        
        n_m = np.arange(Nm)
        n_y = np.arange(Nm, N)
        
        return {
            'T': T_full - P,
            'N': N,
            'Nm': Nm,
            'Ny': Ny,
            'P': P,
            'Y0': Y0,
            'X': X,
            'Y': Y,
            'n_m': n_m,
            'n_y': n_y,
            'names': list(data.columns),
            'K': P * N  # Number of regressors per equation
        }
    
    def _compute_residual_variances(self, Y: np.ndarray, P: int) -> np.ndarray:
        """Compute residual variances from univariate autoregressions
        
        Args:
            Y: Data matrix (T x N)
            P: Number of lags to use
            
        Returns:
            Vector of residual standard deviations
        """
        T, N = Y.shape
        sigma = np.zeros(N)
        
        for n in range(N):
            y_n = Y[:, n]
            valid_idx = ~np.isnan(y_n)
            
            if np.sum(valid_idx) > P + 2:  # Need enough data
                y_valid = y_n[valid_idx]
                # Simple AR(1) for estimating standard deviation
                ar = np.zeros((len(y_valid)-1, 2))
                ar[:, 0] = 1  # constant
                ar[:, 1] = y_valid[:-1]
                beta = np.linalg.lstsq(ar, y_valid[1:], rcond=None)[0]
                sigma[n] = np.std(y_valid[1:] - ar @ beta)
            else:
                sigma[n] = 1.0  # Default if not enough data
                
        return sigma
    
    def _setup_minnesota_prior(self, data: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """Set up Minnesota prior following the original MATLAB implementation
        
        Args:
            data: Prepared VAR data
            
        Returns:
            Dictionary with prior matrices
        """
        N = data['N']
        Nm = data['Nm']
        Ny = data['Ny']
        P = data['P']
        Y = data['Y']
        K = data['K']
        
        # Compute residual variances for scaling
        sigma = self._compute_residual_variances(Y, P)
        
        if self.gs_settings.verbose:
            print(f"Residual standard deviations: {sigma}")
        
        # Create matrices for Minnesota prior
        # Following VAR_withiid1kf.m implementation
        
        # Create P x 1 vector of lag decays
        p_decay = np.power(np.arange(1, P+1), -self.prior.decay).reshape(-1, 1)
        
        # Create N x N matrix of ones
        ones_matrix = np.ones((N, N))
        
        # Kronecker product to get P*N x N matrix
        temp1 = np.kron(p_decay, ones_matrix)
        
        # Create P*N x N matrix where each row is sigma
        temp2 = np.zeros((P*N, N))
        for i in range(P*N):
            temp2[i, :] = sigma
        
        # Create P*N x N matrix for sigma_j^-1
        temp3 = np.zeros((P*N, N))
        for p in range(P):
            for i in range(N):
                for j in range(N):
                    temp3[p*N + i, j] = 1.0/sigma[j]
        
        # Element-wise multiplication
        Q0 = self.prior.tightness * np.multiply(np.multiply(temp1, temp2), temp3)
        
        # Square the values
        Q0 = np.square(Q0)
        
        # Drop the equations for monetary variables
        Q0 = Q0[:, Nm:]
        
        # Create diagonal matrices
        Q_diag = Q0.flatten()
        Qinv_diag = 1.0/Q_diag
        
        Q = sparse.diags(Q_diag, 0, shape=(K*Ny, K*Ny))
        Qinv = sparse.diags(Qinv_diag, 0, shape=(K*Ny, K*Ny))
        
        # Prior mean
        B = np.zeros((K, Ny))
        
        return {
            'Q': Q,
            'Qinv': Qinv,
            'B': B,
            'sigma': sigma
        }
    
    def estimate(self, data_df: pd.DataFrame) -> Dict[str, np.ndarray]:
        """Estimate the BVAR model using Gibbs sampling
        
        Args:
            data_df: DataFrame with datetime index and variables as columns
            
        Returns:
            Dictionary with estimation results
        """
        # Prepare data
        data = self.prepare_data(data_df)
        
        if self.gs_settings.verbose:
            print("\nData dimensions:")
            print(f"Total variables: {data['N']}")
            print(f"Monetary variables: {data['Nm']}")
            print(f"Other variables: {data['Ny']}")
            print(f"Sample size: {data['T']}")
        
        # Setup prior
        prior_matrices = self._setup_minnesota_prior(data)
        
        # Initialize storage for Gibbs sampling results
        n_save = self.gs_settings.n_draws // self.gs_settings.save_every
        
        N = data['N']
        K = data['K']
        T = data['T']
        Nm = data['Nm']
        Ny = data['Ny']
        X = data['X']
        Y = data['Y']
        n_m = data['n_m']
        n_y = data['n_y']
        
        beta_draws = np.zeros((n_save, K, N))
        sigma_draws = np.zeros((n_save, N, N))
        
        # Initialize parameters
        BB = np.zeros((K, N))
        Sigma = np.eye(N)
        
        # Prior for Sigma
        prior_S = np.eye(N)
        prior_v = N + 2
        
        # Gibbs sampler
        for draw in range(self.gs_settings.burnin + self.gs_settings.n_draws):
            # Draw Sigma conditional on B
            U = Y - X @ BB
            S_post = U.T @ U + prior_S
            v_post = T + prior_v
            Sigma = stats.invwishart.rvs(df=v_post, scale=S_post)
            
            # Draw B conditional on Sigma
            Csig = np.linalg.cholesky(Sigma)
            SigmaYY1inv = np.linalg.inv(Csig[n_y][:, n_y].T @ Csig[n_y][:, n_y])
            
            # Kronecker product for posterior precision
            XX = np.kron(SigmaYY1inv, X.T @ X)
            
            # Posterior precision matrix
            A = prior_matrices['Qinv'] + XX
            
            # Adjusted dependent variable
            yst = Y[:, n_y] - Y[:, n_m] @ np.linalg.inv(Csig[n_m][:, n_m].T) @ Csig[n_y][:, n_m].T
            
            # Posterior mean times precision
            a = X.T @ yst @ SigmaYY1inv
            
            # Convert to dense array for Cholesky
            if hasattr(A, 'toarray'):
                A_dense = A.toarray()
            elif isinstance(A, np.matrix):
                A_dense = np.array(A)
            else:
                A_dense = A
            
            # Cholesky decomposition of posterior precision
            C = np.linalg.cholesky(A_dense)
            
            # Draw from posterior
            B = np.linalg.solve(C, np.linalg.solve(C.T, a.flatten()) + np.random.randn(K*Ny))
            
            # Reshape to K x Ny
            B = B.reshape(K, Ny)
            
            # Create full coefficient matrix
            BB = np.zeros((K, N))
            BB[:, Nm:] = B
            
            # Store draws
            if draw >= self.gs_settings.burnin and (draw - self.gs_settings.burnin) % self.gs_settings.save_every == 0:
                idx = (draw - self.gs_settings.burnin) // self.gs_settings.save_every
                beta_draws[idx] = BB
                sigma_draws[idx] = Sigma
                
            # Print progress
            if self.gs_settings.verbose and draw % 100 == 0:
                print(f"Completed draw {draw}/{self.gs_settings.burnin + self.gs_settings.n_draws}")
        
        return {
            'beta_draws': beta_draws,
            'sigma_draws': sigma_draws,
            'data': data,
            'prior': prior_matrices
        }
    
    def compute_impulse_responses(self, beta_draws: np.ndarray, sigma_draws: np.ndarray, 
                                 horizon: int = 40, identification: str = 'chol') -> np.ndarray:
        """Compute impulse response functions"""
        if beta_draws is None or sigma_draws is None:
            print("Warning: Cannot compute impulse responses - draws are None")
            return None
        
        n_draws, K, N = beta_draws.shape
        P = self.prior.lags
        irfs = np.zeros((n_draws, N, N, horizon))
        
        try:
            for d in range(n_draws):
                # Reshape VAR coefficients to companion form
                companion = np.zeros((N * P, N * P))
                companion[:N, :N*P] = beta_draws[d].T
                
                for i in range(1, P):
                    companion[i*N:(i+1)*N, (i-1)*N:i*N] = np.eye(N)
                
                # Cholesky identification
                if identification == 'chol':
                    chol = np.linalg.cholesky(sigma_draws[d])
                else:
                    raise ValueError(f"Identification scheme '{identification}' not implemented")
                
                # Compute impulse responses
                irf = np.zeros((N, N, horizon))
                irf[:, :, 0] = chol
                
                for h in range(1, horizon):
                    # Propagate the shock through the system
                    impact = np.zeros((N * P, N))
                    impact[:N, :] = irf[:, :, h-1]
                    
                    # Next period effect
                    next_impact = companion @ impact
                    irf[:, :, h] = next_impact[:N, :]
                
                irfs[d] = irf
                
            return irfs  # Explicit return
        
        except Exception as e:
            print(f"Error in impulse response calculation: {str(e)}")
            return None 