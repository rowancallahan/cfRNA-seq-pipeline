import pandas as pd
import pystan
import numpy as np


def wrangle_BCL6(rna_f, tfa_f, meta_f, c_line):
    rna = pd.read_csv(rna_f, sep='\t', index_col=0)
    meta = pd.read_csv(meta_f, sep='\t')
    tfa = pd.read_csv(tfa_f, sep='\t', index_col=0)
    meta_filt = meta[(meta['Cell-line'] == c_line)]

    tfa_filt = tfa.loc[meta_filt.SampleID, :].T
    rna_filt = rna.loc[:, meta_filt.SampleID]

    rna_melt = rna_filt.unstack().reset_index()
    rna_melt.columns = ['Sample','Gene','RNA']

    tfa_melt = tfa_filt.unstack().reset_index()
    tfa_melt.columns = ['Sample','Gene','TFA']
    tfa_melt = tfa_melt[tfa_melt.TFA != 0]
    tfa_melt['tf_idx'] = pd.Categorical(tfa_melt['Gene'].astype(str)).codes
    tfa_melt.loc[:, 'tf_idx'] = tfa_melt.loc[:, 'tf_idx'] + 1

    expression_frame = pd.merge(rna_melt, tfa_melt, how='left', on=['Sample', 'Gene'])
    expression_frame['samp_idx'] = pd.Categorical(expression_frame['Sample'].astype(str)).codes
    expression_frame['gene_idx'] = pd.Categorical(expression_frame['Gene'].astype(str)).codes
    expression_frame['gene_idx'] = expression_frame['gene_idx'].astype(np.int64)
    expression_frame['samp_idx'] = expression_frame['samp_idx'].astype(np.int64)
    expression_frame = expression_frame.fillna(0)
    expression_frame.tf_idx = expression_frame.tf_idx.astype(np.int64)
    conditions = [(expression_frame['TFA'] > 0.0), (expression_frame['TFA'] == 0), (expression_frame['TFA'] < 0.0)]
    choices = [1, 0, -1]
    expression_frame['MoA'] = np.select(conditions, choices)

    return expression_frame


def hierarchical_intercept_model(expression_frame):
    hierarchical_intercept = """
    data {
      int<lower=0> J; 
      int<lower=0> N; 
      int<lower=1,upper=J> Gene[N];
      vector[N] u;
      vector[N] x;
      vector[N] y;
    } 
    parameters {
      vector[J] a;
      vector[2] b;
      real mu_a;
      real<lower=0,upper=100> sigma_a;
      real<lower=0,upper=100> sigma_y;
    } 
    transformed parameters {
      vector[N] y_hat;
      vector[N] m;

      for (i in 1:N) {
        m[i] <- a[Gene[i]] + u[i] * b[1];
        y_hat[i] <- m[i] + x[i] * b[2];
      }
    }
    model {
      mu_a ~ normal(0, 1);
      a ~ normal(mu_a, sigma_a);
      b ~ normal(0, 1);
      y ~ normal(y_hat, sigma_y);
    } 
    """
    hierarchical_intercept_data = {'N': expression_frame.shape[0],
                                   'J': expression_frame.tf_idx.max() + 1,
                                   'Gene': expression_frame.tf_idx + 1,  # Stan counts starting at 1
                                   'u': expression_frame.RNA,
                                   'x': expression_frame.MoA,
                                   'y': expression_frame.TFA}

    hierarchical_intercept_fit = pystan.stan(model_code=hierarchical_intercept, data=hierarchical_intercept_data,
                                             iter=1000, chains=2)
    return hierarchical_intercept_fit
