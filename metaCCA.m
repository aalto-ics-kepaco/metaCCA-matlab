%%%% metaCCA %%%

% 15/07/2015 Anna Cichonska
% Helsinki Institute for Information Technology HIIT, Department of Computer Science, Aalto University, Finland
% Institute for Molecular Medicine Finland FIMM, University of Helsinki, Finland
% anna.cichonska@helsinki.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.


% An example of the meta-analysis of 2 studies: 10 traits and 5 SNPs. 
% Test data will be provided here by September 2015
% Univariate summary statistics [100 000 SNPs x 10 traits]
%   S_XY_full_raw_study1 = load();
%   S_XY_full_raw_study2 = load();
% Standard errors SE
%   se_study1 = load();
%   se_study2 = load()
% Number of individuals N
%   N1 = 1000;
%   N2 = 2000;
% Genotypic correlation matrices estimated e.g. from the 1000Genomes database
% (an example of estimating C_XX based on the data for a subset 
%  of individuals from the 1000Genomes database will be provided)
% [5 SNPs x 5 SNPs]
%   S_XX_study1 = load(); 
%   S_XX_study2 = load(); 


% Normalizing regression coefficients 
S_XY_full_study1 = normalize_Sxy(S_XY_full_raw_study1, se_study1, N1);
S_XY_full_study2 = normalize_Sxy(S_XY_full_raw_study2, se_study2, N2);


% Estimating phenotypic correlation matrices from S_XY's
S_YY_study1 = estimate_Syy(S_XY_full_study1);
S_YY_study2 = estimate_Syy(S_XY_full_study2);


% Extracting summary statistics for 5 SNPs (and 10 traits)
S_XY_study1 = S_XY_full_study1( 1:5, :);
S_XY_study2 = S_XY_full_study2( 1:5, :);


% Pooling covariance matrices of the same type
[C_XX, C_YY, C_XY, N_tot] = pool_cov(S_XX_study1,S_XX_study2,  S_YY_study1,S_YY_study2,   S_XY_study1,S_XY_study2,  N1,N2);


% Bulding a full covariance matrix
full_cov  =  [C_XX, C_XY; C_XY', C_YY];


% Shrinkage of the full covariance matrix
[C_XX_out, C_XY_out, C_YY_out] = shrinkPSD(full_cov, size(C_XX,1));


% Canonical Correlation Analysis (CCA)
% Genotype-phenotype association result
[r, ~, ~, ~, ~, ~, pval] = my_cca(C_XX_out, C_YY_out, C_XY_out, N_tot); 

