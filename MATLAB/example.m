%% An example of the meta-analysis of 2 studies using metaCCA and metaCCA+. 
%% Test input files are provided and described below.

% 15/07/2015 Anna Cichonska
% Helsinki Institute for Information Technology HIIT, Department of Computer Science, Aalto University, Finland
% Institute for Molecular Medicine Finland FIMM, University of Helsinki, Finland
% anna.cichonska@helsinki.fi


clear all; close all; clc


%% TEST DATA

%% Number of individuals in each study
N1 = 1000;
N2 = 2000;



%% Univariate summary statistics S_XY_full for estimating phenotypic correlation structure S_YY.

%  Here, [1000 SNPs x 10 traits].
%  In practice, summary statistics of at least one chromosome should be used 
%  in order to ensure good quality of S_YY estimate.
%
%  FORMAT: tab delimited text file containing the following columns:
%          - SNP id, e.g., position or rs_id (type: string or numeric);
%          - allele 0 (type: string);
%          - allele 1 (type: string);
%          - then, 2 columns for each trait: in turn, a column with linear regression coefficients
%            beta, a column with corresponding standard errors SE (type: numeric).
%
%          The file needs to contain a following header line:
%          first 3 entries must be "SNP_id", "allele_0", "allele_1",
%          the remaining ones must correspond to trait ids: 
%          "traitID_b", "traitID_se".
%          Please include "_b"/"_se" after the trait id, as in the provided files.
%          Do not use underscores "_" in trait ids outside "_b"/"_se"
%          in order for the ids to be processed correctly.    
%
% The data can be read to Matlab using Matlab's 'importdata' function.
% The resulting type of the array should be "structure array".

S_XY_full_study1 = importdata('S_XY_full_study1.txt');
S_XY_full_study2 = importdata('S_XY_full_study2.txt'); 

% class(S_XY_full_study1)



%% Univariate summary statistics S_XY corresponding to the variables to be included in the analysis.

%  Here, [10 SNPs x 10 traits].
%
%  FORMAT: the same as above ('S_XY_full').

S_XY_study1 = importdata('S_XY_study1.txt');
S_XY_study2 = importdata('S_XY_study2.txt');



%% Genotypic correlation matrices S_XX estimated, e.g., from the 1000Genomes database.

%  These matrices are needed ONLY in case of multi-SNP analysis.
%  Here, [10 SNPs x 10 SNPs].
%
%  FORMAT: tab delimited text file.
%          The first column needs to contain SNP ids (type: string or numeric).
%          No header line.
%          Genetic variants in this matrix need to correspond to the ones given above in S_XY.

S_XX_study1 = importdata('S_XX_study1.txt'); 
S_XX_study2 = importdata('S_XX_study2.txt'); 





%% Estimation of phenotypic correlation structure S_YY. 

%  Function 'estimate_Syy'.
%  INPUT:  
%          - a matrix containing univariate summary statistics in the format described above.
%  OUTPUT: 
%          - phenotypic correlation structure; the first column contains trait ids.
%
% The function can be used no matter if the univariate analysis has been performed
% on standardised or non-standardised data.
%
% Please keep in mind that in practice, summary statistics of at least one chromosome 
% should be used in order to ensure good quality of S_YY estimate.

 
S_YY_study1 = estimate_Syy(S_XY_full_study1);
S_YY_study2 = estimate_Syy(S_XY_full_study2);


% % Extracting phenotypic correlation matrix (without trait ids):
% S_YY_1 = cell2mat(S_YY_study1(:, 2:end));
% S_YY_2 = cell2mat(S_YY_study2(:, 2:end));
% % Visualization:
% hcol = [-1 1];
% figure; imagesc(S_YY_1, hcol)
% figure; imagesc(S_YY_2, hcol)





%% metaCCA/metaCCA+.

% Function 'metaCCA'      runs metaCCA, 
% function 'metaCCAplus'  runs metaCCA+.
% Both functions require the same inputs, and they have the 
% same output format.
% They accept a varying number of inputs, depending on the
% number of studies being meta-analysed, and the type of the analysis. 
% Analysis of a single GWAS is also possible.

% By default, single-SNP analysis is performed; each given SNP is analysed
% in turn against all given phenotypic features.


% INPUT
%
% Types of REQUIRED inputs, in order:
%
% - Number of studies analysed.
% Next, for each study:
% - Univariate summary statistics corresponding to variables to be analysed;
% - An information if the univariate analysis has been performed
%   on standardised ('1') or non-standardised ('0') data; please keep
%   in mind that most likely the data were not standardised (option '0' 
%   should be used);
% - Phenotypic correlation matrix estimated using 'estimate_Syy' function;
% - Number of individuals.
%
%
% Subsequent OPTIONAL inputs, in order:
%
% 1) Single-SNP analysis of only one selected SNP
% - '1' (indicator of the analysis type);
% - selected SNP id.
%
% 2) Multi-SNP analysis  
% - '2' (indicator of the analysis type);
% - List of SNP ids to be analysed jointly;
% - Genotypic correlation matrix
%   (corresponding to SNPs in the matrix containing summary statistics,
%   and in the same order).




%% Default single-SNP analysis

% metaCCA
metaCCA_result1     = metaCCA(     2,                        ...
                                   S_XY_study1, S_XY_study2, ...
                                   0, 0,                     ...
                                   S_YY_study1, S_YY_study2, ...
                                   N1, N2 );
% metaCCA+
metaCCAplus_result1 = metaCCAplus( 2,                        ...
                                   S_XY_study1, S_XY_study2, ...
                                   0, 0,                     ...
                                   S_YY_study1, S_YY_study2, ...
                                   N1, N2 );
                              
% OUTPUT 
%   A matrix containing three columns:
%   1) SNP id,
%   2) leading canonical correlation value,
%   3) -log10(p-value).
    

                 

%% Single-SNP analysis of one selected SNP

% Provide SNP id in the following string format:
  SNP_id = 'rs80';
% SNP_id = '80';


metaCCA_result2     =     metaCCA( 2,                        ...
                                   S_XY_study1, S_XY_study2, ...
                                   0, 0,                     ...
                                   S_YY_study1, S_YY_study2, ...
                                   N1, N2,                   ...
                                   1,                        ...
                                   SNP_id );
               
metaCCAplus_result2 = metaCCAplus( 2,                        ...
                                   S_XY_study1, S_XY_study2, ...
                                   0, 0,                     ...
                                   S_YY_study1, S_YY_study2, ...
                                   N1, N2,                   ...
                                   1,                        ...
                                   SNP_id );
                               
% OUTPUT 
%   A vector with three entries:
%   1) SNP id,
%   2) leading canonical correlation value,
%   3) -log10(p-value).                    
              



%% Multi-SNP analysis

% Provide SNP ids in the following cell array format:
  SNP_ids = {'rs10', 'rs80', 'rs140', 'rs170', 'rs172'};
% SNP_ids = {'10', '80', '140', '170', '172'};


metaCCA_result3     =     metaCCA( 2,                        ...
                                   S_XY_study1, S_XY_study2, ...
                                   0, 0,                     ...
                                   S_YY_study1, S_YY_study2, ...
                                   N1, N2,                   ...
                                   2,                        ...
                                   SNP_ids,                  ...
                                   S_XX_study1, S_XX_study2 );
               
metaCCAplus_result3 = metaCCAplus( 2,                        ...
                                   S_XY_study1, S_XY_study2, ...
                                   0, 0,                     ...
                                   S_YY_study1, S_YY_study2, ...
                                   N1, N2,                   ...
                                   2,                        ...
                                   SNP_ids,                  ...
                                   S_XX_study1, S_XX_study2 );
                        
               
% OUTPUT 
%   A vector with three entries:
%   1) list of SNP ids,
%   2) leading canonical correlation value,
%   3) -log10(p-value).               
