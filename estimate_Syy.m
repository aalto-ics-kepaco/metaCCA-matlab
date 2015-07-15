function S_YY = estimate_Syy(S_XY)

% This function uses Pearson correlation to compute
% phenotypic correlation matrix based on summary statistics 

% Anna Cichonska
% anna.cichonska@helsinki.fi

S_YY = corr(S_XY);




 