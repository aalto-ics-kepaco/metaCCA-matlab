function S_YY = estimate_Syy(S_XY)

% This function uses Pearson correlation to compute phenotypic 
% correlation matrix S_YY based on univariate summary statistics S_XY.


% Anna Cichonska
% anna.cichonska@helsinki.fi




% Validating if the file format is correct
if isstruct(S_XY)~= 1 
    error('Input data are not in the correct format; make sure that the header line is given.');
end


% Betas - trait ids
trait_ids_betas_temp = S_XY.textdata(1, 4:2:end);
for i = 1:length(trait_ids_betas_temp)
    trait_ids_betas(i) = getfield( strsplit(trait_ids_betas_temp{i},'_'), {1});
end

% SE - trait ids
trait_ids_se_temp    = S_XY.textdata(1, 5:2:end);
for i = 1:length(trait_ids_se_temp)
    trait_ids_se(i) = getfield( strsplit(trait_ids_se_temp{i},'_'), {1});
end


% Validating if trait ids of the corresponding regression coefficients 
% and standard errors match
if isequal(trait_ids_betas, trait_ids_se) == 0
    error('Trait ids of regression coefficients and standard errors do not match');
end


% Validating if trait ids are unique
if length(unique(trait_ids_betas)) ~= length(trait_ids_betas)
   error('Trait ids are not unique!');
end



Betas = S_XY.data(:, 1:2:end);

S_YY = corr(Betas);

S_YY = [trait_ids_betas', num2cell(S_YY)];
