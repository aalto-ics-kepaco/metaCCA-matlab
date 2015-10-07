function result = metaCCAplus(varargin)

% This function performs metaCCA+.
% 'varargin' argument is a cell array that contains all the inputs, 
% where each input is in its own cell.

% Anna Cichonska
% anna.cichonska@helsinki.fi




nr_in  = nargin;            % number of inputs
M      = varargin{1};       % number of studies to be analysed 
def_in = 4*M + 1;           % number of inputs expected in a default single-SNP analysis

% Validating if any optional inputs were given 
if nr_in == def_in
    option = 0; 
else 
    option = varargin{def_in+1}; 
    if option ~= 1 &&  option ~= 2
        error('Wrong indicator of the analysis type. Please provide >>1<< for a single-SNP analysis of selected SNP, >>2<< for a multi-SNP analysis.');  
    end
end




for i = 1:M    % M - #studies
    % Phenotypic correlaton matrix (1st column contains trait ids)
    S_YY_temp       = varargin{1+ 2*M +i};
    % A matrix without trait ids
    S_YY{i}         = cell2mat(S_YY_temp(:, 2:end));
    trait_id_syy{i} = S_YY_temp(:,1)';
    
    % Information if univariate analysis has been performed on standardized
    % or non-standardized data 
    uni_info(i)     = varargin{1+ M +i};
    
    % Number of individuals
    N{i}            = varargin{1+ 3*M +i};
end






% DEFAULT SINGLE-SNP ANALYSIS

if option == 0      
    
    for i = 1:M
        
        % Summary statistics
        S_XY_temp = varargin{1+i};
        
        % Validating if trait ids are provided as strings
        if size(S_XY_temp.textdata,2)-3 ~= size(S_XY_temp.data,2)
            error( strcat('Trait ids in the header line should be given as strings (study ', num2str(i), ').') );
        end
        
        trait_ids_betas{i} = S_XY_temp.textdata(1,4:2:end);
        trait_ids_se{i}    = S_XY_temp.textdata(1,5:2:end);
        
        % Validating if trait ids of the correspondung regression coefficients 
        % and standard errors match
        if isequal(trait_ids_betas{i}, trait_ids_se{i}) == 0
            error( strcat('Trait ids of regression coefficients and standard errors do not match (study ', num2str(i), ').') );
        end
        
        SNPid{i}   =  S_XY_temp.textdata(2:end,1);
        allele0{i} =  S_XY_temp.textdata(2:end,2);
        
        S_XY{i}    =  S_XY_temp.data(:, 1:2:end);   % univariate regression coefficients 
        se{i}      =  S_XY_temp.data(:, 2:2:end);
    end


    % Validating if SNP ids match between different studies
    if isequal(SNPid{:}) ~= 1
       error('SNP ids in summary statistics of different studies do not match.');  
    end
    % Validating if trait ids in S_XY match between different studies
    if isequal(trait_ids_betas{:}) ~= 1
       error('Trait ids in summary statistics of different studies do not match.');  
    end
    % Validating if trait ids in S_YY match between different studies
    if isequal(trait_id_syy{:}) ~= 1
       error('Trait ids in phenotypic correlation structures of different studies do not match.');  
    end
    % Validating if trait ids in S_YY and S_XY match 
    if isequal(trait_id_syy{1}, trait_ids_betas{1}) ~= 1
       error('Trait ids in phenotypic correlation structures and summary statistics do not match.');  
    end
    % Validating if allele_0 match between different studies
    if isequal(allele0{:}) ~= 1
       error('Alleles 0 in summary statistics of different studies do not match.');  
    end
    
    % Validating if SNP ids and trait ids are unique
    if length(unique(SNPid{1})) ~= length(SNPid{1})
        error('SNP ids are not unique! ');
    end
    if length(unique(trait_ids_betas{1})) ~= length(trait_ids_betas{1})
        error('Trait ids are not unique! ');
    end
       
    
    % Normalizing regression coefficients (if the univariate analysis has
    % been performed on non-standardised data)
    for i = 1:M
        if uni_info(i) == 0
            S_XY_norm{i} = normalize_Sxy(S_XY{i}, se{i}, N{i});
        elseif uni_info(i) == 1
            S_XY_norm{i} = S_XY{i};
        else 
            error( strcat('Wrong indicator of the univariate analysis type (study ', num2str(i), '). Please provide >>0<< or >>1<<.') );
        end
    end
        
        
    % Pooling covariance matrices of the same type
    [~, C_YY, C_XY, N_tot] = pool_cov_cell( repmat({1}, 1, M),  S_YY,  S_XY_norm,  N);
    
    
    result = cell(size(C_XY,1),3);
    result(:,1) = SNPid{1};
  
    % Analysing one SNP at a time (against all given traits)
    for i_snp = 1:size(C_XY,1)
       
        % Bulding a full covariance matrix
        full_cov  =  [1, C_XY(i_snp,:); C_XY(i_snp,:)', C_YY];
        
        % Ensuring the PSD property of the full covariance matrix
        [~, C_XY_out, C_YY_out] = shrink_plus(full_cov, 1);
        
        % Canonical Correlation Analysis (CCA)
        % genotype-phenotype association result
        [r, ~, ~, ~, ~, ~, pval] = my_cca(1, C_YY_out, C_XY_out, N_tot); 
        result{i_snp,2} = r;
        result{i_snp,3} = -log10(pval);
    end
    
    
    
    
    
  
% SINGLE-SNP ANALYSIS OF ONE SELECTED SNP    

elseif option == 1  
    
    if nr_in ~= def_in+2
        error('Wrong number of inputs. ')
    end
    
    selected_SNPid = varargin{end};
       
    
    % All SNP ids
    for i = 1:M
        SNPid{i}   = varargin{1+i}.textdata(2:end,1);
    end

    
    % Summary statistics for a given SNP
    for i = 1:M
        
        S_XY_temp = varargin{1+i};
        
        % Validating if trait ids are provided as strings
        if size(S_XY_temp.textdata,2)-3 ~= size(S_XY_temp.data,2)
            error( strcat('Trait ids in the header line should be given as strings (study ', num2str(i), ').') );
        end
                      
        trait_ids_betas{i} = S_XY_temp.textdata(1,4:2:end);
        trait_ids_se{i}    = S_XY_temp.textdata(1,5:2:end);
        
        % Validating if trait ids of the correspondung regression coefficients 
        % and standard errors match
        if isequal(trait_ids_betas{i}, trait_ids_se{i}) == 0
            error( strcat('Trait ids of regression coefficients and standard errors do not match (study ', num2str(i), ').') );
        end

        % match SNP id
        h_id        = find( strcmp(selected_SNPid, SNPid{i}) == 1 );
        if length(h_id) > 1
            error( strcat('There is more than one SNP with id "', selected_SNPid, '" in study ', num2str(i), '.') );
        end
        S_XY{i}     = S_XY_temp.data(h_id, 1:2:end);   % univariate regression coefficients 
        se{i}       = S_XY_temp.data(h_id, 2:2:end);     
        allele0{i}  = S_XY_temp.textdata(h_id+1 ,2);   % h_id+1 --> first entry is a header 'allele_0'
    end
    
    
    % Validating if trait ids in S_XY match between different studies
    if isequal(trait_ids_betas{:}) ~= 1
       error('Trait ids in summary statistics of different studies do not match.');  
    end
    % Validating if trait ids in S_YY match between different studies
    if isequal(trait_id_syy{:}) ~= 1
       error('Trait ids in phenotypic correlation structures of different studies do not match.');  
    end
    % Validating if trait ids in S_YY and S_XY match 
    if isequal(trait_id_syy{1}, trait_ids_betas{1}) ~= 1
       error('Trait ids in phenotypic correlation structures and summary statistics do not match.');  
    end
    % Validating if allele_0 match between different studies
    if isequal(allele0{:}) ~= 1
       error('Alleles 0 in summary statistics of different studies do not match.');  
    end
    
    % Validating if trait ids are unique
    if length(unique(trait_ids_betas{1})) ~= length(trait_ids_betas{1})
        error('Trait ids are not unique! ');
    end
    
    
    % Normalizing regression coefficients (if the univariate analysis has
    % been performed on non-standardised data)
    for i = 1:M
        if uni_info(i) == 0
            S_XY_norm{i} = normalize_Sxy(S_XY{i}, se{i}, N{i});
        elseif uni_info(i) == 1
            S_XY_norm{i} = S_XY{i};
        else 
            error( strcat('Wrong indicator of the univariate analysis type (study ', num2str(i), '). Please provide >>0<< or >>1<<.') );
        end
    end

    
    % Pooling covariance matrices of the same type
    [~, C_YY, C_XY, N_tot] = pool_cov_cell( repmat({1}, 1, M),  S_YY,  S_XY_norm,  N);
    
    result = cell(1,3);
    result{1,1} = selected_SNPid;
    
    
    % Bulding a full covariance matrix
    full_cov  =  [1, C_XY; C_XY', C_YY];

    
    % Ensuring the PSD property of the full covariance matrix
    [~, C_XY_out, C_YY_out] = shrink_plus(full_cov, 1);

    
    % Canonical Correlation Analysis (CCA)
    % genotype-phenotype association result
    [r, ~, ~, ~, ~, ~, pval] = my_cca(1, C_YY_out, C_XY_out, N_tot); 
    result{1,2} = r;
    result{1,3} = -log10(pval);
    
    
    
    
    
    
    
% MULTI-SNP ANALYSIS   

elseif option == 2  
    
    if nr_in ~= def_in+2+M
        error('Wrong number of inputs. ')
    end
    
    selected_SNPid = varargin{def_in+2}; 
    
    % All SNP ids
    for i = 1:M
        SNPid{i}   = varargin{1+i}.textdata(2:end,1);
        
        % ids in S_XX
        S_XX_temp  = varargin{def_in+2+i};
        if  strcmp( class(S_XX_temp), 'double') == 1      % then, SNP ids are given as numbers
            SNPid_inSxx{i} = strread( num2str(S_XX_temp(:,1)'), '%s');
            S_XX_all{i}    = S_XX_temp(:, 2:end);
        else                                              % SNP ids are given as text
            SNPid_inSxx{i} = S_XX_temp.textdata;
            S_XX_all{i}    = S_XX_temp.data;
        end 
        
    end
    

    
    % Validating if SNP ids in S_XY match between different studies
    if isequal(SNPid{:}) ~= 1
       error('SNP ids in summary statistics of different studies do not match.');  
    end
    % Validating if SNP ids in S_XX match between different studies
    if isequal(SNPid_inSxx{:}) ~= 1
       error('SNP ids in S_XX of different studies do not match.');  
    end
    % Validating if SNP ids in S_XX match those in S_XY
    if isequal(SNPid{1}, SNPid_inSxx{1}) ~= 1
       error('SNP ids in S_XX and S_XY do not match. SNP ids can be given as numbers or strings; please ensure that the type is the same in S_XX and S_XY.');  
    end
    % If all above is ture, then SNP ids in S_XX and S_XY are the same and in the same
    % order.
    
    
    
   
    % Summary statistics and genotypic correlation structures correspinding to given SNPs
    for i = 1:M
        
        S_XY_temp = varargin{1+i};
        
        % Validating if trait ids are provided as strings
        if size(S_XY_temp.textdata,2)-3 ~= size(S_XY_temp.data,2)
            error( strcat('Trait ids in the header line should be given as strings (study ', num2str(i), ').') );
        end

        trait_ids_betas{i} = S_XY_temp.textdata(1,4:2:end);
        trait_ids_se{i}    = S_XY_temp.textdata(1,5:2:end);
        
        % Validating if trait ids of the correspondung regression coefficients 
        % and standard errors match
        if isequal(trait_ids_betas{i}, trait_ids_se{i}) == 0
            error( strcat('Trait ids of regression coefficients and standard errors do not match (study ', num2str(i), ').') );
        end
        
          
        % match SNP ids
        h_id = find(ismember(SNPid{i}, selected_SNPid));
        if length(h_id) > length(selected_SNPid)
            error( strcat('SNP ids are not unique (study ', num2str(i), ').') );
        elseif length(h_id) < length(selected_SNPid)
            error( 'At least one of the given SNPs not found in the data');
        end
  
        S_XY{i}     = S_XY_temp.data(h_id, 1:2:end);       % univariate regression coefficients 
        se{i}       = S_XY_temp.data(h_id, 2:2:end); 
        allele0{i}  = S_XY_temp.textdata(h_id+1 ,2);       % h_id+1 --> first entry is a header 'allele_0'
   
        S_XX{i}     = S_XX_all{i}(h_id, h_id);  
    end
    
    
     % Validating if trait ids in S_XY match between different studies
    if isequal(trait_ids_betas{:}) ~= 1
       error('Trait ids in summary statistics of different studies do not match.');  
    end
    % Validating if trait ids in S_YY match between different studies
    if isequal(trait_id_syy{:}) ~= 1
       error('Trait ids in phenotypic correlation structures of different studies do not match.');  
    end
    % Validating if trait ids in S_YY and S_XY match 
    if isequal(trait_id_syy{1}, trait_ids_betas{1}) ~= 1
       error('Trait ids in phenotypic correlation structures and summary statistics do not match.');  
    end
    % Validating if allele_0 match between different studies
    if isequal(allele0{:}) ~= 1
       error('Alleles 0 in summary statistics of different studies do not match.');  
    end
    
    % Validating if trait ids are unique
    if length(unique(trait_ids_betas{1})) ~= length(trait_ids_betas{1})
        error('Trait ids are not unique! ');
    end
    
    
    % Normalizing regression coefficients (if the univariate analysis has
    % been performed on non-standardised data)
    for i = 1:M
        if uni_info(i) == 0
            S_XY_norm{i} = normalize_Sxy(S_XY{i}, se{i}, N{i});
        elseif uni_info(i) == 1
            S_XY_norm{i} = S_XY{i};
        else 
            error( strcat('Wrong indicator of the univariate analysis type (study ', num2str(i), '). Please provide >>0<< or >>1<<.') );
        end
    end
    
    
    
    % Pooling covariance matrices of the same type
    [C_XX, C_YY, C_XY, N_tot] = pool_cov_cell( S_XX,  S_YY,  S_XY_norm,  N);
    
   
    result{1,1} = selected_SNPid;
    
    
    % Bulding a full covariance matrix
    full_cov  =  [C_XX, C_XY; C_XY', C_YY];

    
    % Ensuring the PSD property of the full covariance matrix
    [C_XX_out, C_XY_out, C_YY_out] = shrink_plus(full_cov, size(C_XX,1) );

    
    % Canonical Correlation Analysis (CCA)
    % genotype-phenotype association result
    [r, ~, ~, ~, ~, ~, pval] = my_cca(C_XX_out, C_YY_out, C_XY_out, N_tot); 
    result{1,2} = r(1);
    result{1,3} = -log10(pval(1));
   
end