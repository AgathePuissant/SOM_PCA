function[sMap, estimated, observed, BMUS, idx_obs]=validation_k_fold(T,k,nb_pigment_variables,select_varsat,munits)
%{
Function performing the cross validation k fold of the SOM_PCA method with munits number of neurons on
the data T. The satellite data used for the inversion can be selected with
the select_varsat parameter which lists the variables tha t will be used
for the inversion.
------
Input
------
T: table
k: int
nb_pigment_variables: int
select_varsat: string vector
munits: int
------
Output
------
sMap: som map struct from the SOM Toolbox
estimated: matrix
observed: matrix
BMUS: vector
idx_obs: vector

%}
    
    %Initialize a list of the observation index to choose from during the
    %folds
    idx = [1:size(T,1)];
    
    %Initializing the variables to be filled during the folds
    estimated = [] ;
    observed = [] ;
    BMUS = [] ;
    idx_obs = []; % List of the index chosen at each fold to keep track of the observations
    
    for i=1:k
        sprintf("Fold %d",i)
        
        %Selecting a subset of the observations
        k_idx = randsample(idx,fix(size(T,1)/k));
        %Removing them from the pool to choose from
        idx_obs=[idx_obs,k_idx];
        idx(ismember(idx,k_idx)) = [] ;

        validation_set = T(k_idx,:); % the chosen observations are used as a validation set
        training_set = T(~ismember([1:size(T,1)],k_idx),:) ; % While the rest is used as a training setS

        %Training of the map
        [sMap, datanorm, ~] = training(training_set,munits);

        %Keeping the real observed values and masking the pigment variables
        %values
        real = validation_set{:,1:nb_pigment_variables};
        validation_set{:,~ismember(datanorm.comp_names,select_varsat)} = NaN;

        %Completing the masked values from the selected satellite variables with the SOM
        [completed,bmus] = completion_data(sMap,validation_set,false) ;
        completed = completed(:,1:nb_pigment_variables);
        completed = completed{:,:};
        

        estimated = [estimated ; completed] ;
        observed = [observed ; real] ;
        BMUS = [BMUS ; bmus] ;      
        
    end
    
    %If some observation remains, do a final round to use them
    if (fix(size(T,1)/k)~=0)
        k_idx = idx;

        validation_set = T(k_idx,:);
        training_set = T(~ismember([1:size(T,1)],k_idx),:) ;

        [sMap, datanorm, ~] = training(training_set,munits);

        real = validation_set{:,1:nb_pigment_variables};
        validation_set{:,~ismember(datanorm.comp_names,select_varsat)} = NaN;

        [completed,bmus] = completion_data(sMap,validation_set,false) ;
        completed = completed(:,1:nb_pigment_variables);
        completed = completed{:,:};
        

        estimated = [estimated ; completed] ;
        observed = [observed ; real] ;
        BMUS = [BMUS ; bmus] ; 

        
    end
    
end