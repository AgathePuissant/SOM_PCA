function[completed,BMUS] =  completion_data_itcompsom(sMap,T,liste_corr,save)
%{ 
Completion of the T dataset by finding the corresponding BMUs from the SOM sMap using the between variable squared correlation (liste_corr) weighted euclidean distance. Returns the completed datasat and the list of corresponding BMUS. The output can be saved to a csv file by changing the save parameter to true.
-------
Input
-------
sMap: som struct from SOM Toolbox
T: table
liste_corr: matrix
save: boolean
-------
Output
------
completed: table
BMUS: vector
%}

% Normalization of the data we want to complete according to the SOM
% normalization parameters
data = som_data_struct(T{:,:},'comp_names',T.Properties.VariableNames) ;
data_norm = som_normalize(data,sMap) ;

%Find the location of the NaN values
stock_nan = isnan(T{:,:}) ;

%SOM denormalization to complete with denormalized values
sMap_denorm = som_denormalize(sMap);

for row = 1:size(data_norm.data,1)
    
    %Create a vector containing the correlation weight between NaN
    %variables and the others
    nan_row = find(isnan(data_norm.data(row,:))) ;
    corr_row = ones(1,size(data_norm.data,2));
    for nan_index = 1:length(nan_row)
        nan_corr = liste_corr(nan_row(nan_index),:) ;
        corr_row = corr_row + nan_corr ;
    end

    %Finding the BMUs with the weighted distance
    corr_matrix = repmat(corr_row,size(sMap.codebook(:,:),1),1);
    BMUS(row)=som_bmus(corr_matrix.*sMap.codebook(:,:),corr_row.*data_norm.data(row,:));

end 

%Completing the missing observations with the corresponding components from
%the BMU's referent vector
for i=1:size(T,1)
    if (~isnan(BMUS(i)))
        T(i,find(stock_nan(i,:)>0))=array2table(sMap_denorm.codebook(BMUS(i),find(stock_nan(i,:)>0)));
    end
end

completed = T ;

if (save==true)
    writeTrix(completed{:,:},'completed.csv','Delimiter',';');
end

end