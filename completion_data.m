function[completed, BMUS] = completion_data(sMap, T, save)
%{ 
Completion of the T dataset by finding the corresponding BMUs from the SOM
sMap. Returns the completed datasat and the list of corresponding BMUS. The output can be saved to a csv file by changing the save parameter
to true.
-------
Input
-------
sMap: som struct from SOM Toolbox
T: table
save: boolean
-------
Output
------
completed: table
BMUS: vector
%}

    % Normalization of the data we want to complete according to the SOM
    % normalization parameters
    data=som_data_struct(T{:,:});
    data_norm=som_normalize(data,sMap);
    
    %Find the BMU for each observation of the data
    BMUS = som_bmus(sMap,data_norm);
    
    %SOM denormalization to complete with denormalized values
    sMap_denorm = som_denormalize(sMap) ;
    
    %Find the location of the NaN values
    stock_nan = isnan(T{:,:}) ;
    
    for row=1:size(BMUS)
        
        %Check if a BMU could be found and if there are NaN values in the
        %corresponding observation
        if (~isnan(BMUS(row)) && ~isempty(find(stock_nan(row,:)>0)))
            
            T{row,find(stock_nan(row,:)>0)} = sMap_denorm.codebook(BMUS(row),find(stock_nan(row,:)>0)) ;
        
        end
    end
    
    completed = T;

    if (save==true)
        writetable(completed,'completed.csv','Delimiter',';');
    end

end