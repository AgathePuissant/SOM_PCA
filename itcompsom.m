function[sMap,completed,datanorm,liste_corr] = itcompsom(T,n_it,munits,save)
%{ 
Iterative Completion SOM: train a SOM sMap on the dataset T with a number munits 
of neurons using the ITCOMPSOM algorithm (ref) with a number n_it of iterations. 
Returns the trained SOM sMap, the completed data, the normalized training data datanorm
and the matrix of squared correlation between the variables of the table T.
The completed data can be saved by changing the save parameter to true.
-------
Input
-------
T: table
n_it: int
munits: int
save: boolean
-------
Output
------
sMap: som struct from the SOM Toolbox
completed: table
datanorm: data struct from the SOM Toolbox
liste_corr: matrix
%}

%Find the location of the NaN values
stock_nan = isnan(T{:,:}) ;

%Sort the dataset by ascending order of number of missinf values
nb_nan = sum(stock_nan,2) ;
[~,index] = sort(nb_nan) ;

%Creating and normalizing the data struct and computing the squared
%correlation matrix
S1 = som_data_struct(T{:,:},'comp_names',T.Properties.VariableNames);
S1n = som_normalize(S1,'var') ;
liste_corr = (corr(S1n.data,'Type','Spearman','rows','complete')).^2;

for turn = 1:n_it
    
    %Selection of a subset of the dataset
    subset = index(1:ceil(length(index)/(n_it+1-turn)));
    
    %Creating the data struct from this subset and normalizing it according
    %to the normalization parameters for the entire dataset
    Sobs = som_data_struct(T{subset,:},'comp_names',T.Properties.VariableNames);
    Sobsn = som_normalize(Sobs,S1n) ;
    
    %Creating and training the SOM on the subset
    nneurons = round(munits*(turn/n_it));
    sMapobs=som_make(Sobsn,'munits',nneurons,'training','long','tracking',0);
    sMapobs=som_batchtrain(sMapobs,Sobsn,'radius',[100 1],'trainlen',80,'tracking',0);
    sMapobs=som_batchtrain(sMapobs,Sobsn,'radius',[1 0.5],'trainlen',250,'tracking',0);
    sMapobs=som_batchtrain(sMapobs,Sobsn,'radius',[0.5 0.1],'trainlen',200,'tracking',0);
    sMapobs=som_batchtrain(sMapobs,Sobsn,'radius',[0.1 0.001],'trainlen',100,'tracking',0);
    
    %SOM denormalization to complete with denormalized values
    sMapobs_denorm=som_denormalize(sMapobs);
    
    %Create a vector containing the correlation weight between NaN
    %variables and the others
    for row = 1:size(S1n.data,1)
        nan_row = find(isnan(S1n.data(row,:))) ;
        corr_row = ones(1,size(S1n.data,2));
        for nan_index = 1:length(nan_row)
            nan_corr = liste_corr(nan_row(nan_index),:) ;
            corr_row = corr_row + nan_corr ;
        end
        
        %Finding the BMUs with the weighted distance
        corr_matrix = repmat(corr_row,size(sMapobs.codebook(:,:),1),1);
        bmusobs(row)=som_bmus(corr_matrix.*sMapobs.codebook(:,:),corr_row.*S1n.data(row,:));
    end
   
    %Completing the missing observations with the corresponding components from
    %the BMU's referent vector
    for i=1:size(T,1)
        if (~isnan(bmusobs(i)))
           T(i,find(stock_nan(i,:)>0))=array2table(sMapobs_denorm.codebook(bmusobs(i),find(stock_nan(i,:)>0)));
        end
    end
    
end

completed = T ;

if (save==true)
        writetable(completed,'completed.csv','Delimiter',';');
end

datanorm = S1n ;
sMap = sMapobs ;

end