function [sMap, datanorm, data] = training(T, munits)
%{ 
Creation and training of the SOM sMap trained on the T dataset with munits
number of neurons. Returns the SOM sMap, the variance normalized dataset
datanorm and the non-normalized dataset data.
-------
Input
-------
T: table
munits: int
-------
Output
------
sMap: som struct from the SOM Toolbox
datanorm: data struct from the SOM Toolbox
data: data struct from the SOM Toolbox
%}

%Creation of the data struct from the table T
x = T{:,:};
data=som_data_struct(x,'label_names',T.Properties.VariableNames,'comp_names',T.Properties.VariableNames);
%Variance normalization
datanorm=som_normalize(data,'var');
%Creation of the map sMap
sMap=som_make(datanorm,'munits',munits,'tracking',0,'training','Long');

%Batch training of the map
sMap=som_batchtrain(sMap,datanorm,'radius',[100 1],'trainlen',80,'tracking',0); 
sMap=som_batchtrain(sMap,datanorm,'radius',[1 0.5],'trainlen',250,'tracking',0); 
sMap=som_batchtrain(sMap,datanorm,'radius',[0.5 0.1],'trainlen',200,'tracking',0); 
sMap=som_batchtrain(sMap,datanorm,'radius',[0.1 0.001],'trainlen',100,'tracking',0); 

end