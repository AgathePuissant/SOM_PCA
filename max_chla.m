function[] = max_chla(sMap, datanorm,depths,pigment_names,PCA,loadings,scale_parameter,center_parameter,axes_per_pigment,pigment_depths_names,sat_names)
%{
Plot function to plot the DCM concentration and the DCM depth for the
neurons of the SOM sMap. If the SOM was trained on PCA pretreated data, the
PCA parameter should be changed to true and the PCA information should be
provided to allow the reconstruction of the profiles for the neurons of the
map.
------
Input
------
sMap: som map struct from the SOM Toolbox
datanorm : som data struct from the SOM Toolbox
depths: vector
pigment_names: string vector
PCA: boolean
loadings_pig: table
scale_parameter: table
center_parameter: table
axes_per_pigment: int
pigment_depths_names: string vector
sat_names:string vector
%}

    hits = som_hits(sMap,datanorm); %To indicate the empty neurons
    
    if (~PCA==true)
        Ss = som_denormalize(sMap) ;
    else
        %Reconstruction of the neurons profiles from the PCA information
    
        %Denormalization of the SOM codebook to correctly reconstruct the
        %profiles
        sMap = som_denormalize(sMap) ;
        
        codebook_interm=[];
        liste=axes_per_pigment;
        nb=cumsum(axes_per_pigment);
        for i=1:length(pigment_names)     
        codebook_interm(:,(i-1)*9+1:(i-1)*9+9) = (sMap.codebook(:,nb(i)-liste(i)+1:nb(i))*loadings{:,nb(i)-liste(i)+1:nb(i)}'.*repmat(scale_parameter{:,i}',size(sMap.codebook,1),1))+repmat(center_parameter{:,i}',size(sMap.codebook,1),1);
        end
        depths = ["5.0","8.35","13.92","23.23","38.75","64.63","107.81","179.85","300.0"];
        codebook = codebook_interm;
        sMap.codebook = [codebook,sMap.codebook(:,sum(axes_per_pigment)+1:end)];
        sMap.comp_names=[pigment_depths_names,sat_names];
        Ss=sMap;
    end

    %Get the maximum value of the first 9 variables corresponding to the
    %chla variables and its index corresponding to its depth
    [M,I]= max(Ss.codebook(:,1:length(depths)),[],2);
    
    %Fill the empty neurons with an outside value
    I(hits==0) = length(depths)+1;
    
    som_cplane(Ss,I);
    
    title('DCM depth')
    shading flat
    C=colorbar;
    C.Label.String ='m';
    C.FontSize=8;
    C.TickLabels = [depths,"Empty"] ;
    C.Direction = 'reverse';
    colormap([jet(9);1,1,1])
    
    saveas(gcf,'max_chla_depth.eps','epsc')
    
    close all
    
    %Fill the empty neurons with an outside value
    M(hits == 0) = max(M)+1 ;
    
    %Plotting the log values to have a better colorbar spreading
    M = log10(M);
    
    som_cplane(Ss,M);
    
    title('DCM concentration')
    shading flat
    C=colorbar;
    C.Label.String ='mg.m^-^3';
    C.FontSize=8;
    C.TickLabels=round(10.^C.Ticks,2);
    colormap([jet;1,1,1])
    
    saveas(gcf,'max_chla_con.eps','epsc')
    
    
    close all
    
end
