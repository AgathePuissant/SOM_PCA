function[] = component_maps(sMap,nb_pigments,unit,pigment_names,depths,PCA,loadings,scale_parameter,center_parameter,axes_per_pigment,pigment_depths_names,sat_names)
%{
Plot function to plot the value of each component for each referent vector
of the SOM sMap. If the SOM was performed on PCA pretreated data, the PCA 
parameter should be changed to true and the PCA information should be proemptyd 
to allow the reconstruction of the profiles for the referent vectors.
------
Input
------
sMap: som map struct from the SOM Toolbox
nb_pigments: int
unit: string
pigment_names: string vector
depths: vector
PCA: boolean
loadings: table
scale_parameter: table
center_parameter: table
axes_per_pigment: int
pigment_depths_names: string vector
sat_names:string vector
%}

if (~PCA==true)
    Ss = som_denormalize(sMap) ;
else
    %Reconstruction of the neurons profiles from the PCA information
    
    %Denormalization of the SOM codebook to correctly reconstruct the
    %profiles
    sMap = som_denormalize(sMap) ;
    
    codebook_interm=[];
    nbaxes_cumsum=cumsum(axes_per_pigment);
    for i=1:length(pigment_names)     
    codebook_interm(:,(i-1)*9+1:(i-1)*9+9) = (sMap.codebook(:,nbaxes_cumsum(i)-axes_per_pigment(i)+1:nbaxes_cumsum(i))*loadings{:,nbaxes_cumsum(i)-axes_per_pigment(i)+1:nbaxes_cumsum(i)}'.*repmat(scale_parameter{:,i}',size(sMap.codebook,1),1))+repmat(center_parameter{:,i}',size(sMap.codebook,1),1);
    end
    depths = ["5.0","8.35","13.92","23.23","38.75","64.63","107.81","179.85","300.0"];
    codebook = codebook_interm;
    sMap.codebook = [codebook,sMap.codebook(:,sum(axes_per_pigment)+1:end)];
    sMap.comp_names=[pigment_depths_names,sat_names];
    Ss=sMap;
end


for i = 0:length(depths):length(depths)*nb_pigments-1
    figure;
    [ha,~]=tight_subplot(mod(length(depths),3)+fix(length(depths)/3),3,[.1 .1],[.1 .1],[.1 .1]);
    title_map = strrep(Ss.comp_names(i+1:i+length(depths)),'_','\_');
    for j = 1:length(depths)
        axes(ha(j));
        som_cplane(Ss,Ss.codebook(:,i+j));
        title((title_map(j)),'FontSize',8)
        shading flat
        caxis([min(Ss.codebook(:,i+1:i+length(depths)),[],'all') max(Ss.codebook(:,i+1:i+length(depths)),[],'all')])
        C=colorbar;
        C.Label.String=unit;
        C.FontSize=8;
        colormap(jet)
    end
    
    for empty = length(depths)+1:length(ha)
        subplot(ha(empty)); 
        set(gca,'Color','none','XColor','none','YColor','none');
    end
    str = sprintf('image_%d.eps',i);
    saveas(gcf,str,'epsc')
    close all
end


figure;
[ha,~]=tight_subplot(4,3,[.1 .1],[.1 .1],[.1 .1]);
title_map = strrep(Ss.comp_names(nb_pigments*length(depths)+1:end),'_','\_');
for i=1:size(sMap.codebook,2)-nb_pigments*length(depths)
    axes(ha(i));
    som_cplane(Ss,Ss.codebook(:,nb_pigments*length(depths)+i));
    title((title_map(i)),'FontSize',8)
    shading flat

    C=colorbar;
    colormap(jet)
end

for empty = size(sMap.codebook,2)-nb_pigments*length(depths)+1:length(ha)
    subplot(ha(empty)); 
    set(gca,'Color','none','XColor','none','YColor','none');
end

str = 'image_varsat.eps';
saveas(gcf,str,'epsc')
close all

end