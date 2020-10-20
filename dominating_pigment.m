function[] = dominating_pigment(sMap,datanorm,depths,pigment_names,PCA,loadings,scale_parameter,center_parameter,axes_per_pigment,pigment_depths_names,sat_names)
%{
Plot function to plot the dominating pigment for the
neurons of the SOM sMap. If the SOM was trained on PCA pretreated data, the
PCA parameter should be changed to true and the PCA information should be
provided to allow the reconstruction of the diles for the neurons of the
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
    %Reconstruction of the neurons diles from the PCA information
    
    %Denormalization of the SOM codebook to correctly reconstruct the
    %diles
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

pigment_names(strcmp(pigment_names,'total chla'))=[]; %Removing chla

figure;
[ha,~]=tight_subplot(mod(length(depths),3)+fix(length(depths)/3),3,[.1 .1],[.1 .15],[.1 .1]);

for d = 1:length(depths)
    
    %Getting the index of the max value for each depth, which correspond to the dominating
    %pigment
    [~,I] = max(Ss.codebook(:,d+length(depths):length(depths):length(depths)*length(pigment_names)+length(depths)),[],2) ;
    
    %Fill the empty neurons with an outside value
    I(hits==0) = length(pigment_names)+1;
    
    axes(ha(d));    
    
    som_cplane(Ss,I);
    
    title(sprintf('%.2f m',depths(d)))
    shading flat
    C=colorbar;
    C.Label.String ='Pigment';
    C.FontSize=8;
    caxis([1,length(pigment_names)+1])
    C.Ticks = [1:length(pigment_names)+1];
    C.TickLabels = [pigment_names,'Vide'] ;
    C.Direction = 'reverse';
    colormap([jet(length(pigment_names));1,1,1])
    
end

sgtitle('Dominating pigment','FontSize',12)
saveas(gcf,'max_pigment.eps','epsc')
close all

%Decomment the following part to get the details of each pigment value for
%each depth and each neuron

% for d = 1:length(depths)
%     figure;
%     pies = Ss.codebook(:,d+length(depths):length(depths):length(depths)*length(pigment_names)+length(depths));
%     pies(pies<0)=0;
%     pies(:,size(pies,2)+1) = 0 ;
%     pies(hits==0,:)=0;
%     pies(hits==0,end)=1; 
%     
%     som_pieplane(Ss,pies,[jet(length(pigment_names));1,1,1]);
%     title(sprintf('%.2f m',depths(d)))
%     shading flat
%     C=colorbar;
%     C.Label.String ='Pigment';
%     C.FontSize=8;
%     caxis([1,length(pigment_names)+1])
%     C.Ticks = [1:length(pigment_names)+1];
%     C.TickLabels = [pigment_names,'Empty neuron'] ;
%     C.Direction = 'reverse';
%     colormap([jet(length(pigment_names));1,1,1])
%    
%     str = sprintf('dpigpies_%d.eps',d);
%     saveas(gcf,str,'epsc')
%     close all
% end

end