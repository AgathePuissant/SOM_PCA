function [] = umat_hitmap(sMap, datanorm)
%{
Plot function to plot simultaneously the U-Mat and the hitmap 
from the SOM sMap and the normalized data associated with the SOM.
------
Input
------
sMap: som map struct from the SOM Toolbox
datanorm : som data struct from the SOM Toolbox
%}

%U-Matrix
figure; som_show(sMap,'umat',1:size(sMap.codebook,2));%pour visualiser la matrice des distance (UMatrice)
saveas(gcf,'UMAT.eps','epsc')
close all

%Hitmap
figure; som_show(sMap,'empty','Hitmap');
hits = som_hits(sMap,datanorm);
som_show_add('hit',hits,'Marker','lattice','MarkerColor','m','Text','on','TextColor','k');
saveas(gcf,'Hitmap.eps','epsc')
close all

end