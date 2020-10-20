function[X_temp,Y_temp]=reconstruction(observed,estimated,depths,pigment_names,nbaxes_cumsum,axes_per_pigment,loadings,center_parameter,scale_parameter)
Y_temp=[];
X_temp=[];

for i=1:length(pigment_names)
Y_temp(:,(i-1)*length(depths)+1:(i-1)*length(depths)+length(depths)) = (observed(:,nbaxes_cumsum(i)-axes_per_pigment(i)+1:nbaxes_cumsum(i))*loadings{:,nbaxes_cumsum(i)-axes_per_pigment(i)+1:nbaxes_cumsum(i)}'.*repmat(scale_parameter{:,i}',size(observed,1),1))+repmat(center_parameter{:,i}',size(observed,1),1);
X_temp(:,(i-1)*length(depths)+1:(i-1)*length(depths)+length(depths)) = (estimated(:,nbaxes_cumsum(i)-axes_per_pigment(i)+1:nbaxes_cumsum(i))*loadings{:,nbaxes_cumsum(i)-axes_per_pigment(i)+1:nbaxes_cumsum(i)}'.*repmat(scale_parameter{:,i}',size(estimated,1),1))+repmat(center_parameter{:,i}',size(estimated,1),1);
end

end