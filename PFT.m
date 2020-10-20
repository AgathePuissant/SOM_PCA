function[micro,diatoms,dino,prok,prochloro,hapto]=PFT(chla,dvchla,X19hex,fucox,perid,zeax,coeff)
    
DPW = reshape(repmat(coeff,[1,9])',[1,9*5]).*[dvchla,X19hex,fucox,perid,zeax];
% SDPW = sum(DPW,2);
SDPW = chla;

micro = DPW(:,19:27)+DPW(:,28:36)./SDPW;
diatoms = DPW(:,19:27)./SDPW;
dino = DPW(:,28:36)./SDPW;
prok = DPW(:,37:45)./SDPW;
prochloro = DPW(:,1:9)./SDPW;
hapto = (Xn(chla).*DPW(:,10:18))./SDPW;

end