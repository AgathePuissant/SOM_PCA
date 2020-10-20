function[res]=Xn(chla)
    vec1=linspace(0.001,0.08);
    vec2=linspace(0,1);
    mdl = fitlm(vec2,vec1);
    coeff=mdl.Coefficients.Estimate;
    
    res=coeff(2)*chla+coeff(1); 
    res(chla<0.001)=0;
    res(chla>0.08)=1;
    
end