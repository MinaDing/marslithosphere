function [adm,cor,sd_adm,sd_cor]=admcor(la,Hilm,faailm,degrees)
Hpower = degrees*0;
faapower = degrees*0;
crosspower = degrees*0;
for i=1:length(degrees)
    s = find(la==degrees(i));
    crosspower(i) = sum(Hilm(s).*faailm(s));
    Hpower(i) = sum(Hilm(s).*Hilm(s));
    faapower(i) = sum(faailm(s).*faailm(s));
    
end
adm= crosspower./Hpower;
cor = crosspower./sqrt(Hpower)./sqrt(faapower);
sd_adm = adm./abs(cor).*sqrt((1-cor.^2)/2/degrees);
sd_cor = (1-cor.^2)./sqrt(2*degrees);

end
%admerr = sqrt(faapower./Hpower.*(1-cor.^2)./degrees/2);