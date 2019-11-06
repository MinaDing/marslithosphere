function Mij = MatrixM(lmax,Lwin,Sww)

% lmax = 90;
% Lwin = 17;
% Sww = [0:Lwin]*0+1;

% the localization matrix defined in [Wieczorek and Simons, 2005, localized
% spectal analysis on the sphere]
% input parameters:
%   lmax: the maximum spherical harmonic degree for global spectra Sfg (Sfg
%   is the global spectra with spherical harmonic degree 0-lmax, therefore
%   length is lmax+1
%   Lwin: the maximum spherical harmonic degree for the localization window
%   Sww: power spectra of the localization window with degrees from 0 to
%   Lwin, therefore the length of Sww must be Lwin+1

addpath('/Users/dingmin/Desktop/MarsCodes_Results/Polar&plain/NorthPole1_finished/Subroutines','-end')

if length(Sww)~=Lwin+1
    error('The length of Shh does not match the given maximum degree Lwin')
end


[iii,jjj,kkk]=ndgrid(0:lmax-Lwin,0:lmax,0:Lwin);
N = numel(iii);
[Ni,Nj,Nk]=size(iii);
ivec = reshape(iii,N,1);
jvec = reshape(jjj,N,1);
kvec = reshape(kkk,N,1);

C2vec = 0*jvec;
for i=1:N
    C2vec(i) = ClebschGordan(jvec(i),kvec(i),ivec(i),0,0,0)^2;
end
C2matrix = reshape(C2vec,Ni*Nj,Nk);
%imatrix = reshape(ivec,Ni*Nj,Nk);
%jmatrix = reshape(jvec,Ni*Nj,Nk);
Sww = reshape(Sww,Nk,1);
Mvec = C2matrix*Sww;
Mij = reshape(Mvec,Ni,Nj);
%iout = reshape(imatrix(:,1),Ni,Nj); 
%jout = reshape(jmatrix(:,1),Ni,Nj); 
%should be consistent with [iout,jout]=ndgrid(0:lmax-Lwin,0:lmax);
% I found that ndgrid is better than meshgrid olz

end



