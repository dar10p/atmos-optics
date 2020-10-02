function[uout,L2]=propFF(uin,L1,lambda,z);
% propagation - Fraunhofer pattern
% assumes uniform sampling
% uin - source plane field
% L1 - source plane side length
% lambda - wavelength
% z - propagation distance
% L2 - observation plane side length
% uout - observation plane field
%
[M,N]=size(uin); %get input field array size
dx1=L1/M; %source sample interval
k=2*pi/lambda; %wavenumber
%
L2=lambda*z/dx1; %obs sidelength
dx2=lambda*z/L1; %obs sample interval
x2=-L2/2:dx2:L2/2-dx2; %obs coords
[X2,Y2]=meshgrid(x2,x2);
%
c=1/(j*lambda*z)*exp(j*k/(2*z)*(X2.^2+Y2.^2));
uout=c.*ifftshift(fft2(fftshift(uin)))*dx1^2;
end