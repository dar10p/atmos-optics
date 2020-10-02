function[uout]=propTF(uin,L,lambda,z);
% propagation - transfer function approach
% assumes same x and y side lengths and
% uniform sampling
% uin - source plane field
% L - source and observation plane side length
% lambda - wavelength
% z - propagation distance
% uout - observation plane field

[M,N]=size(uin); %get input field array size
dx=L/M; %sample interval
k=2*pi/lambda; %wavenumber

fx=-1/(2*dx):1/L:1/(2*dx)-1/L; %freq coords
[FX,FY]=meshgrid(fx,fx);

H=exp(-j*pi*lambda*z*(FX.^2+FY.^2)); %trans func
H=fftshift(H); %shift trans func
Uin=fft2(fftshift(uin)); %shift, fft src field
Uout=H.*Uin; %multiply
uout=ifftshift(ifft2(Uout)); %inv fft, center obs field