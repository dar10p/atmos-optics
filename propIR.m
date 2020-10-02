function[uout]=propIR(uin,L,lambda,z);
% propagation - impulse response approach
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

x=-L/2:dx:L/2-dx; %spatial coords
[X,Y]=meshgrid(x,x);

h=1/(j*lambda*z)*exp(j*k/(2*z)*(X.^2+Y.^2)); %impulse
H=fft2(fftshift(h))*dx^2; %create trans func
Uin=fft2(fftshift(uin)); %shift, fft src field
Uout=H.*Uin; %multiply
uout=ifftshift(ifft2(Uout)); %inv fft, center obs field
end