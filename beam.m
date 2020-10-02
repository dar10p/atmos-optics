% display the field
%
L1=0.5; %side length
M=250; %number of samples
dx1=L1/M; %src sample interval
x1=-L1/2:dx1:L1/2-dx1; %src coords
y1=x1;
lambda=0.5*10^-6; %wavelength
k=2*pi/lambda; %wavenumber
w=0.051; %source half width (m) 
%w=0.011 % for  Fraunhofer
w1=1e-3; %for circular case
z=2000; %propagation dist (m) % 50 for circulr case Fraunhofer


[X1,Y1]=meshgrid(x1,y1);
uin=(rectangularPulse(-1/2,1/2,(X1/(2*w)))).*(rectangularPulse(-1/2,1/2,(Y1/(2*w)))); %src field
%uin=circshift(sqrt(X1.^2+Y1.^2)/w1.*exp(i.*2*pi*w1),2)
Iin=abs(uin.^2); %src irradiance
%
figure(1)
imagesc(x1,y1,Iin);
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title('z= 0 m');

%%
%apply tilt 
deg=pi/180;
alpha=5.0e-5; %rad
theta=65*deg;
[uin]=tilt(uin,L1,lambda,alpha,theta);

%%
%apply focus 
zf=2000;
[uin]=focus(uin,L1,lambda,zf);

%%
% Propagator
%uout=propTF(uin,L1,lambda,z); % Transference function
uout=propIR(uin,L1,lambda,z); % Fresnel Impulse response
x2=x1; %obs coords
y2=y1;

Iout=abs(uout.^2); %obs irrad 
%%
% propagation - Fraunhofer pattern
[uout,L2]=propFF(uin,L1,lambda,z);
dx2=L2/M;
x2=-L2/2:dx2:L2/2-dx2; %oBs coords
y2=x2;
I2=abs(uout.^2); %obs irrad
imagesc(x2,y2,nthroot(I2,3));%stretch image contrast

%%
figure(2) %display obs irrad
imagesc(x2,y2,Iout);
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title(['z= ',num2str(z),' m']);
%%
figure(3) %irradiance profile
plot(x2,Iout(M/2+1,:));
xlabel('x (m)'); ylabel('Irradiance');
title(['z= ',num2str(z),' m']);
%
figure(4) %plot obs field mag
plot(x2,abs(uout(M/2+1,:)));
xlabel('x (m)'); ylabel('Magnitude');
title(['z= ',num2str(z),' m']);
%
figure(5) %plot obs field phase
plot(x2,unwrap(angle(uout(M/2+1,:))));
xlabel('x (m)'); ylabel('Phase (rad)');
title(['z= ',num2str(z),' m']);