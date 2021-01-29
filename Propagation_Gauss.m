% Propagador de Fourier
clear all; close all;
set(0,'defaultTextInterpreter','latex');

N = 2^9;                % Numero de puntos, escogemos una potencia de 2
                        % para que el calculo de la transformada sea más
                        % rápido

% Creamos el vector de indices
NV = (-N/2:1:N/2-1);     % Con esto nos aseguramos de tener la frecuencia 0
                         % y de que NV tiene N elementos
L = 1e-3;                % Las unidades son metros
dx = 2*L/N;
kmax = pi/dx;
[X,Y] = meshgrid(NV*dx);
%r = sqrt(X.^2+Y.^2);
[Kx,Ky] = meshgrid(kmax*2/N*NV);

lambda = 633e-9;
k = 2*pi/lambda;
kt = sqrt(Kx.^2+Ky.^2);
for n=1:4
Nf = 1;

if n==1
    Nf = 10;
elseif n==2
    Nf = 5;
elseif n==3
    Nf = 1;
elseif n==4
    Nf = 0.1;
end
D = N/10;
z = (L/10).^2/(4*lambda*Nf);      %Distancia de propagación maxima
nz = 300;
dz = z/nz;

% Propagador 
Prop = (exp(-1i*0.5*dz*(kt.^2)/k));      %Full paraxial propagator
%Prop = exp(i*dz*sqrt(k^2-kt.^2));       %Non paraxial propagator

% Perfil Inicial 
Tlens = 1;        %Función de transmitacia

U0=zeros(N);
U0(floor(N/2-D/2):floor(N/2+D/2),floor(N/2-D/2):floor(N/2+D/2))=1;
U0 = U0.*Tlens;
Ur(:,1) = U0(:,N/2+1);

% Inicia la propagacion 
F = fftshift(fft2(U0));
for ii=1:nz
    F = F.*Prop;
    A = ifft2(F);
    Ur(:,ii+1)=A(:,N/2+1);
end

%  Grafica del campo inicial 
figure(1),subplot(2,2,n),surf(X(N/2+1,:)/L,Y(:,N/2+1)/L,abs(A));shading interp,lighting phong, colormap hot, view(2)
ejes = gca;
ejes.FontSize = 13;
if n==1
    title('$Nf = 10$')
elseif n==2
    title('$Nf = 5$')
elseif n==3
    title('$Nf = 1$')
elseif n==4
    title('$Nf = 0.1$')
end
xlabel('$x$/L','FontSize',20);
ylabel('$y$/L','FontSize',20);
axis square;

end

