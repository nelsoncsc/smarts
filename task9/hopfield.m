clear all
close all
clc

f = [0.30 0.32];
a = [1.0 0.8];
phi = [pi/4 pi/8];

n = [0: 1: 10000];

x1 = a(1)*cos(2*pi*f(1)*n+phi(1));
x2 = a(2)*cos(2*pi*f(2)*n+phi(2));
xn = x1+x2;

  
% figure(1), subplot (2,2,1), plot(n,x1),
%            subplot (2,2,2), plot(n,x2),
%            subplot (2,2,3), plot(n,xn);



ak = 1.25*a;
fk = 1.25*f;
phik = 1.5*phi;

y1 = ak(1)*cos(2*pi*fk(1)*n+phik(1));
y2 = ak(2)*cos(2*pi*fk(2)*n+phik(2));
yn=y1+y2;

dn=xn-yn;
for i=1:length(n)-1 
    
  y1 = ak(1)*cos(2*pi*fk(1)*n+phik(1));
  y2 = ak(2)*cos(2*pi*fk(2)*n+phik(2));
  yn=y1+y2;
  dn = xn-yn; 
  deltaE = -2*pi*sum( (fk(1)*ak(1)*sin(2*pi*fk(1)+phik(1))+fk(2)*ak(2)*sin(2*pi*fk(2)+phik(2))).*dn );
  
  fk(1) = fk(1)-deltaE/sum( 2*pi.*n*ak(1).*sin(2*pi*fk(1)*n+phik(1)).*dn);
  fk(2) = fk(2)-deltaE/sum( 2*pi.*n*ak(2).*sin(2*pi*fk(2)*n+phik(2)).*dn);
  
  ak(1) = ak(1)+deltaE/sum( cos(2*pi*fk(1)*n+phik(1)).*dn);
  ak(2) = ak(2)+deltaE/sum( cos(2*pi*fk(2)*n+phik(2)).*dn);
  
  phik(1) = phik(1)-deltaE/sum( ak(1)*sin(2*pi*fk(1)*n+phik(1)).*dn );
  phik(2) = phik(2)-deltaE/sum( ak(2)*sin(2*pi*fk(2)*n+phik(2)).*dn );
end

a, ak
f, fk
phi, phik