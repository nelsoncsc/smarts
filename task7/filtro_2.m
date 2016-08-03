clear all
close all
clc

rng(1);
%
an=1;
wn=2*pi;
wa=2^(6)*wn;
h=wn/wa;

%
x1=0;
x2=1.25*wn;

x=[x1;x2];

%
N=2^(9)-1;
ks=zeros(N+1,1);
sm=zeros(N+1,1);
sp=zeros(N+1,1);
se=zeros(N+1,1);
we=zeros(N+1,1);

%
I2=eye(2,2);
P=1e6*I2;
sigv=1e-2;
R=sigv^2;
sigw=1e-2;

%
for i=0:N,
    ks(i+1)=i;
    %
    s=an*sin(h*wn*i);
    sp(i+1) = s;
    
    w = x2;
    wp(i+1) = w;
   
    ye = an*sin(x1);
    yep(i+1) = ye;
   
    ym = s + sigv*randn();
    ymp(i+1) = ym;

    x1=x(1);
    x2=x(2);

    F=[1,  h
       0,  1];
    
    % Q
    Q=[0,     0
       0,  sigw^2];
    
   % C
    C=[an*cos(x1), 0];
    
    %
    P=F*P*F'+Q;
    
    K=P*C'*inv(C*P*C'+R);
    
    % f(x)
    x1=x1+h*x2;
    x2=x2+sigw*randn();
    
    %
    x=[x1;x2];

    x=x+K*(ym-ye);
    P=P-K*C*P;
end
%
figure(1)
subplot(2,1,1)
plot(ks',wp,'b',ks',wn,'r','Linewidth',2)
set(gca,'FontSize',20)
xlabel('k')
grid
axis([0 length(ks) 0 10])
legend('wp','wn')
%
subplot(2,1,2)
plot(ks',sp,'b',ks',ymp,'k',ks',yep,'r','Linewidth',2)
set(gca,'FontSize',20)
xlabel('k')
grid
axis([0 length(ks) -1.5 1.5])
legend('s','ymp','yep')