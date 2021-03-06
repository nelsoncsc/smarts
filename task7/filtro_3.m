clear all
close all
rng(1);
%
an=1;
wn=2*pi;
wa=2^(6)*wn;
h=wn/wa;
%
x1=0;
x2=1.25*an;
x3=1.25*wn;
x=[x1;x2;x3];
%
N=2^(9)-1;
ks=zeros(N+1,1);
sm=zeros(N+1,1);
sp=zeros(N+1,1);
se=zeros(N+1,1);
we=zeros(N+1,1);
%
I3=eye(3,3);
P=1e6*I3;
sigv=1e-2;
R=sigv^2;
sigw=1e-2;
%
for i=0:N,
    ks(i+1)=i;
    %
    s=an*sin(h*wn*i);
    sp(i+1)=s;
    ym=s+sigv*randn();
    sm(i+1)=ym;
    %
    x1=x(1);x2=x(2);x3=x(3);
    ye=x3*x1;
    se(i+1)=ye;
    we(i+1)=x3;
    % F
    F=[1,       h,        0
       -h*x3^2, 1,    -2*h*x1*x3
       0,       0,        1];
    % Q
    Q=[0,0,0
       0,0,0
       0,0,sigw^2];
    % C
    %C=[x3,0,x1];
     %C=[1,0,0];
     C=[x3,0,x1];
    %
    P=F*P*F'+Q;
    K=P*C'*inv(C*P*C'+R);
    % f(x)
    x1=x1+h*x2;
    x2=x2-h*x1*x3^2;
    x3=x3+sigw*randn();
    %
    x=[x1;x2;x3];
    x=x+K*(ym-ye);
    P=P-K*C*P;
end
%
figure(1)
subplot(2,1,1)
stairs([ks ks ks],[sp sm se])
set(gca,'FontSize',20)
set(gca,'defaultLineLineWidth',2)
xlabel('k')
grid
legend('s(k)','y_{m}(k)','y_{e}(k)')
%
subplot(2,1,2)
stairs([ks ks],[wn*ones(size(we)) we])
set(gca,'FontSize',20)
set(gca,'defaultLineLineWidth',2)
xlabel('k')
grid
legend('\omega_{n}(k)','\omega_{e}(k)')