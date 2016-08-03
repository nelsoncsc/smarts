%
clear all
close all
%
global N V p w medidas
%
Rr=3.84;
Xr=6.789;
Xs=1.658;
Rs=1.93;
Xm=38.7;
%
N=25;
V=220;
f=60;
w=2*pi*f;
p=2;
In=5.8/50;
Pn=1.5e3/50;
Tn=8.0/50;
s=[0.02:(1-0.02)/N:1.0]';
for i=1:length(s),
%
A=Rs*(1+Xr/Xm)+(1+Xs/Xm)*Rr/s(i);
B=Xr+Xs*(1+Xr/Xm)-Rs*Rr/Xm/s(i);
C=1+Xr/Xm;
D=Rr/Xm/s(i);
Im(i,1)=V*sqrt((C*C+D*D)/(A*A+B*B))+In*randn();
Tm(i,1)=3*V*V*(p/w)*(Rr/s(i))/(A*A+B*B)+Tn*randn();
end
%
medidas=[s Im Tm];
%
lb = [0  0  0 0 0]';
ub = [10 10 5 5 50]';
%
FitnessFunction = @fcnls;
numberOfVariables = 5;
options = gaoptimset('PopulationSize',200, 'PlotFcn',@gaplotpareto);
x = gamultiobj(FitnessFunction,numberOfVariables,[],[],[],[],lb,ub,options);
msg = ['ideal - Rr = ' num2str(Rr) ', Xr = ' num2str(Xr) ', Xs = ' num2str(Xs) ', Rs = ' num2str(Rs) ', Xm = ' num2str(Xm)];
disp(msg)
%
linha=size(x,2);
Rr=x(linha,1);
Xr=x(linha,2);
Xs=x(linha,3);
Rs=x(linha,4);
Xm=x(linha,5);
%
msg = ['final - Rr = ' num2str(Rr) ', Xr = ' num2str(Xr) ', Xs = ' num2str(Xs) ', Rs = ' num2str(Rs) ', Xm = ' num2str(Xm)];
disp(msg)
%
for i=1:length(s)
A=Rs*(1+Xr/Xm)+(1+Xs/Xm)*Rr/s(i);
B=Xr+Xs*(1+Xr/Xm)-Rs*Rr/Xm/s(i);
C=1+Xr/Xm;
D=Rr/Xm/s(i);

%
Ie(i,1)=V*sqrt((C*C+D*D)/(A*A+B*B));
Te(i,1)=3*V*V*(p/w)*(Rr/s(i))/(A*A+B*B);
end
%
figure(1)
subplot(2,1,1)
plot(s,Im,'ro',s,Ie,'b-'), grid, title('I \times s'),title('Ixs for an Induction Machine'), ylabel('Ixs'), xlabel('slip'), legend('Estimated', 'Measured')
subplot(2,1,2)
plot(s,Tm,'ro',s,Te,'b-'), grid, title('T \times s'), title('Txs for an Induction Machine'), ylabel('Txs'), xlabel('slip'), legend('Estimated', 'Measured')
