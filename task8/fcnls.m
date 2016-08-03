function f = fcnls(x)
   global N V p w medidas
   
   Rr=x(1);
   Xr=x(2);
   Xs=x(3);
   Rs=x(4);
   Xm=x(5);

   s = medidas(:,1);
   y=zeros(size(s,2),2);
   for i=1:length(s)
     A=Rs*(1+Xr/Xm)+(1+Xs/Xm)*Rr/s(i);
     B=Xr+Xs*(1+Xr/Xm)-Rs*Rr/Xm/s(i);
     C=1+Xr/Xm;
     D=Rr/Xm/s(i);
     y(i,1) = V*sqrt((C*C+D*D)/(A*A+B*B));  %Ic
     y(i,2) = 3*V*V*(p/w)*(Rr/s(i))/(A*A+B*B); %Tc
   end
   
   f(1)=sum((y(:,1)-medidas(:,2)).^2);
   f(2)=sum((y(:,2)-medidas(:,3)).^2);
end
