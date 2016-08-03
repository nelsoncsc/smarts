%clear all
%close all

N = 1000;
t = 1:N;
u = zeros(N,1);
std_e = 0.5;
mean_e = 0;
e = sqrt(std_e)*randn(N,1)+mean_e;
a = -0.8;
b = 0.5;
c = -0.5;
y = zeros(N,1);

T = 200;
u = zeros(N,1);
y = zeros(N,1);

for i=1:N/T
  for j=1:T
      index = (i-1)*T+j;
      if j <= T/2
          u(index) = 1;
      else u(index) = 0;
      end
  end
end

for index=2:N
    y(index) = -a*y(index-1)+b*u(index-1)+e(index)+c*e(index-1);
end

P = 100*eye(3);
phi = zeros(3,1);
theta = zeros(3,1);
error_pred = zeros(N,1);
est_a = zeros(N,1);
est_b = zeros(N,1);
est_c = zeros(N,1);

for i=2:N
  phi = [-y(i-1) u(i-1) error_pred(i-1)]';
  error_pred(i) = y(i)-(phi')*theta;  
  K = P*phi/(1+phi'*P*phi);
  theta = theta + K *(y(i) - (phi')* theta); 
  P = (eye(3) - K * (phi'))* P; 
  est_a(i) = theta(1);
  est_b(i) = theta(2);
  est_c(i) = theta(3);
end

figure(4), plot(t,a,t, est_a, t,b,t,est_b,t,c,t,est_c), legend('real a', 'estimated a', 'real b', 'estimated b','real c', 'estimated c'), xlabel('time'), ylabel('Estimated parameters'),
title('Example 2.13 (b)')
