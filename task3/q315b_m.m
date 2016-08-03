%% Question 3.15 part II
%%

close all
clear all

Ki = 1.5;
K0 = 1;
deltaK = -0.5;
K = K0+deltaK;
a0 = 1.4;
delta_a = 2;
a = 1;
a = a0+delta_a;

sim('q315b.mdl')

uc = uc.signals.values;
u = u.signals.values;
y = y.signals.values;
t = 1:size(uc,1);

% Estimation goes here
N = 4;
theta = zeros(N, size(t,2));
phi = zeros(N, size(t,2));
error = zeros(1,size(t,2));
K = zeros(N,1);
P = zeros(N);

%The initial conditions
%theta = [a1 a2 b0 b1]
theta(1,1) = 0;
theta(2,1) = 0;
theta(3,1) = 0.01;
theta(4,1) = 0.2;
I = eye(N);
P(1,1) = 100;
P(2,2) = 100;
P(3,3) = 1;
P(4,4) = 1;

lambda =  1;
est_y = zeros(1, size(t,2));

for i = 3 : size(t,2)
          
phi(:,i) = [-y(i-1) -y(i-2) u(i-1) u(i-2)]';
K = P * phi(:,i)* inv(lambda + phi(:,i)'* P * phi(:,i));
P = ((I - K * phi(:,i)')* P)/lambda;
est_y(i) = phi(:,i)'*theta(:,i-1);
error(i) = y(i)-est_y(i);

theta(:,i) = theta(:,i-1) + K * error(i);

end

figure(1), subplot(2,2,1), plot(uc), xlabel('Time'), ylabel('Amplitude'), title('input u_c(t): step')
subplot(2,2,2), plot(u), xlabel('Time'), ylabel('Amplitude'), title('signal of the controller u(t)')
subplot(2,2,3),plot(y), xlabel('Time'), ylabel('Amplitude'), title('y(t) with noise and step input')
subplot(2,2,4),plot(est_y), xlabel('Time'), ylabel('Amplitude'), title('Estimation of y(t)')

%% References
% [1] Adaptive Control 2nd Edition Astrom & Wittenmark  
%%
