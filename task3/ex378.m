%% Task 3 - Part II: Examples 7 and 8
% Nelson Campos
%%

%% Example 3.7: Direct self-tuner with d0 = 1
% Consider the system in Example 3.1. Since deg A = 2 and deg B = 1,
% we have deg Am = 2 and deg Ao = 0. Hence Ao = 1, and we will choose
% Bm = qAm(1). Equation (3.29) in Algorithm 3.3 then gives T = qAm(1). The
% controller structure is given by deg R = deg S = deg T = deg A ? 1 = 1. The
% model given by Eq. (3.27) therefore becomes
%
% $$ y(t) = r_0u_f(t-1)+r_1u_f(t-2)+s_0y_f(t-1)+s_1y_f(t-2) $$
%
% where:
%
% $$  u _f(t)+a_{m1}u f(t-1)+a_{m1}u_f(t-2)=u(t) $$
%
% $$ y_f(t)+a_{m1}y_f(t-1)+a_{m1}y_f(t-2) = y(t) $$
%
% Figure 1 shows the simulation of the results. 
% In a practical case the time delay and the order of the process that we
% would like to control are not known. It is therefore natural to consider these
% variables as design parameters that are chosen by the user. The parameter d0
% is of particular importance for a direct algorithm. In the next example we show
% that “ringing” can be avoided simply by increasing the value of d0.
%%

close all
clear all

a1 = -1.6065;
a2 = 0.6065;
b0 = 0.1065;
b1 = 0.0902;
%---------------
bm0 = 0.1761;
am1 = -1.3205;
am2 = 0.4966;

t_step = 1/10;
t = [0:t_step:100-t_step];
uc = zeros(1,size(t,2));
y = zeros(1, size(t,2));
T = 50; 

for i=1:size(t,2)/T
    for j=1:T
        index = (i-1)*T+j;
        if j <= T/2
            uc(index) = 1;
        else uc(index) = -1;
        end
    end
end

for i=3:size(t,2)
  y(i) = -am1*y(i-1)-am2*y(i-2)+bm0*uc(i-1);
end

% R(q) = B+ = q+b1/b0, S(q) = s0q+s1, T(q) = AoBm' = bm0q/b0
s0 = (am1-a1)/b0;
s1 = (am2-a2)/b0;
R = [1 b1/b0];
S = [s0 s1];
T = [bm0/b0 0];

u = zeros(1, size(t,2));

for i=2:size(t,2)
  u(i) = -R(2)*u(i-1)+T(1)*uc(i)-S(1)*y(i)-S(2)*y(i-1);
end

figure(1), subplot(2,1,1), plot(t,uc, t, y), xlabel('Time'), ylabel('Amplitude'), title('y(t) vs t and uc(t) vs t')
           subplot(2,1,2), plot(t,u), xlabel('Time'), ylabel('Amplitude'), title('u(t) vs t with zero cancellation')

           
uf = zeros(1,size(t,2));
yf = zeros(1, size(t,2));

for i=3:size(t,2)
  uf(i) = -am1*uf(i-1)-am2*uf(i-2)+u(i);
  yf(i) = -am1*yf(i-1)-am2*yf(i-2)+y(i);
end

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
est_r0 = zeros(1,size(t,2));
est_r1 = zeros(1,size(t,2));
est_s0 = zeros(1,size(t,2));
est_s1 = zeros(1,size(t,2));

for i = 3 : size(t,2)
          
phi(:,i) = [uf(i-1) uf(i-2) yf(i-1) yf(i-2)]';
K = P * phi(:,i)* inv(lambda + phi(:,i)'* P * phi(:,i));
P = ((I - K * phi(:,i)')* P)/lambda;
est_y = phi(:,i)'*theta(:,i-1);
error(i) = y(i)-est_y;

theta(:,i) = theta(:,i-1) + K * error(i);
est_r0(i) = theta(1,i);
est_r1(i) = theta(2,i);
est_s0(i) = theta(3,i);
est_s1(i) = theta(4,i);

end

t0 = 1+am1+am2;
r0 = est_r0(size(t,2));
r1 = est_r1(size(t,2));
s0 = est_s0(size(t,2));
s1 = est_s1(size(t,2));

for i=3:size(t,2)
  y(i) = r0*uf(i-1)+r1*uf(i-2)+s0*yf(i-1)+s1*yf(i-2);
  u(i) = (-r1*u(i-1)+t0*uc(i)-s0*y(i)-s1*y(i-1))/r0;
end

figure(2), subplot(2,1,1), plot(t,uc, t, y), xlabel('Time'), ylabel('Amplitude'), title('y(t) vs t and uc(t) vs t direct self-tuning d=1')
           subplot(2,1,2), plot(t,u), xlabel('Time'), ylabel('Amplitude'), title('u(t) vs t with zero cancellation direct self-tuning d=1')

figure(3), plot(t,est_r1./est_r0, 'b--', t, est_s0./est_r0, 'r', t, est_s1./est_r0, 'g', t, t0./est_r0, 'k'),xlabel('Time'), ylabel('Amplitude'),
title('Estimated parameters r1/r0, s0/r0, s1/s0, t0/s0 and d=1'), legend('r1/r0', 's0/r0', 's1/r0', 't0/r0')

%% Example 3.8: Direct self-tuner with d0 = 2
% In the derivation of the direct algorithm the parameter d0 was the pole excess
% of the plant. Assume for a moment that we do not know the value of d0 and that
% we treat it as a design parameter instead. Figure 2 shows the simulation of 
% the results, but with d0 = 2 instead of d0 = 1.
% We thus find the interesting and surprising result that cancellation of the
% process zero can be avoided by increasing the parameter d0. This observation
% will be explained later when we will be analyzing the algorithms.

%%

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
est_r0 = zeros(1,size(t,2));
est_r1 = zeros(1,size(t,2));
est_s0 = zeros(1,size(t,2));
est_s1 = zeros(1,size(t,2));

%d=2
for i = 4 : size(t,2)
          
phi(:,i) = [uf(i-2) uf(i-3) yf(i-2) yf(i-3)]';
K = P * phi(:,i)* inv(lambda + phi(:,i)'* P * phi(:,i));
P = ((I - K * phi(:,i)')* P)/lambda;
est_y = phi(:,i)'*theta(:,i-1);
error(i) = y(i)-est_y;

theta(:,i) = theta(:,i-1) + K * error(i);
est_r0(i) = theta(1,i);
est_r1(i) = theta(2,i);
est_s0(i) = theta(3,i);
est_s1(i) = theta(4,i);

end

t0 = 1+am1+am2;
r0 = est_r0(size(t,2));
r1 = est_r1(size(t,2));
s0 = est_s0(size(t,2));
s1 = est_s1(size(t,2));

for i=4:size(t,2)
  y(i) = r0*uf(i-2)+r1*uf(i-3)+s0*yf(i-2)+s1*yf(i-3);
  u(i) = (-r1*u(i-1)+t0*uc(i)-s0*y(i)-s1*y(i-1))/r0;
end

figure(4), subplot(2,1,1), plot(t,uc, t, y), xlabel('Time'), ylabel('Amplitude'), title('y(t) vs t and uc(t) vs t direct self-tuning d=2')
           subplot(2,1,2), plot(t,u), xlabel('Time'), ylabel('Amplitude'), title('u(t) vs t with zero cancellation direct self-tuning d=2')

figure(5), plot(t,est_r1./est_r0, 'b--', t, est_s0./est_r0, 'r', t, est_s1./est_r0, 'g', t, t0./est_r0, 'k'),xlabel('Time'), ylabel('Amplitude'),
title('Estimated parameters r1/r0, s0/r0, s1/s0, t0/s0 and d=2'), legend('r1/r0', 's0/r0', 's1/r0', 't0/r0')
