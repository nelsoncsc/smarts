%% Task 3 - Part III: Example 3.9
% Nelson Campos
%%

%%
% Example 3.9: Effect of Load Disturbances
% Consider the system in Example 3.5, that is, an indirect self-tuning regulator
% with no zero cancellation. We will now make a simulation that is identical to
% the one shown in Fig. 3.6 except that the load disturbance will be v(t) = 0.5 for
% t $$ \geq $$ 40. A forgetting factor $$ \lambda $$ = 0.98 has also been introduced; otherwise, the
% conditions are identical to those in Example 3.5. The behavior of the system
% is shown in Fig. 3. A load disturbance may be disastrous. Notice that the response is strongly asymmetric.
% The reason for this is that the controller parameters change rapidly when the
% control signal changes. Rapid changes of the estimates in response to command signals indicates that
% the model structure is not correct. The parameter estimates also change significantly at the step in the load disturbance. When the command signal is
% constant, the parameters appear to settle at constant values that are far from
% the true parameters.
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

t_step = 1;
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
           subplot(2,1,2), plot(t,u),, xlabel('Time'), ylabel('Amplitude'), title('u(t) vs t with zero cancellation')
           
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

for i = 3 : size(t,2)
          
phi(:,i) = [-y(i-1) -y(i-2) u(i-1) u(i-2)]';
K = P * phi(:,i)* inv(lambda + phi(:,i)'* P * phi(:,i));
P = ((I - K * phi(:,i)')* P)/lambda;
est_y = phi(:,i)'*theta(:,i-1);
error(i) = y(i)-est_y;

theta(:,i) = theta(:,i-1) + K * error(i);

end

figure(2), subplot(2,1,1), plot(t,a1,t,theta(1,:),t,a2,t,theta(2,:)), xlabel('Time'), ylabel('Amplitude'), title('Estimated a1 and a2 with zero cancellation'), legend('a1', 'est a1', 'a2', 'est a2')
           subplot(2,1,2), plot(t,b0,t,theta(3,:),t,b1,t,theta(4,:)), xlabel('Time'), ylabel('Amplitude'), title('Estimated b0 and b1 with zero cancellation'), legend('b0', 'est b0', 'b1', 'est b1')

%%Example 3.5
%%

% Hm(q)= beta*B(q)/Am(q) = Bm(q)/Am(q), bm0 = beta*b0
% beta = (1+am1+am2)/(b0+b1)

a0 = 0; % deg(A(q)) = 1
beta = (1+am1+am2)/(b0+b1);
bm0 = beta*b0;
bm1 = beta*b1;

for i=3:size(t,2)
  y(i) = -am1*y(i-1)-am2*y(i-2)+bm0*uc(i-1)+bm1*uc(i-2);
end

y_noise = y;

for i=40:size(t,2)
  y_noise(i) = y(i)+0.5;
end

r1 = b1/b0 + (b1^2-am1*b0*b1+am2*b0^2)*(-b1+a0*b0)/(b0*(b1^2-a1*b0*b1+a2*b0^2));
s0_2 = b1*(a0*am1-a2-am1*a1+a1^2+am2-a1*a0)/(b1^2-a1*b0*b1+a2*b0^2)+...
    b0*(am1*a2-a1*a2-a0*am2+a0*a2)/(b1^2-a1*b0*b1+a2*b0^2);

s1_2 = b1*(a1*a2-am1*a2+a0*am2-a0*a2)/(b1^2-a1*b0*b1+a2*b0^2)+...
     b0*(a2*am2-a2^2-a0*am2*a1+a0*a2*am1)/(b1^2-a1*b0*b1+a2*b0^2);

R2 = [1 r1];
S2 = [s0_2 s1_2];
T2 = [bm0/b0];

for i=2:size(t,2)
  u(i) = -R2(2)*u(i-1)+T2(1)*uc(i)-S2(1)*y(i)-S2(2)*y(i-1);
end

figure(3), subplot(2,1,1), plot(t,uc, t, y_noise), xlabel('Time'), ylabel('Amplitude'), title('y_noise(t) vs t and uc(t) vs t')
           subplot(2,1,2), plot(t,u),, xlabel('Time'), ylabel('Amplitude'), title('u(t) vs t without zero cancellation')

%%%
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

lambda =  0.98;

for i = 3 : size(t,2)
          
phi(:,i) = [-y_noise(i-1) -y_noise(i-2) u(i-1) u(i-2)]';
K = P * phi(:,i)* inv(lambda + phi(:,i)'* P * phi(:,i));
P = ((I - K * phi(:,i)')* P)/lambda;
est_y = phi(:,i)'*theta(:,i-1);
error(i) = y(i)-est_y;

theta(:,i) = theta(:,i-1) + K * error(i);

end

figure(4), subplot(2,1,1), plot(t,a1,t,theta(1,:),t,a2,t,theta(2,:)), xlabel('Time'), ylabel('Amplitude'), title('Estimated a1 and a2 without zero cancellation'), legend('a1', 'est a1', 'a2', 'est a2')
           subplot(2,1,2), plot(t,b0,t,theta(3,:),t,b1,t,theta(4,:)), xlabel('Time'), ylabel('Amplitude'), title('Estimated b0 and b1 without zero cancellation'), legend('b0', 'est b0', 'b1', 'est b1')

% The controller obtained from the estimated parameters
est_A = [theta(1,size(t,2)) theta(2,size(t,2))];
est_B = [theta(3,size(t,2)) theta(4,size(t,2))];
h = 0.5; % The sampling time

est_Hq = tf(est_B, est_A, h, 'variable', 'q');
est_beta  = (1+am1+am2)/(est_B(1)+est_B(2));

b0 = est_B(1);
b1 = est_B(2);

est_r1 = b1/b0 + (b1^2-am1*b0*b1+am2*b0^2)*(-b1+a0*b0)/(b0*(b1^2-a1*b0*b1+a2*b0^2));
est_s0_2 = b1*(a0*am1-a2-am1*a1+a1^2+am2-a1*a0)/(b1^2-a1*b0*b1+a2*b0^2)+...
    b0*(am1*a2-a1*a2-a0*am2+a0*a2)/(b1^2-a1*b0*b1+a2*b0^2);

est_s1_2 = b1*(a1*a2-am1*a2+a0*am2-a0*a2)/(b1^2-a1*b0*b1+a2*b0^2)+...
     b0*(a2*am2-a2^2-a0*am2*a1+a0*a2*am1)/(b1^2-a1*b0*b1+a2*b0^2);

est_R2 = tf([1 est_r1], 1, h,'variable', 'q');
est_S2 = tf([est_s0_2 est_s1_2], 1, h,'variable', 'q');
est_T2 = tf([bm0/b0], 1, h,'variable', 'q');

