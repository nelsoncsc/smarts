%% Task 3: Examples 3.{4;5;7;8;9;10} and Questions 3.{5;11;12;15}  
% Nelson Campos
%
%%

%% Example 3.4: Indirect self-tuner with cancellation of process zero
%
% Let the process be the same as in Example 3.1 and assume that the process 
% zero is canceled. The specifications are the same as in Example 3.1, that is, 
% to obtain a closed-loop characteristic polynomial Am. The parameters of the 
% model $$  y(t) + a_1y(t-1) + a_2y(t-2) = b_0u(t-1) + b_1 u(t-2) $$
% which has the same structure as Eq. (3.17), are estimated by using the leastsquares 
% algorithm. Algorithm 3.2 is used for the self-tuning regulator. The
% calculations, which were done in Example 3.1, give the control law
%
% $$ u(t)+r_1u(t-1) = t_0uc(t)-s_0y(t)-s_1y(t-1) $$
%
% The controller parameters were expressed as functions of the model parameters
% and the specifications. Figure 1 shows the process output and the control
% signal in a simulation of the process with the self-tuner when the command
% signal is a square wave. Notice that the estimates converge quickly.
% The system in Example 3.4 behaves quite well, apart from the “ringing”
% control signal. This can be avoided by using a design in which the process zero
% is not canceled. The consequences of this are illustrated in the next example.
%%


close all
clear all

w = warning ('off','all');

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

%% Example 3.5: Indirect self-tuner without cancellation of process zero
% 
% Consider the same process as in Example 3.4, but use a control design in
% which there is no cancellation of the process zero. The parameters are estimated
% in the same way as in Example 3.4, but the control law is now computed as in
% Example 3.2. Figure 2 shows the process output and the control
% signal in a simulation of the process with the indirect self-tuner when the command
% signal is a square wave. Notice that the behavior of the
% process output is quite similar to that in Fig. 1 but that there is no “ringing”
% in the control signal. 
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

figure(3), subplot(2,1,1), plot(t,uc, t, y), xlabel('Time'), ylabel('Amplitude'), title('y(t) vs t and uc(t) vs t')
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

lambda =  1;

for i = 3 : size(t,2)
          
phi(:,i) = [-y(i-1) -y(i-2) u(i-1) u(i-2)]';
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

