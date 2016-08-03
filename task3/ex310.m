%% Task 3 - Part IV: Example 3.10
% Nelson Campos
%%

%%
% Example 3.10: Load disturbances: Modified estimator and controller
%
% We now show that the difficulties found in Example 3.9 can be avoided by
% modifying the estimator and the controller. We first introduce a controller
% that has integral action by applying the design procedure that we have just
% described. To do this, we consider the same system as in Example 3.5 where
% the controller was defined by
%
% $$ R^0 = q+r_1 $$
%
% $$ S^0 = s_0q+s_1 $$
%
% The closed-loop characteristic polynomial Ac has degree three. To obtain a
% controller with integral action, the order of the closed-loop system is increased
% by introducing an extra closed-loop pole at q = ?x0 = 0.
% It then follows from Eq. (3.41) that
%
% $$ y_0 = -\frac{1+r_1}{b_0+b_1} $$
%
% The estimates are based on the model (3.42) with Ad = q ? 1 to reduce the
% effects of the disturbances. Figure 2 shows a simulation corresponding to
% Fig. 2 with the modified self-tuning regulator. A comparison with the
% results of example 3.9 shows a significant improvement. The load disturbance is reduced quickly.
% Because of the integral action the control will decrease with a magnitude
% corresponding to the load disturbance shortly after t = 40. The
% parameters estimates are shown in Figure 2, which indicates the advantages in using the
% modified estimator. Notice in particular that there is a very small change in
% the estimates when the load disturbance occurs.
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


u = zeros(1, size(t,2));

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

yf = zeros(1,size(t,2)); % The filtered output 

for i=2:size(t,2)
  yf(i) = y(i)-y(i-1);
end

r1 = b1/b0 + (b1^2-am1*b0*b1+am2*b0^2)*(-b1+a0*b0)/(b0*(b1^2-a1*b0*b1+a2*b0^2));
s0_2 = b1*(a0*am1-a2-am1*a1+a1^2+am2-a1*a0)/(b1^2-a1*b0*b1+a2*b0^2)+...
    b0*(am1*a2-a1*a2-a0*am2+a0*a2)/(b1^2-a1*b0*b1+a2*b0^2);

s1_2 = b1*(a1*a2-am1*a2+a0*am2-a0*a2)/(b1^2-a1*b0*b1+a2*b0^2)+...
     b0*(a2*am2-a2^2-a0*am2*a1+a0*a2*am1)/(b1^2-a1*b0*b1+a2*b0^2);

 y0 = -(1+r1)/(b0+b1);

% R2 = [1 r1];
% S2 = [s0_2 s1_2];
R2 = [1 -(y0*b1+1) y0*b1];
S2 = [(s0_2-y0) (s1_2-y0*a1) (-y0*a2)];

T2 = [bm0/b0];

uf = zeros(1,size(t,2)); % The filtered input 
for i=3:size(t,2)
  u(i) = -R2(2)*u(i-1)-R2(3)*u(i-2)+T2(1)*uc(i)-S2(1)*y_noise(i)-S2(2)*y_noise(i-1)-S2(3)*y_noise(i-2);
  uf(i) = u(i)-u(i-1);
end

figure(1), subplot(2,1,1), plot(t,uc, t, y_noise), xlabel('Time'), ylabel('Amplitude'), title('y_noise(t) vs t and uc(t) vs t')
           subplot(2,1,2), plot(t,uf),, xlabel('Time'), ylabel('Amplitude'), title('u(t) vs t without zero cancellation')

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
          
phi(:,i) = [-yf(i-1) -yf(i-2) uf(i-1) uf(i-2)]';
K = P * phi(:,i)* inv(lambda + phi(:,i)'* P * phi(:,i));
P = ((I - K * phi(:,i)')* P)/lambda;
est_y = phi(:,i)'*theta(:,i-1);
error(i) = yf(i)-est_y;

theta(:,i) = theta(:,i-1) + K * error(i);

end

figure(2), subplot(2,1,1), plot(t,a1,t,theta(1,:),t,a2,t,theta(2,:)), xlabel('Time'), ylabel('Amplitude'), title('Estimated a1 and a2 without zero cancellation'), legend('a1', 'est a1', 'a2', 'est a2')
           subplot(2,1,2), plot(t,b0,t,theta(3,:),t,b1,t,theta(4,:)), xlabel('Time'), ylabel('Amplitude'), title('Estimated b0 and b1 without zero cancellation'), legend('b0', 'est b0', 'b1', 'est b1')
