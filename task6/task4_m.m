% Nelson Campos

clear all
close all

s = tf('s');

H = tf(1, [100 20 1], 'InputDelay', 1);

Kc = 15;
Ti = 6.75;
Td = 1.69;

PID = [Kc Ti Td];
Gc = PID(1)*(1+1/(PID(2)*s)+PID(3)*s);

disp('Features = [PO; OR; damping; T; riseTime]')
F = feedback(Gc*H,1);
[y,t] = step(F);

Features = calculateFeatures(F,y,t);

% From the Ziegler-Nichols Tuning Rule
Ku = 25; 
Tu = 13.5;

Kpmin = 0.55*Ku;
Kpmax = 0.65*Ku;
Kdmin = 0.070*Ku*Tu;
Kdmax = 0.075*Ku*Tu;
Kimin = 1.15*Ku/Tu;
Kimax = 1.30*Ku/Tu;

Kc_clp = load('Kc.txt');
Ti_clp = load('Ti.txt');
Td_clp = load('Td.txt');

Kc2 = Kc_clp(length(Kc_clp));
Ti2 = Ti_clp(length(Ti_clp));
Td2 = Td_clp(length(Td_clp));

Kp_new = Kpmin + (Kpmax-Kpmin)*Kc2;
Ti_new = Kimin + (Kimax-Kimin)*Ti2;
Td_new = Kdmin + (Kdmax-Kdmin)*Td2;

PID_new = [Kp_new Kp_new/Ti_new Td_new/Kp_new];
Gc_new = PID_new(1)*(1+1/(PID_new(2)*s)+PID_new(3)*s);
F_new = feedback(Gc_new*H, 1);

figure(1), step(F, F_new), legend('Original Response', 'Clips Response')
