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
figure(1), step(F)


%sim('task4.mdl')