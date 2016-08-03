function Features = calculateFeatures(F, y, t)
N = size(t,1);
dy = zeros(N,1);
info = stepinfo(F);
zero_cross = 0;
t_cross = 0;
index = 0;
 

for i=2:N
  dy(i) = y(i)-y(i-1);
end


for i=1:N
    if abs(y(i)-1) < abs(dy(i))/2
        t_cross(i) = t(i);
        index(i) = i;
        zero_cross = zero_cross+1;
    end
end

%Remove all the elements which values are 0
t_cross = t_cross(t_cross ~=0);
index = index(index ~=0);

P1 = max(y);
PO = 100*(P1-1);
riseTime = info.RiseTime;
  
if size(index,2) == 2 || size(index,2) == 3 
    T = NaN;
    P2 = 0;
    OR = 0;
    damping = 0;
    disp('System without oscillation')
  
else if size(index,2) > 3 % If the system has oscillation
  T = t_cross(3)-t_cross(1);
  t0 = index(1);
  t1 = index(3);
  size(index)
  t2 = index(4);
  
  V1 = min(y(t0:t1));  
  P2 = max(y(t1:t2));
  OR = (P2-1)/(P1-1);
  damping = (P2-V1)/(P1-V1);
  
  if OR < 1e-6
      disp('System Unstable')
  else disp('System with oscillation')
  end

else
    T = NaN;
    PO = NaN;
    OR = NaN;
    damping = NaN;
    riseTime = NaN;
    disp('Unstable System')
end
end
Features = [PO; OR; damping; T; riseTime];
%csvwrite('data.csv', Features);
