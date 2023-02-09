%% generating parameters
clear all
close all

for i = 1:450
   phi(i)=0.25*rand(1,1)+0.1;
   M(i)=0.30*rand(1,1)+0.05; 
   T(i) = 1150+rand(1,1)*80;
   
   A = rand(1,1)*0.25;
   B = rand(1,1)*0.70;
   D = rand(1,1)*0.05;
   
   O(i) = A/(A+B+D);
   P(i) = B/(A+B+D);
   C(i) = D/(A+B+D);
end

Results(:,1) = T;
Results(:,2) = phi;
Results(:,3) = M;
Results(:,4) = O;
Results(:,5) = P;
Results(:,6) = C;

save('input_scenario3.mat', 'Results')