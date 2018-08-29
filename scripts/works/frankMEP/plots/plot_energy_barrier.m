close all;
clear;
clc;
figure(2);
hold on;

for II = 3 : 3
if II == 1 
dirname = '../../../runs/frankMEP/test_db_0.04/';
spec = '*-';
end
if II == 2 
dirname = '../../../runs/frankMEP/test_db_0.05/';
spec = 'o-';
end
if II == 3 
dirname = '../../../runs/frankMEP/test_db_0.06/';
spec = 's-';
end
		
epot2_file = strcat(dirname,'/potential_new.dat');
%epot1_file = strcat(dirname,'/EPOT_1.dat');
EPOT_2=load(epot2_file);
%EPOT_1=load(epot1_file);

d_void_epot = 4.63 ; % 2.64; % 4.63; %9.04101350863513;4.63
tmp = load(strcat(dirname,'/nucleus-', num2str(0),'.dat'));
numvoids_0 = size(tmp,1);  

valid_disl_index_0 = [1:length(EPOT_2)];

EPOT_2_0 = EPOT_2(valid_disl_index_0);

for i = 1 : length(EPOT_2_0)
  id = valid_disl_index_0(i);
  tmp = load(strcat(dirname,'/nucleus-', num2str(id-1),'.dat'));
  numvoids = size(tmp,1);  
  NumVoids_0(i) = numvoids;
  EPOT_2_offset_0(i) = EPOT_2_0(i)-(numvoids)*d_void_epot;
  display(EPOT_2_0(i));
  display(numvoids);
end

%plot(NumVoids_0, EPOT_2_offset_0(1,:)'-EPOT_2_offset_0(1,1),spec, 'LineWidth', 2, 'markers',10);
plot(NumVoids_0, EPOT_2-EPOT_2(1),spec, 'LineWidth', 2, 'markers',10);

end

xlim([0 220])
legend('0.04','0.05','0.06', 'Location','southwest')
xlabel('Nucleus size')
ylabel('Relative potential energy')
set(gca, 'FontSize', 15) 

print('improvedExample','-depsc2','-r300');
