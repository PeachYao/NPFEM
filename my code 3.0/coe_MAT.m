function [alpha_m,beta_k]=coe_MAT(omega_1,omega_2,epsilon_1,epsilon_2)
%% �ڼ���alpha beta ��������ϵ��
% By ��� at 2605-750 Bay ST. Tronto at 2017/10/06
%      ***      revised *** at ������
alpha_m=2*(epsilon_1*omega_2-epsilon_2*omega_1)*omega_1*omega_2/(omega_2^2-omega_1^2); %�������� ϵ��alpha_m�� beta_k��
beta_k=2*(epsilon_2*omega_2-epsilon_1*omega_1)/(omega_2^2-omega_1^2);  %�������� ϵ��alpha_m�� beta_k��
end