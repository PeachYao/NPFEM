function [alpha_m,beta_k]=coe_MAT(omega_1,omega_2,epsilon_1,epsilon_2)
%% 节计算alpha beta 材料阻尼系数
% By 李朝峰 at 2605-750 Bay ST. Tronto at 2017/10/06
%      ***      revised *** at ？？？
alpha_m=2*(epsilon_1*omega_2-epsilon_2*omega_1)*omega_1*omega_2/(omega_2^2-omega_1^2); %材料阻尼 系数alpha_m， beta_k。
beta_k=2*(epsilon_2*omega_2-epsilon_1*omega_1)/(omega_2^2-omega_1^2);  %材料阻尼 系数alpha_m， beta_k。
end