clc
% global E_case nu_case rho_case g_case % 弹性模量、泊松比、密度
% global num_node Dnum_pnode Dnum_Anode num_elem Nnum_pelem Dnum_pelem  %单元，节点信息 Elem为单元
% global ke_all me_all Ge_all fk_all % Arr_amp Vec_fac_F% 局部坐标系下
% global alpha_m beta_k %材料阻尼比例系数
% 材料及环境参数：
E_case=113e9; rho_case=4430; nu_case=0.35; %g_case=9.8*1; v_g=[0 -1 0]*0;% 重力向量
%
omega_1=50; omega_2=1500; epsilon_1=0.05; epsilon_2=0.08;  % for the alpha_m & beta_k
[alpha_m,beta_k]=coe_MAT(omega_1,omega_2,epsilon_1,epsilon_2);  beta_k=1*beta_k; alpha_m=1*alpha_m;
%% ready for 3D element
I_forQ=eye(3); A=zeros(8); A(:,1)=1; Q=kron(A,I_forQ);
Node=importdata('NODE.txt');  Elem=importdata('ELEM.txt'); NODE_side=importdata('NODE_side.txt'); Gauss_points=importdata('Gauss points.txt');
[num_node,Dnum_pnode]=size(Node); Dnum_pnode=Dnum_pnode-1; Dnum_Anode=num_node*Dnum_pnode; %节点数；每节点自由度数; 总自由度数
[num_elem,Nnum_pelem]=size(Elem); Nnum_pelem=Nnum_pelem-1; Dnum_pelem=Nnum_pelem*Dnum_pnode; %单元数；每单元节点数；每单元自由度数；

[ke_all,me_all,fk_all,Ge_all]=Nummat_3D(Elem,Node,Q,Gauss_points,E_case,nu_case,rho_case,Dnum_pnode,Dnum_Anode,num_elem,num_node,Nnum_pelem,Dnum_pelem);

%% Hamiltonian solution
% initial parameter
% Rotating installation
alpha=0*pi/180; % 有转角时记得开始实时转换矩阵运算 in MCKF
A_alpha=[1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
X0=zeros(size(Node)-[0,1]);
for i=1:size(Node,1)
    X0(i,:)=(A_alpha*Node(i,2:4)')';
end
X0=X0.'; X0=X0(:);
Xd0=X0*0; P0=Xd0;
y0=[X0;P0];
tspan=[0 0.1];
t0=tic;
% Linear
[t,y] = ode45(@(t,y) NPFEM_H(t,y,X0,Dnum_pnode,Dnum_Anode,num_elem,Nnum_pelem,Elem,ke_all,me_all,Ge_all,fk_all,alpha_m,beta_k), tspan, y0);
% Nonlinear
% [t,y] = ode45(@(t,y) NPFEM_NL(t,y,X0,Elem,Node,Q,Gauss_points,E_case,nu_case,rho_case,num_node,Dnum_pnode,Dnum_Anode,Nnum_pelem,alpha_m,beta_k), tspan, y0);
tCount0=toc(t0);
fprintf(['use ',num2str(tCount0),' seconds.\n'])

% save all points
step=fix(size(t,1)/10000); %fix(size(t,1)/10);
output=zeros(size(y(1:step:end,:))+[0,1]);
output(:,1)=t(1:step:end,1);
output(:,2:end)=y(1:step:end,:);
save('ALL_output.txt','output')

% Get output(5 middle points 3 9 15 21 27)
% step=fix(size(t,1)/1000);
% output=zeros(size(t(1:step:end,1),1),11);
% output(:,1)=t(1:step:end,1);
% for i=1:5
% output(:,i*2:i*2+1)=y(1:step:end,18*(i-1)*(num_elem/8)+7:18*(i-1)*(num_elem/8)+8);
% end
% for i=1:size(output,1)
%     theta=atan2(output(i,3),output(i,2));
%     A=[cos(theta) sin(theta);-sin(theta) cos(theta)]; %from global to local
%     for j=1:4
%         output(i,j*2+2:j*2+3)=(A*(output(i,j*2+2:j*2+3)-output(i,2:3))')';
%     end
%     output(i,2:3)=[0,0];
% end
% save('output_NPFEM.txt','output')