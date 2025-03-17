function dydt=NPFEM_H(t,y,X0,Dnum_pnode,Dnum_Anode,num_elem,Nnum_pelem,Elem,ke_all,me_all,Ge_all,fk_all,alpha_m,beta_k)
% Hamiltonian function of NPFEM
% by Fuzhen Yao at 2024/04/17
% t1=tic;

X=y(1:Dnum_Anode); 
P=y(Dnum_Anode+1:end); 
% constrains
dw=10000; % acceleration and time step
if t<0.05
theta=dw*t^2/2;
dtheta=dw*t;
else 
dtheta=dw*0.05;
theta=dw*0.05^2/2+dtheta*(t-0.05);
end
% dtheta=100;
% theta=dtheta*t;
A_theta=[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];
Te=A_theta';

% assemble
[M,C,K,FK]=MCKF(X,Te,Dnum_pnode,Dnum_Anode,num_elem,Nnum_pelem,Elem,ke_all,me_all,Ge_all,fk_all,alpha_m,beta_k);

% 初始化集中质量矩阵
M_lumped = zeros(size(M));
% 计算每个节点的集中质量
for i = 1:length(M)
    M_lumped(i,i) = sum(M(i,:));
end
M= M_lumped;

dX=M\P;%zeros(Dnum_Anode,1);%
for i=1:6
    X(i*3-2:i*3)=A_theta*X0(i*3-2:i*3);
    dX(i*3-2:i*3)=cross([0;0;dtheta],X(i*3-2:i*3));
end

dydt=[dX;
      (-K*X+FK)];%+C*dX+9.8*repmat(v_g',Dnum_Anode/3,1)

% tCount1=toc(t1);
% fprintf(['use ',num2str(tCount1),' seconds.\n'])
end