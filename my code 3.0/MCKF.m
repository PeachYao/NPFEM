function [M,C,K,FK]=MCKF(X,Te,Dnum_pnode,Dnum_Anode,num_elem,Nnum_pelem,Elem,ke_all,me_all,Ge_all,fk_all,alpha_m,beta_k)
%用于生成计算过程中的整体质量、刚度、张力扩展向量以及个各单元的局部-整体坐标转换矩阵
% global Dnum_pnode Dnum_Anode num_elem Nnum_pelem Elem
% global ke_all me_all Ge_all fk_all 
% global alpha_m beta_k
M=zeros(Dnum_Anode); K=zeros(Dnum_Anode); FK=zeros(Dnum_Anode,1);   
T= kron(eye(Nnum_pelem), Te); % 生成整体转换矩阵 
for m=1:num_elem   %当前单元编号

    %求转换矩阵：
    % code_nodes=Elem(m,2:end); %取出第m个单元的节点号；获得第m个单元的所有节点的坐标向量
    % nodeA_posi=(code_nodes(1)-1)*Dnum_pnode+(1:Dnum_pnode); nodeB_posi=(code_nodes(5)-1)*Dnum_pnode+(1:Dnum_pnode); nodeC_posi=(code_nodes(3)-1)*Dnum_pnode+(1:Dnum_pnode); %用1 5 3点计算转换矩阵
    % pA=X( nodeA_posi ) ;  pB=X( nodeB_posi ) ;  pC=X( nodeC_posi ) ;
    % Te=transMatr(pA,pB,pC); T= kron(eye(Nnum_pelem), Te); % 生成整体转换矩阵   

       
    %
    K_elm=T'*ke_all(:,:,m)*T ;    
    Fk_elm=T'*fk_all(:,m) ; 
    M_elm=T'*me_all(:,:,m)*T;
    Ge=Ge_all(:,:,m);%获取当前转换矩阵
    M=M+Ge'*M_elm*Ge;
    K=K+Ge'*K_elm*Ge;
    FK=FK+Ge'*Fk_elm;
end
C=alpha_m*M+beta_k*K;

end