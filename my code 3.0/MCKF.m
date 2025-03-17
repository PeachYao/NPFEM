function [M,C,K,FK]=MCKF(X,Te,Dnum_pnode,Dnum_Anode,num_elem,Nnum_pelem,Elem,ke_all,me_all,Ge_all,fk_all,alpha_m,beta_k)
%�������ɼ�������е������������նȡ�������չ�����Լ�������Ԫ�ľֲ�-��������ת������
% global Dnum_pnode Dnum_Anode num_elem Nnum_pelem Elem
% global ke_all me_all Ge_all fk_all 
% global alpha_m beta_k
M=zeros(Dnum_Anode); K=zeros(Dnum_Anode); FK=zeros(Dnum_Anode,1);   
T= kron(eye(Nnum_pelem), Te); % ��������ת������ 
for m=1:num_elem   %��ǰ��Ԫ���

    %��ת������
    % code_nodes=Elem(m,2:end); %ȡ����m����Ԫ�Ľڵ�ţ���õ�m����Ԫ�����нڵ����������
    % nodeA_posi=(code_nodes(1)-1)*Dnum_pnode+(1:Dnum_pnode); nodeB_posi=(code_nodes(5)-1)*Dnum_pnode+(1:Dnum_pnode); nodeC_posi=(code_nodes(3)-1)*Dnum_pnode+(1:Dnum_pnode); %��1 5 3�����ת������
    % pA=X( nodeA_posi ) ;  pB=X( nodeB_posi ) ;  pC=X( nodeC_posi ) ;
    % Te=transMatr(pA,pB,pC); T= kron(eye(Nnum_pelem), Te); % ��������ת������   

       
    %
    K_elm=T'*ke_all(:,:,m)*T ;    
    Fk_elm=T'*fk_all(:,m) ; 
    M_elm=T'*me_all(:,:,m)*T;
    Ge=Ge_all(:,:,m);%��ȡ��ǰת������
    M=M+Ge'*M_elm*Ge;
    K=K+Ge'*K_elm*Ge;
    FK=FK+Ge'*Fk_elm;
end
C=alpha_m*M+beta_k*K;

end