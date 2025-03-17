function [ke_all,me_all,fk_all,Ge_all]=Nummat_3D(Elem,Node,Q,Gauss_points,E_case,nu_case,rho_case,Dnum_pnode,Dnum_Anode,num_elem,num_node,Nnum_pelem,Dnum_pelem)
% Substitute parameters into symbolic matrix for every element
ke_all=zeros(Dnum_pelem,Dnum_pelem,num_elem); me_all=zeros(Dnum_pelem,Dnum_pelem,num_elem);
Ge_all=zeros(Dnum_pelem,Dnum_Anode,num_elem);fk_all=zeros(Dnum_pelem,num_elem);   
for m=1:num_elem
    [ke_all_m,me_all_m,Ge_all_m,fk_all_m]=Fun_me_ke_fk(Elem,Node,m,Gauss_points,E_case,nu_case,rho_case,num_node,Dnum_pnode,Nnum_pelem);
    ke_all(:,:,m)=(eye(24)-Q)'*ke_all_m*(eye(24)-Q); 
    me_all(:,:,m)=me_all_m; 
    Ge_all(:,:,m)=Ge_all_m; 
    fk_all(:,m)=(eye(24)-Q)'*fk_all_m;      % equivalent nodal force
    % lump ele
%     me_all(:,:,m)=diag(sum(me_all_m, 2));
end

end

