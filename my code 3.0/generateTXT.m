% Node=fopen('NODE.txt','wt'); ELEM=fopen('ELEM.txt','wt');
% set
n_x=20; n_y=2; n_z=1; % number of elements on x y z axis
l_x=0.2/n_x; l_y=0.01/n_y; l_z=0.1/n_z; % length of the element on x y z axis
x0=0.3; y0=-0.005; z0=-0.05; % initial position
% calculate
N_node=(n_x+1)*(n_y+1)*(n_z+1); % total node
N_elem=(n_x)*(n_y)*(n_z); % total element
Node=zeros(N_node,4); Elem=zeros(N_elem,9); NODE_side=zeros(N_node,4);
% NODE (x-y-z)
account=1;
for i=1:n_x+1
    for j=1:n_y+1
        for k=1:n_z+1
            n_node=(n_y+1)*(n_z+1)*(i-1)+(n_z+1)*(j-1)+k;
            Node(n_node,:)=[n_node,x0+l_x*(i-1),y0+l_y*(j-1),z0+l_z*(k-1)]; 
            if i==1; NODE_side=n_node; end
        end
    end
end
% NODE_side=[1:NODE_side];
save('NODE.txt','Node')
save('NODE_side.txt','NODE_side')
% ELEM (x-y-z)
for i=1:n_x
    for j=1:n_y
        for k=1:n_z
            n_elem=(n_y)*(n_z)*(i-1)+(n_z)*(j-1)+k;
            Elem(n_elem,:)=[n_elem;
                            (n_y+1)*(n_z+1)*(i-1)+(n_z+1)*(j-1)+k;
                            (n_y+1)*(n_z+1)*(i-1)+(n_z+1)*(j-1)+k+1;
                            (n_y+1)*(n_z+1)*(i-1)+(n_z+1)*(j-1)+k+(n_z+1);
                            (n_y+1)*(n_z+1)*(i-1)+(n_z+1)*(j-1)+k+(n_z+1)+1;
                            (n_y+1)*(n_z+1)*(i-1)+(n_z+1)*(j-1)+k+(n_y+1)*(n_z+1);
                            (n_y+1)*(n_z+1)*(i-1)+(n_z+1)*(j-1)+k+(n_y+1)*(n_z+1)+1;
                            (n_y+1)*(n_z+1)*(i-1)+(n_z+1)*(j-1)+k+(n_y+1)*(n_z+1)+(n_z+1);
                            (n_y+1)*(n_z+1)*(i-1)+(n_z+1)*(j-1)+k+(n_y+1)*(n_z+1)+(n_z+1)+1]';   % (x-y-z)
            plot3(Node(Elem(n_elem,2:9),2),Node(Elem(n_elem,2:9),3),Node(Elem(n_elem,2:9),4),'o')
            hold on
        end
    end
end
axis equal
save('ELEM.txt','Elem')