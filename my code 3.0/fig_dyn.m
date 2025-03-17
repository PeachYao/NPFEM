%%% dynamics vs rotation beam
clc,clear
Node=importdata('NODE.txt'); X0=Node(:,2:end); X0=X0.'; X0=X0(:);
Elem=importdata('ELEM.txt'); [num_elem,Nnum_pelem]=size(Elem);  %单元数；每单元节点数；每单元自由度数；
P_npfem=importdata('ALL_output.txt');
alpha=0*pi/180;
num_node=(size(P_npfem,2)-1)/6;
P_NPFEM=P_npfem(1:end,1:num_node*3+1);

A_alpha=[1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];

dw=10000;
% dtheta=100;
for i=1:size(P_NPFEM,1)
    if P_NPFEM(i,1)<=0.05
        theta=dw*P_NPFEM(i,1)^2/2;
    else
        theta=dw*0.05^2/2+dw*0.05*(P_NPFEM(i,1)-0.05);
    end
%     theta=dtheta*P_NPFEM(i,1);
    A_theta=[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];
    for j=1:num_node
        P_NPFEM(i,j*3-1:j*3+1)=(A_alpha'*A_theta'*(P_NPFEM(i,j*3-1:j*3+1)-(A_theta*A_alpha*[0.5;0;-0.05])')')';
    end
end

figure(1)
for i=1:4
plot(P_NPFEM(:,1),P_NPFEM(:,i*(2+1)*6*(4/4)+3*(2+1)),'r-'); % i*(2+1)*6*(ny/4)+3*(2/2*2+1)
hold on
end

%% angular velocity
X1=0:0.001:0.05;
dw=1000;
Y1=dw.*X1;
X2=0.05:0.001:0.1;
Y2=dw*0.05*ones(size(X2));


%% frequncy
% theory
load MK.mat
M=M_40_C3D8I(6*3+1:end,6*3+1:end);
K=K_40_C3D8I(6*3+1:end,6*3+1:end);

% 计算特征值和特征向量
[V, D] = eig(M\K);

% 计算振动频率
[omega_n,w_order] = sort(sqrt(diag(D))); 
omega_n(1:6)./2./pi

%%
% numerical
NN=2;
for j=1:NN
N=90; % need change with number of ele
dw=1000;
t=P_npfem(:,1);
y=P_NPFEM(:,4*(2+1)*6*(4/4)+3*(2+1)+j-2);  % i*(2+1)*6*(ny/4)+3*(2/2*2+1)

% TIME domain
figure(2)
subplot(NN,1,j)
yyaxis left
plot(t,y); hold on
% axis([0,0.1,-8e-4,8e-4]);
xlabel('t (s)');
ylabel('y_b (m)');

yyaxis right
plot(X1,Y1,'r-.'); hold on
plot(X2,Y2,'r-.');
axis([0,0.1,0,100]);
ylabel('\Omega (rad/s)');

% Frequency domain
Fs = 1e6;           % 采样频率
T = 1/Fs;          % 采样周期
t_query=0.06:T:0.099;     % 时间向量
L = length(t_query)-1;   % 信号长度
y_query = interp1(t, y, t_query, 'spline');

figure(3)
plot(t,P_npfem(:,4*(2+1)*6*(4/4)+3*(2+1)+j-2)); hold on
xlabel('t (s)');
ylabel('y (m)');

% 计算信号的FFT
Y = fft(y_query);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% 计算频率向量
f = Fs*(0:(L/2))/L;

% 绘制FFT的复数模
figure(4)
subplot(NN,1,j)
plot(f(1:300), P1(1:300), 'LineWidth', 1); hold on
plot(f(9), P1(9),'bo'); hold on
plot(f(51), P1(51),'bo'); hold on
title('Complex Magnitude of FFT Spectrum');
xlabel('f (Hz)');
ylabel('|FFT(X)|');
grid on
% plot(ans,zeros(6,1),'ro')
end


%% ABAQUS
% % numerical
% P_ABAQUS=importdata('ABAQUS_N_8.txt');
% t=P_ABAQUS(:,1);
% X=zeros(size(t,1),2);
% dw=1000;
% for i=1:size(t,1)
%     if t(i)<=0.05
%         theta=dw*t(i)^2/2;
%     else
%         theta=dw*0.05^2/2+dw*0.05*(t(i)-0.05);
%     end
%     A_theta=[cos(theta) -sin(theta);sin(theta) cos(theta);];
%     X(i,:)=(A_theta'*(P_ABAQUS(i,2:3)+[0.5,0]-(A_theta*[0.5;0])')')';
% end
% y=X(:,2);  
% 
% % TIME domain
% figure(3)
% yyaxis left
% plot(t,y,'m--'); hold on
% axis([0,0.1,-8e-4,8e-4]);
% xlabel('t (s)');
% ylabel('y_b (m)');
% legend('NPFEM','ABAQUS','\Omega')
% 
% % Frequency domain
% Fs = 1e6;           % 采样频率
% T = 1/Fs;          % 采样周期
% t_query=0.06:T:0.099;     % 时间向量
% L = length(t_query);   % 信号长度
% y_query = interp1(t, y, t_query, 'linear');
% 
% % 计算信号的FFT
% Y = fft(y_query);
% 
% % 计算频率向量
% f = Fs/L * (0:L-1);
% 
% % 绘制FFT的复数模
% figure(4)
% plot(f(1:100), abs(Y(1:100)),'m--', 'LineWidth', 1);
% title('Complex Magnitude of FFT Spectrum');
% xlabel('f (Hz)');
% ylabel('|FFT(X)|');
% legend('NPFEM','ABAQUS')

%% rotation period
% nn=0.02;
% figure(2)
% % 定义圆柱体的参数
% radius = 0.3; % 圆柱体的半径
% height = 0.1; % 圆柱体的高度
% n = 100; % 圆柱体周围的分段数
% k = 1;
% % 使用 cylinder 函数创建圆柱体
% [X, Y, Z] = cylinder(radius, n); hold on
% % 将圆柱体高度拉伸为指定高度
% Z = Z * height/2;
% Z(1,:)=-Z(2,:);
% % 绘制圆柱体
% surf(X, Y, Z);
% 
% for i=1:size(P_NPFEM,1)
% if P_npfem(i,1)>nn
% nn=nn+0.02; k=k+1;
% for j=1:size(Elem,1)
% plot3(P_npfem(i,(Elem(j,[2,4])-1)*3+2),P_npfem(i,(Elem(j,[2,4])-1)*3+3),P_npfem(i,(Elem(j,[2,4])-1)*3+4),'r');hold on
% plot3(P_npfem(i,(Elem(j,[3,5])-1)*3+2),P_npfem(i,(Elem(j,[3,5])-1)*3+3),P_npfem(i,(Elem(j,[3,5])-1)*3+4),'r');hold on
% plot3(P_npfem(i,(Elem(j,[6,8])-1)*3+2),P_npfem(i,(Elem(j,[6,8])-1)*3+3),P_npfem(i,(Elem(j,[6,8])-1)*3+4),'r');hold on
% plot3(P_npfem(i,(Elem(j,[7,9])-1)*3+2),P_npfem(i,(Elem(j,[7,9])-1)*3+3),P_npfem(i,(Elem(j,[7,9])-1)*3+4),'r');hold on
% 
% plot3(P_npfem(i,(Elem(j,[2,3])-1)*3+2),P_npfem(i,(Elem(j,[2,3])-1)*3+3),P_npfem(i,(Elem(j,[2,3])-1)*3+4),'r');hold on
% plot3(P_npfem(i,(Elem(j,[4,5])-1)*3+2),P_npfem(i,(Elem(j,[4,5])-1)*3+3),P_npfem(i,(Elem(j,[4,5])-1)*3+4),'r');hold on
% plot3(P_npfem(i,(Elem(j,[6,7])-1)*3+2),P_npfem(i,(Elem(j,[6,7])-1)*3+3),P_npfem(i,(Elem(j,[6,7])-1)*3+4),'r');hold on
% plot3(P_npfem(i,(Elem(j,[8,9])-1)*3+2),P_npfem(i,(Elem(j,[8,9])-1)*3+3),P_npfem(i,(Elem(j,[8,9])-1)*3+4),'r');hold on
% 
% plot3(P_npfem(i,(Elem(j,[2,6])-1)*3+2),P_npfem(i,(Elem(j,[2,6])-1)*3+3),P_npfem(i,(Elem(j,[2,6])-1)*3+4),'r');hold on
% plot3(P_npfem(i,(Elem(j,[3,7])-1)*3+2),P_npfem(i,(Elem(j,[3,7])-1)*3+3),P_npfem(i,(Elem(j,[3,7])-1)*3+4),'r');hold on
% plot3(P_npfem(i,(Elem(j,[4,8])-1)*3+2),P_npfem(i,(Elem(j,[4,8])-1)*3+3),P_npfem(i,(Elem(j,[4,8])-1)*3+4),'r');hold on
% plot3(P_npfem(i,(Elem(j,[5,9])-1)*3+2),P_npfem(i,(Elem(j,[5,9])-1)*3+3),P_npfem(i,(Elem(j,[5,9])-1)*3+4),'r');hold on
% end
% axis equal
% end
% if P_npfem(i,1)>0.1; break; end
% end

%% Von Mises stress
% nn=0.02; k=1; 
% n_x=4; n_y=2; n_z=1; % number of elements on x y z axis
% l_x=0.2/n_x; l_y=0.01/n_y; l_z=0.1/n_z; % length of the element on x y z axis
% I_forQ=eye(3); A=zeros(8); A(:,1)=1; Q=kron(A,I_forQ); E=113e9;
% Gauss_points=importdata('Gauss points.txt');
% x0=Node(Elem(1,2:9),2:4)-Node(repmat(Elem(1,2),1,8),2:4);  %获得每个单元的x0
% figure(2)
% for i=1:size(P_NPFEM,1)
% if P_NPFEM(i,1)>nn
% nn=nn+0.02; k=k+1;
% node_VonMisesStress=zeros(num_node,1);
% for j=1:num_elem
% %求转换矩阵：
% nodes=Elem(j,2:end);
% nodeA_posi=(nodes(1)-1)*3+(1:3); nodeB_posi=(nodes(5)-1)*3+(1:3); nodeC_posi=(nodes(3)-1)*3+(1:3); %用1 3 5点计算转换矩阵
% pA=P_npfem( i,nodeA_posi+1 ) ;  pB=P_npfem( i,nodeB_posi+1 ) ;  pC=P_npfem( i,nodeC_posi+1 ) ;
% Te=transMatr(pA,pB,pC);    Te=reshape(Te,3,3)'; % 生成整体转换矩阵
% pe=zeros(24,1);
% for jj=1:8; pe(jj*3-2:jj*3)=Te*(P_npfem(i,(nodes(jj)-1)*3+(1:3)+1)-P_npfem(i,(nodes(1)-1)*3+(1:3)+1))'; end
% u=pe-[0;0;0;0;0;l_z;0;l_y;0;0;l_y;l_z;l_x;0;0;l_x;0;l_z;l_x;l_y;0;l_x;l_y;l_z];
% VonMisesStress=zeros(8,1);
% for ii=1:8
% % Linear B
% for gi=1:length(Gauss_points)
% x=Gauss_points(gi,2);
% y=Gauss_points(gi,3);
% z=Gauss_points(gi,4);
% % jacob
% J=[(x0(2,1)*(y-1)*(z+1))-(x0(1,1)*(y-1)*(z-1))+(x0(3,1)*(y+1)*(z-1))-(x0(4,1)*(y+1)*(z+1))+(x0(5,1)*(y-1)*(z-1))-(x0(6,1)*(y-1)*(z+1))-(x0(7,1)*(y+1)*(z-1))+(x0(8,1)*(y+1)*(z+1)),(x0(2,2)*(y-1)*(z+1))-(x0(1,2)*(y-1)*(z-1))+(x0(3,2)*(y+1)*(z-1))-(x0(4,2)*(y+1)*(z+1))+(x0(5,2)*(y-1)*(z-1))-(x0(6,2)*(y-1)*(z+1))-(x0(7,2)*(y+1)*(z-1))+(x0(8,2)*(y+1)*(z+1)),(x0(2,3)*(y-1)*(z+1))-(x0(1,3)*(y-1)*(z-1))+(x0(3,3)*(y+1)*(z-1))-(x0(4,3)*(y+1)*(z+1))+(x0(5,3)*(y-1)*(z-1))-(x0(6,3)*(y-1)*(z+1))-(x0(7,3)*(y+1)*(z-1))+(x0(8,3)*(y+1)*(z+1));
%    (x0(2,1)*(x-1)*(z+1))-(x0(1,1)*(x-1)*(z-1))+(x0(3,1)*(x-1)*(z-1))-(x0(4,1)*(x-1)*(z+1))+(x0(5,1)*(x+1)*(z-1))-(x0(6,1)*(x+1)*(z+1))-(x0(7,1)*(x+1)*(z-1))+(x0(8,1)*(x+1)*(z+1)),(x0(2,2)*(x-1)*(z+1))-(x0(1,2)*(x-1)*(z-1))+(x0(3,2)*(x-1)*(z-1))-(x0(4,2)*(x-1)*(z+1))+(x0(5,2)*(x+1)*(z-1))-(x0(6,2)*(x+1)*(z+1))-(x0(7,2)*(x+1)*(z-1))+(x0(8,2)*(x+1)*(z+1)),(x0(2,3)*(x-1)*(z+1))-(x0(1,3)*(x-1)*(z-1))+(x0(3,3)*(x-1)*(z-1))-(x0(4,3)*(x-1)*(z+1))+(x0(5,3)*(x+1)*(z-1))-(x0(6,3)*(x+1)*(z+1))-(x0(7,3)*(x+1)*(z-1))+(x0(8,3)*(x+1)*(z+1));
%    (x0(2,1)*(x-1)*(y-1))-(x0(1,1)*(x-1)*(y-1))+(x0(3,1)*(x-1)*(y+1))-(x0(4,1)*(x-1)*(y+1))+(x0(5,1)*(x+1)*(y-1))-(x0(6,1)*(x+1)*(y-1))-(x0(7,1)*(x+1)*(y+1))+(x0(8,1)*(x+1)*(y+1)),(x0(2,2)*(x-1)*(y-1))-(x0(1,2)*(x-1)*(y-1))+(x0(3,2)*(x-1)*(y+1))-(x0(4,2)*(x-1)*(y+1))+(x0(5,2)*(x+1)*(y-1))-(x0(6,2)*(x+1)*(y-1))-(x0(7,2)*(x+1)*(y+1))+(x0(8,2)*(x+1)*(y+1)),(x0(2,3)*(x-1)*(y-1))-(x0(1,3)*(x-1)*(y-1))+(x0(3,3)*(x-1)*(y+1))-(x0(4,3)*(x-1)*(y+1))+(x0(5,3)*(x+1)*(y-1))-(x0(6,3)*(x+1)*(y-1))-(x0(7,3)*(x+1)*(y+1))+(x0(8,3)*(x+1)*(y+1));]/8;
% 
% DN_XI=[-((y - 1)*(z - 1))/8, ((y - 1)*(z + 1))/8, ((y + 1)*(z - 1))/8, -((y + 1)*(z + 1))/8, ((y - 1)*(z - 1))/8, -((y - 1)*(z + 1))/8, -((y + 1)*(z - 1))/8, ((y + 1)*(z + 1))/8;
%        -((x - 1)*(z - 1))/8, ((x - 1)*(z + 1))/8, ((x - 1)*(z - 1))/8, -((x - 1)*(z + 1))/8, ((x + 1)*(z - 1))/8, -((x + 1)*(z + 1))/8, -((x + 1)*(z - 1))/8, ((x + 1)*(z + 1))/8;
%        -((x - 1)*(y - 1))/8, ((x - 1)*(y - 1))/8, ((x - 1)*(y + 1))/8, -((x - 1)*(y + 1))/8, ((x + 1)*(y - 1))/8, -((x + 1)*(y - 1))/8, -((x + 1)*(y + 1))/8, ((x + 1)*(y + 1))/8];
% 
% DN_X = J\ DN_XI;
% 
% B=[DN_X(1,1),        0,        0,DN_X(1,2),        0,        0,DN_X(1,3),        0,        0,DN_X(1,4),        0,        0,DN_X(1,5),        0,        0,DN_X(1,6),        0,        0,DN_X(1,7),        0,        0,DN_X(1,8),        0,        0;
%             0,DN_X(2,1),        0,        0,DN_X(2,2),        0,        0,DN_X(2,3),        0,        0,DN_X(2,4),        0,        0,DN_X(2,5),        0,        0,DN_X(2,6),        0,        0,DN_X(2,7),        0,        0,DN_X(2,8),        0;
%             0,        0,DN_X(3,1),        0,        0,DN_X(3,2),        0,        0,DN_X(3,3),        0,        0,DN_X(3,4),        0,        0,DN_X(3,5),        0,        0,DN_X(3,6),        0,        0,DN_X(3,7),        0,        0,DN_X(3,8);
%     DN_X(2,1),DN_X(1,1),        0,DN_X(2,2),DN_X(1,2),        0,DN_X(2,3),DN_X(1,3),        0,DN_X(2,4),DN_X(1,4),        0,DN_X(2,5),DN_X(1,5),        0,DN_X(2,6),DN_X(1,6),        0,DN_X(2,7),DN_X(1,7),        0,DN_X(2,8),DN_X(1,8),        0;
%             0,DN_X(3,1),DN_X(2,1),        0,DN_X(3,2),DN_X(2,2),        0,DN_X(3,3),DN_X(2,3),        0,DN_X(3,4),DN_X(2,4),        0,DN_X(3,5),DN_X(2,5),        0,DN_X(3,6),DN_X(2,6),        0,DN_X(3,7),DN_X(2,7),        0,DN_X(3,8),DN_X(2,8);
%     DN_X(3,1),        0,DN_X(1,1),DN_X(3,2),        0,DN_X(1,2),DN_X(3,3),        0,DN_X(1,3),DN_X(3,4),        0,DN_X(1,4),DN_X(3,5),        0,DN_X(1,5),DN_X(3,6),        0,DN_X(1,6),DN_X(3,7),        0,DN_X(1,7),DN_X(3,8),        0,DN_X(1,8)];
% end
% % Nonlinear B
% % B=B+BN;
% 
% Stress=E*B*u;
% VonMisesStress(ii)=sqrt(Stress(1)^2-Stress(1)*Stress(2)+Stress(2)^2+3/4*Stress(4)^2);
% end
% node_VonMisesStress(nodes)=node_VonMisesStress(nodes)+VonMisesStress; % get stress
% end
% result = tabulate(reshape(Elem(:,2:end),[],1));
% node_VonMisesStress=node_VonMisesStress./(result(:,2));
% subplot(2,2,k-1)
% v=reshape(P_NPFEM(i,2:end),3,num_node)'; v=v(:,1:2);
% f=Elem(:,[2 4 8 6]);
% col=node_VonMisesStress;
% patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','interp');
% colorbar
% % axis([0,0.225,-0.008,0.002])
% title(['Von Mises Stress at ', num2str(nn-0.02), ' s']);
% xlabel('local u-axis');
% ylabel('local v-axis');
% % axis equal
% end
% if P_NPFEM(i,1)>0.1; break; end
% end