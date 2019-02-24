%% GDA geometric dynamics algorithm
%
%  Authors: Mohammad Safeea
%  University of Coimbra, Coimbra, Portugal
%  Ensam, ParisTech, Lille, France
%  2015-11-09
%  ------------
% 
%  This function is used to calculate the joint space Coriolis matrix (B) of
%  the equation ((t= A.q" + B.dq')) of a robotic manipulator with revolute 
%  joints using GDA algorithm.
%  modified Denavit Hartenberg convention was utilized

function [Bt]=GetCoriolisMatrixGDA(T,Pcii,Icii,mcii,dq)
% Bt: is the joint space coriolis matrix
% T: is 4x4xn matrix, representing the homogeneous transformations of the
% link frames
% Thus matrix T(:,:,i) represents 4x4 homogeneous transform of frame i with
% respect to reference frame
% Pcii: is 3xn matrix, each column Pcii(:,i) represent the coordinate
% vector of the center of mass of link i in the local frame of that link.
% Icii: is 3x3xn matrix, thus Icii(:,:,i) matrix represents the
% inertia of link i, around its center of mass represented in frame i
% mcii: is nx1 column vector, while each element mcii(i) represents 
% the mass of link i
% dq: is nx1 column vector, while each element dq(i) represents the angular
% velocity of joint i.

n=max(size(mcii));
%% Initialization 
Bi=zeros(3,n,n);
Di=zeros(3,n,n);
L=zeros(3,3);
Lj_1=zeros(3,1);
wj=zeros(3,n);
half_wj=zeros(3,n);
%% Calculating some auxuliary variables
Pcii_A=zeros(3,n);
w=zeros(3,n);
vci=zeros(3,n);
vi=zeros(3,n);
mcii_Pcii_A=zeros(3,n);
Pcii_A(:,1)=T(1:3,1:3,1)*Pcii(:,1);
w(:,1)=T(1:3,3,1)*dq(1);
wj(:,1)=w(:,1);
half_wj(:,1)=0.5*wj(:,1);
mcii_Pcii_A(:,1)=mcii(1)*Pcii_A(:,1);
double_kj=zeros(3,n);
double_kj(:,1)=2*T(1:3,3,1);
for i=2:n
        Pcii_A(:,i)=T(1:3,1:3,i)*Pcii(:,i);
        wj(:,i)=T(1:3,3,i)*dq(i);
        half_wj(:,i)=0.5*wj(:,i);
        w(:,i)=w(:,i-1)+wj(:,i);
        double_kj(:,i)=2*T(1:3,3,i);
        mcii_Pcii_A(:,i)=mcii(i)*Pcii_A(:,i);
end
%% Run the algorithm
for i=1:n
    Pci=Pcii_A(:,i)+T(1:3,4,i);
    L=T(1:3,1:3,i)*(trace(Icii(:,:,i))*eye(3)-2*Icii(:,:,i))*T(1:3,1:3,i)';
    Lj_1=L*T(1:3,3,i);
    for j=i:-1:2
    Bi(:,j,i)=Bi(:,j,i)+cross1(Lj_1,half_wj(:,j));
    Lj_1=L*T(1:3,3,j-1);
    Bi(:,j-1,i)=Bi(:,j-1,i)+cross1(Lj_1,w(:,i)-w(:,j-1));   
    end
    j=1;    
    Bi(:,j,i)=Bi(:,j,i)+cross1(Lj_1,half_wj(:,j)); 
    vr=zeros(3,1);
    for j=i:-1:2
        Pcij=Pci-T(1:3,4,j);
        Di(:,j,i)=Di(:,j,i)+(dq(j))*(T(1:3,3,j)*(T(1:3,3,j)'*Pcij)-Pcij);
        vr=vr+cross1(wj(:,j),Pcij);
        Di(:,j-1,i)=Di(:,j-1,i)+cross1(double_kj(:,j-1),vr);
    end
    j=1;
    Pcij=Pci-T(1:3,4,j);
    Di(:,j,i)=Di(:,j,i)+(dq(j))*(T(1:3,3,j)*(T(1:3,3,j)'*Pcij)-Pcij);
end
Bt=zeros(n,n);
Fac_D=zeros(3,n);
Mac_B=zeros(3,n);
Pjp1_j=zeros(3,1);

start=n-1;
j=n;
%% recursive rprocedure on moments and forces
% of the last link
        for k=1:n % on moments
            Mac_B(:,k)=Bi(:,k,j)+cross1(mcii_Pcii_A(:,j),Di(:,k,j));            
        end
        Fac_D=mcii(j)*Di(:,:,j); % on forces     
        Bt(j,:)=T(1:3,3,j)'*Mac_B;
        
% of the remaining links
for j=start:-1:1 % iterate through the links
    Pjp1_j=T(1:3,4,j+1)-T(1:3,4,j);
        for k=1:j % on moments  
            Mac_B(:,k)=Mac_B(:,k)+Bi(:,k,j)+cross1(Pjp1_j,Fac_D(:,k))+cross1(mcii_Pcii_A(:,j),Di(:,k,j));            
        end
        for k=j+1:n % on moments
            Mac_B(:,k)=Mac_B(:,k)+cross1(Pjp1_j,Fac_D(:,k));
        end
    Fac_D(:,1:j)=Fac_D(:,1:j)+mcii(j)*Di(:,1:j,j); % on forces          
    Bt(j,:)=T(1:3,3,j)'*Mac_B;
    
end
end

%% Cross product calculation
function c=cross1(a,b)
c = [a(2,:).*b(3,:)-a(3,:).*b(2,:);
     a(3,:).*b(1,:)-a(1,:).*b(3,:);
     a(1,:).*b(2,:)-a(2,:).*b(1,:)];
end
