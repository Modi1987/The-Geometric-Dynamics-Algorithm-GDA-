%% GDA geometric dynamics algorithm
%
%  Authors: Mohammad Safeea
%  University of Coimbra, Coimbra, Portugal
%  Ensam, ParisTech, Lille, France
%  2015-11-09
%  ------------
% 
%  This function is used to calculate the joint space Centrifugal matrix    
%  using GDA algorithm.
%  Modified Denavit Hartenberg convention was utilized

function [Cent]=GetCentrifugalMatrixGDA(T,Pcii,Icii,mcii)
% Cent: is nxn joint space Centrifugal matrix, the vector Cent*(dq.^2) is
% the joints torques due to centrifugal effect.
% T: is 4x4xn matrix, representing the homogeneous transformations of the
% link frames
% Thus, matrix T(:,:,i) represents 4x4 homogeneous transform of frame i with
% respect to reference frame
% Pcii: is 3xn matrix, each column Pcii(:,i) represent the coordinate
% vector of the center of mass of link i in the local frame of that link.
% Icii: is 3x3xn matrix, thus each Icii(:,:,i) matrix represents the
% inertia tensor of link i, around its center of mass represented in frame i
% mcii: is nx1 column vector, while each element mcii(i) represents 
% the mass of link i

n=max(size(mcii));
Bi=zeros(3,n,n);
Di=zeros(3,n,n);
Kj=zeros(3,n);
half_Kj=zeros(3,n);
Pcii_A=zeros(3,n);
mcii_Pcii_A=zeros(3,n);
Pcii_A(:,1)=T(1:3,1:3,1)*Pcii(:,1);
mcii_Pcii_A(:,1)=mcii(1)*Pcii_A(:,1);
Kj(:,1)=T(1:3,3,1);
half_Kj(:,1)=0.5*Kj(:,1);
for i=2:n
        Pcii_A(:,i)=T(1:3,1:3,i)*Pcii(:,i);
        mcii_Pcii_A(:,i)=mcii(i)*Pcii_A(:,i);
        Kj(:,i)=T(1:3,3,i);
        half_Kj(:,i)=0.5*Kj(:,i);
end

for i=1:n
    Pci=Pcii_A(:,i)+T(1:3,4,i);
    L=T(1:3,1:3,i)*(trace(Icii(:,:,i))*eye(3)-2*Icii(:,:,i))*T(1:3,1:3,i)';
    for j=1:i
        Bi(:,j,i)=cross1(L*half_Kj(:,j),Kj(:,j));
        Pcij=Pci-T(1:3,4,j);
        Di(:,j,i)=Kj(:,j)*(Kj(:,j)'*Pcij)-Pcij;
    end
end
Cent=zeros(n,n);
Fac_D=zeros(3,n);
Mac_B=zeros(3,n);
Pjp1_j=zeros(3,1);
start=n-1;
j=n;
        for k=1:n 
            Mac_B(:,k)=Bi(:,k,j)+cross1(mcii_Pcii_A(:,j),Di(:,k,j));            
        end
        
        Fac_D=mcii(j)*Di(:,:,j);       
        Cent(j,:)=T(1:3,3,j)'*Mac_B;
    
for j=start:-1:1 
    Pjp1_j=T(1:3,4,j+1)-T(1:3,4,j);
        for k=1:j 
            Mac_B(:,k)=Mac_B(:,k)+Bi(:,k,j)+cross1(Pjp1_j,Fac_D(:,k))+cross1(mcii_Pcii_A(:,j),Di(:,k,j));            
        end
        for k=j+1:n
            Mac_B(:,k)=Mac_B(:,k)+cross1(Pjp1_j,Fac_D(:,k));
        end
    Fac_D(:,1:j)=Fac_D(:,1:j)+mcii(j)*Di(:,1:j,j);        
    Cent(j,:)=T(1:3,3,j)'*Mac_B;
    
end
end
%% Cross product calculation
function c=cross1(a,b)
c = [a(2,:).*b(3,:)-a(3,:).*b(2,:);
     a(3,:).*b(1,:)-a(1,:).*b(3,:);
     a(1,:).*b(2,:)-a(2,:).*b(1,:)];
end
