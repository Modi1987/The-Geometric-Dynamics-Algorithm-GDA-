%% GDA geometric dynamics algorithm
%
%  Authors: Mohammad Safeea
%  University of Coimbra, Coimbra, Portugal
%  Ensam, ParisTech, Lille, France
%  2015-11-09
%  ------------
% 
%  This function is used to calculate the time derivative of the joint
%  space inertia matrix
%  GDA algorithm was used to perform the calculations.
%  modified Denavit Hartenberg convention was utilized

function [dM]=GetDerivativeOfMassMatrixGDA(T,Pcii,Icii,mcii,dq)
% dM: is the time derivative of joint space inertia matrix
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
dAi=zeros(3,n,n);midCi=zeros(3,n,n);
Ici=zeros(3,3,n);
Pcii_A=zeros(3,n);
Skew=zeros(3,3);
wj_corss_kj=zeros(3,n);
w=zeros(3,n);
vj=zeros(3,n);
vci=zeros(3,n);
w(:,1)=T(1:3,3,1)*dq(1);
wj_corss_kj(:,1)=zeros(3,1);
vj(:,1)=zeros(3,1);
Pcii_A(:,1)=T(1:3,1:3,1)*Pcii(:,1);
vci(:,1)=vj(:,1)+cross1(w(:,1),Pcii_A(:,1));
for j=2:n
    Pcii_A(:,j)=T(1:3,1:3,j)*Pcii(:,j);
    w(:,j)=w(:,j-1)+T(1:3,3,j)*dq(j);
    wj_corss_kj(:,j)=cross1(w(:,j),T(1:3,3,j));
    vj(:,j)=vj(:,j-1)+cross1(w(:,j-1),T(1:3,4,j)-T(1:3,4,j-1));
    vci(:,j)=vj(:,j)+cross1(w(:,j),Pcii_A(:,j));
end

for i=1:n
    Pci=Pcii_A(:,i)+T(1:3,4,i);
    Ici(:,:,i)=T(1:3,1:3,i)*Icii(1:3,1:3,i)*T(1:3,1:3,i)';
    Skew(1,2)=-w(3,i);Skew(1,3)=w(2,i);
    Skew(2,1)=w(3,i);Skew(2,3)=-w(1,i);
    Skew(3,1)=-w(2,i);Skew(3,2)=w(1,i);
    SkewIci=Skew*Ici(:,:,i);
    SymmIci=SkewIci+SkewIci'; 
    for j=1:i
        dAi(:,j,i)=SymmIci*T(1:3,3,j)+Ici(:,:,i)*wj_corss_kj(:,j);
         Pcij=Pci-T(1:3,4,j);
        Vcij=vci(:,i)-vj(:,j);
        midCi(:,j,i)=mcii(i)*(cross1(wj_corss_kj(:,j),Pcij)+...
        cross1(T(1:3,3,j),Vcij));
    end
end


[miCi,dM1]=GetFirstTermDerivative(T,mcii,Pcii_A,Ici,wj_corss_kj);

dM2=zeros(n,n);
dM2_vec=zeros(3,n);
midCi_acc=zeros(3,n);
miCi_acc=zeros(3,n);
PcimidCi_acc=zeros(3,n);
VcimiCi_acc=zeros(3,n);
PjmidCi=zeros(3,n);
VjmiCi=zeros(3,n);
dAi_acc=zeros(3,n);

for j=n:-1:1
    Pci=Pcii_A(:,j)+T(1:3,4,j);
    midCi_acc(:,1:j)=midCi_acc(:,1:j)+midCi(:,1:j,j);
    miCi_acc(:,1:j)=miCi_acc(:,1:j)+miCi(:,1:j,j);
    dAi_acc(:,1:j)=dAi_acc(:,1:j)+dAi(:,1:j,j);
    for k=1:j
        PcimidCi_acc(:,k)=PcimidCi_acc(:,k)+cross1(Pci,midCi(:,k,j));
        VcimiCi_acc(:,k)=VcimiCi_acc(:,k)+cross1(vci(:,j),miCi(:,k,j));
        PjmidCi(:,k)=-cross1(T(1:3,4,j),midCi_acc(:,k));
        VjmiCi(:,k)=-cross1(vj(:,j),miCi_acc(:,k));        
    end   
        dM2_vec(:,1:j)=dAi_acc(:,1:j)+VcimiCi_acc(:,1:j)+VjmiCi(:,1:j)...
        +PcimidCi_acc(:,1:j)+PjmidCi(:,1:j);
        dM2(j,1:j)=T(1:3,3,j)'*dM2_vec(:,1:j);
end

dM=dM1+dM2;
for i=1:n
    for j=1:i
        dM(j,i)=dM(i,j);
    end
end
end

%% Function that calculates miCi and dM1=(dK'*sigma(Ai+mi*Pcik*Ci))
function [miCi,dM1]=GetFirstTermDerivative(T,mcii,Pcii_A,Ici,wj_corss_kj)
n=max(size(mcii));
Ai=zeros(3,n,n);miCi=zeros(3,n,n);
temp=zeros(3,3);
for i=1:n
    Pci=Pcii_A(:,i)+T(1:3,4,i);
    for j=1:i
    Ai(:,j,i)=Ici(:,:,i)*T(1:3,3,j);
    Pcij=Pci-T(1:3,4,j);
    miCi(:,j,i)=mcii(i)*cross1(T(1:3,3,j),Pcij);
    end
end

At=zeros(n,n);Bt=zeros(n,1);
Fac_C=zeros(3,n);
Mac_A=zeros(3,n);
Pjp1_j=zeros(3,1);

start=n-1;
j=n;
for k=1:n 
    Mac_A(:,k)=Ai(:,k,j)+cross1(Pcii_A(:,j),miCi(:,k,j));            
end
    Fac_C=miCi(:,:,j);
        
At(j,:)=wj_corss_kj(:,j)'*Mac_A;
    
for j=start:-1:2 
    Pjp1_j=T(1:3,4,j+1)-T(1:3,4,j);
        for k=1:j
            Mac_A(:,k)=Mac_A(:,k)+Ai(:,k,j)+cross1(Pjp1_j,Fac_C(:,k))+cross1(Pcii_A(:,j),miCi(:,k,j));            
        end
    Fac_C(:,1:j)=Fac_C(:,1:j)+miCi(:,1:j,j);        
    At(j,1:j)=wj_corss_kj(:,j)'*Mac_A(:,1:j);
end
dM1=At;
end

%% Cross product calculation
function c=cross1(a,b)
c = [a(2,:).*b(3,:)-a(3,:).*b(2,:);
     a(3,:).*b(1,:)-a(1,:).*b(3,:);
     a(1,:).*b(2,:)-a(2,:).*b(1,:)];
end

