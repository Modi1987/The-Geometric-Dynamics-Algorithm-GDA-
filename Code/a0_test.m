% Copyright: Mohammad SAFEEA, 24th-Feb-2019
% Script on the implementation of GDA dynamics algorithms.
n=8;
%% DH parameters
a=rand(n,1);
alfa=rand(n,1);
d=rand(n,1);
%% State of robot
q=rand(n,1);
dq=rand(n,1);
%% Inertial parameters
mcii=abs(rand(n,1));
pcii=rand(3,n);
Icii=zeros(3,3,n);
for i=1:n
    eigenVals=diag(abs(rand(3,1)));
    R=rotx(rand())*roty(rand())*rotz(rand());
    Icii(:,:,i)=R*eigenVals*R';
end
%% Calculate transofrmation matrices
TefTool=eye(4);
[T]=forwardKinematics(q,TefTool,alfa,d,a,n);
%% Calculate dynamic matrices
[Cent]=GetCentrifugalMatrixGDA(T,pcii,Icii,mcii);

[Bt]=GetCoriolisMatrixGDA(T,pcii,Icii,mcii,dq);

[At]=GetMassMatrixGDA(T,pcii,Icii,mcii);

[dM]=GetDerivativeOfMassMatrixGDA(T,pcii,Icii,mcii,dq);

[dM,Btdq]=GetCTdqGDA(T,pcii,Icii,mcii,dq);

%% Comparison of mass matrix calculated using GDA and GDAHJ
[M_GDAHJ]=GetInertiaMatrixGDAHJ(T,pcii,Icii,mcii);
total_error=sum(sum(abs(At-M_GDAHJ)));
disp('Total error in calculating Mass matrix using GDA and GDAHJ');
disp(total_error);
%% Numerical calculation of derivative of mass matrix
dt=0.000001;
q2=q+dq*dt;
T2=forwardKinematics(q2,TefTool,alfa,d,a,n);
dAt_dt=(GetMassMatrixGDA(T2,pcii,Icii,mcii)-At)/dt;
% Comparison with analytical solution
dM-dAt_dt
