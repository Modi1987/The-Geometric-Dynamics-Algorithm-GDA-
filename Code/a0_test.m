% Copyright: Mohammad SAFEEA, 24th-Feb-2019
% Script on the implementation of GDA dynamics algorithms.
n=8;
%% DH parameters, (generated randomly)
a=rand(n,1);
alfa=rand(n,1);
d=rand(n,1);
%% State of robot, (generated randomly)
q=rand(n,1);
dq=rand(n,1);
%% Inertial parameters, (generated randomly)
mcii=abs(rand(n,1));
pcii=rand(3,n);
Icii=zeros(3,3,n); % must be symmetric positive-definite
for i=1:n
    eigenVals=diag(abs(rand(3,1)));
    R=rotx(rand())*roty(rand())*rotz(rand());
    Icii(:,:,i)=R*eigenVals*R';
end
%% Calculate transofrmation matrices
TefTool=eye(4);
[T]=forwardKinematics(q,TefTool,alfa,d,a,n);
%% Calculate dynamic matrices
[Cent]=GetCentrifugalMatrixGDA(T,pcii,Icii,mcii); % Centrifugal matrix

[Bt]=GetCoriolisMatrixGDA(T,pcii,Icii,mcii,dq); % Coriolis matrix

[At]=GetMassMatrixGDA(T,pcii,Icii,mcii); % Mass matrix

[dM]=GetDerivativeOfMassMatrixGDA(T,pcii,Icii,mcii,dq); % Time derivative of mass matrix

[dM,Btdq]=GetCTdqGDA(T,pcii,Icii,mcii,dq);

%% Comparison of mass matrix calculated using GDA and GDAHJ
[M_GDAHJ]=GetInertiaMatrixGDAHJ(T,pcii,Icii,mcii);
total_error=sum(sum(abs(At-M_GDAHJ)));
disp('Total error in calculating Mass matrix using GDA and GDAHJ');
disp(total_error);
%% Numerical calculation of time derivative of mass matrix
dt=0.000001;
q2=q+dq*dt;
T2=forwardKinematics(q2,TefTool,alfa,d,a,n);
dAt_dt=(GetMassMatrixGDA(T2,pcii,Icii,mcii)-At)/dt;
% Comparison with analytical solution
disp('The error in calculating the time derivative of the mass matrix using:')
disp('1- GDA (analytical solution)')
disp('2- Numerical differentiation');
disp('is ===>')
dM-dAt_dt
