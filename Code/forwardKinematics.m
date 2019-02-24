function [T]=forwardKinematics(q,TefTool,alfa,d,a,n)
%% Calculates the forward kinematics of robot defined by the modefied DH parameters

%% Areguments:
% q: 1xn vector of joints angles of the robot
% DH PARAMETERS FOR THE ROBOT:
% alfa: 1xn vector
% d: 1xn vector
% a: 1xn vector
% TefTool: 4x4 transformation matrix from EEF frame to last frame of the
% robot
% n: an integer of the number of degrees of freedom of the robot

%% Returning value:
% T: is 4x4x(n+1) transformation matrices of the joints and the EEF.

% Copyright: Mohammad SAFEEA

%% Calculating the forward Kinematics
n_plus_1=n+1;
T=zeros(4,4,n_plus_1);
i=1;
T(:,:,i)=getDHMatrix(alfa(i),q(i),d(i),a(i));
    for i=2:n
        T(:,:,i)=T(:,:,i-1)*getDHMatrix(alfa(i),q(i),d(i),a(i));
        T(:,:,i)=normalizeColumns(T(:,:,i));
    end
        T(:,:,n_plus_1)=T(:,:,n)*TefTool;
        T(:,:,n_plus_1)=normalizeColumns(T(:,:,n_plus_1));
end

function T=getDHMatrix(alfa,theta,d,a)
T=zeros(4,4);

calpha=cos(alfa);
sinalpha=sin(alfa);
coshteta=cos(theta);
sintheta=sin(theta);

T(1,1)=coshteta;
T(2,1)=sintheta*calpha;
T(3,1)=sintheta*sinalpha;				

T(1,2)=-sintheta;
T(2,2)=coshteta*calpha;
T(3,2)=coshteta*sinalpha;

T(2,3)=-sinalpha;
T(3,3)=calpha;

T(1,4)=a;
T(2,4)=-sinalpha*d;
T(3,4)=calpha*d;
T(4,4)=1;

end

function normalizedT=normalizeColumns(T)
%% This function is used to normalize the columns of a rotation matrix with 
% some numerical errors resulting from matrix multiplication problems
r=zeros(4,3); % corrected rotatio matrix, with zero badding row at the end
for j=1:3
    r(1:3,j)=T(1:3,j)/norm(T(1:3,j));
end
normalizedT=[r,T(:,4)];
end
