function [K,rho,feas]=LMI_CT_DeDicont_H2(A,B,C,N,ContStruc)
% Computes, using LMIs, the distributed "state feedback" control law for the continuous-time system, with reference to the control
% information structure specified by 'ContStruc'.
%
% Inputs:
% - A: system matrix.
% - B: input matrices (i.e., B{1},..., B{N} are the input matrices of the decomposed system, one for each channel).
% - C: output matrices  (i.e., C{1},..., C{N} are the output matrices of the decomposed system, one for each channel, where [Cdec{1}',...,
% Cdec{N}']=I).
% - N: number of subsystems.
% - ContStruc: NxN matrix that specifies the information structure
% constraints (ContStruc(i,j)=1 if communication is allowed between channel
% j to channel i, ContStruc(i,j)=0 otherwise).
%
% Output:
% - K: structured control gain
% - rho: spectral abscissa of matrix (A+B*K) - note that [C{1}',...,
% C{N}']=I
% - feas: feasibility of the LMI problem (=0 if yes)

Btot=[];
for i=1:N
    m(i)=size(B{i},2);
    n(i)=size(C{i},1);
    Btot=[Btot,B{i}];
end
ntot=size(A,1);
mtot=sum(m);

yalmip clear

if ContStruc==ones(N,N)
    % Centralized design
    Y=sdpvar(ntot);
    L=sdpvar(mtot,ntot);
else
    % Decentralized/distributed design
    Y=[];
    L=sdpvar(mtot,ntot);
    minc=0;
    for i=1:N
        Y=blkdiag(Y,sdpvar(n(i)));
        ninc=0;
        for j=1:N
            if ContStruc(i,j)==0
                L(minc+1:minc+m(i),ninc+1:ninc+n(j))=zeros(m(i),n(j));
            end
            ninc=ninc+n(j);
        end
        minc=minc+m(i);
    end  
end

S=sdpvar(ntot+1,ntot+1);

Q = eye(ntot)*10;
R = eye(2)*0.01;

C = [sqrt(Q); 0,0,0,0,0];
D = [0,0; 0,0; 0,0; 0,0; sqrt(R)];

Bw = eye(ntot);

LMIconstr=[Y*A'+A*Y+Btot*L+L'*Btot'+Bw*Bw']<=1e-2*eye(ntot);

LMIconstr2 = [S C*Y+D*L; L'*D'+Y*C' Y]>=1e-2*eye(ntot*2+1);

Yconstr = Y>=1e-2*eye(ntot);

options=sdpsettings('solver','sedumi');
J=optimize([LMIconstr,LMIconstr2, Yconstr],[trace(S)],options);
feas=J.problem;
L=double(L);
Y=double(Y);

K=L/Y;
rho=max(real(eig(A+Btot*K)));
