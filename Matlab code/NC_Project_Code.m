close all
clear
clc
clearvars
yalmip clear

%% Model of the system

% Continuous time (given)
Atot = [0 1 0 3.2 1.98;
        0 0 1 -14.72 0.49;
       -8.86 8.50 9.39 -7.92 36.01;
        1.69 1.26 0.08 0 1;
       -7.52 -5.23 0.49 32.32 -1.36];
         
Bdec{1} = [0 0 1 0 0]';
Bdec{2} = [0 0 0 0 1]';
Btot = [Bdec{1}, Bdec{2}];
    
Ctot = eye(5);
Cdec{1} = Ctot(1:3,:);
Cdec{2} = Ctot(4:5,:);

sysCT = ss(Atot,Btot,Ctot,0); % State-space repres. of the CT system

n = 5;      % Number of states (dimension of matrix A)
N = 2;      % Number of subsystems

% Discrete time
h = 0.1;    % Sampling time

sysDT = c2d(sysCT,h); % State-space repres. of the DT system
Ftot = sysDT.A;
Gtot = sysDT.B;
Htot = sysDT.C;

% Ftot = expm(Atot*h);
% Gtot = Atot\(expm(Atot*h)-eye(size(Ftot,1)))*Btot;
% Htot = Ctot;

Gdec{1} = Gtot(:,1);
Gdec{2} = Gtot(:,1);

Hdec{1} = Htot(1:3,:);
Hdec{2} = Htot(4:5,:);

% Random initial conditions
% x0 = [0; 0; 1; 0; 0];
x0 = rand(5,1); % pick 5 random values between 0 and 1

Tmax = 5;
T = 0:0.01:Tmax;

%% Open-loop eigenvalues and the spectral abscissa

% Continuous time
eigA = eig(Atot)
disp(['The eigenvalue of the c.t. system with greater real part is: ', num2str(max(real(eig(Atot))))]);

% Discrete time
eigF = eig(Ftot)
disp(['The eigenvalue of the d.t. system with maximum abs values is: ', num2str(max(abs(eig(Ftot))))]);

%% Centralized controllers

% Centralized control structure
CentrStruct = ones(N,N);                                                

% CT and DT fixed modes
[fm_CT_centralized] = di_fixed_modes(Atot,Bdec,Cdec,N,CentrStruct,3);
[fm_DT_centralized] = di_fixed_modes(Ftot,Gdec,Hdec,N,CentrStruct,3);

% CT and DT STABILIZING state-feedback control laws
[K_CT_stabilizing_centralized, rho_CT_stabilizing_centralized, feas_CT_stabilizing_centralized] = LMI_CT_DeDicont(Atot,Bdec,Cdec,N,CentrStruct);
[K_DT_stabilizing_centralized, rho_DT_stabilizing_centralized, feas_DT_stabilizing_centralized] = LMI_DT_DeDicont(Ftot,Gdec,Hdec,N,CentrStruct);

sysCTtraj = ss((Atot+Btot*K_CT_stabilizing_centralized),zeros(5,2),Ctot,[]);
sysDTtraj = ss((Ftot+Gtot*K_DT_stabilizing_centralized),zeros(5,2),Htot,[],h);
x_traj_CT_stabilizing = initial(sysCTtraj,x0,T); % forced response to input signal
x_traj_DT_stabilizing = initial(sysDTtraj,x0,[0:h:T(end)]);


% Disk (discrete-time)

[K_DT_disk_centralized, rho_DT_disk_centralized, feas_DT_disk_centralized] = LMI_DT_DeDicont_Disk(Ftot,Gdec,Hdec,N,CentrStruct);
sysDTtraj = ss((Ftot+Gtot*K_DT_disk_centralized),zeros(5,2),Htot,[],h);
x_traj_DT_disk = initial(sysDTtraj,x0,[0:h:T(end)]);


% Disk (continuous-time)

[K_CT_disk_centralized, rho_CT_disk_centralized, feas_CT_disk_centralized] = LMI_CT_DeDicont_Disk(Atot,Bdec,Cdec,N,CentrStruct);
sysCTtraj = ss((Atot+Btot*K_CT_disk_centralized),zeros(5,2),Ctot,[]);
x_traj_CT_disk = initial(sysCTtraj,x0,T);


% Disk + Control Effort Minimization (continuous-time)

[K_CT_DCE_centralized, rho_CT_DCE_centralized, feas_CT_DCE_centralized] = LMI_CT_DeDicont_Disk_ControlEffort(Atot,Bdec,Cdec,N,CentrStruct);
sysCTtraj = ss((Atot+Btot*K_CT_DCE_centralized),zeros(5,2),Ctot,[]);
x_traj_CT_DCE = initial(sysCTtraj,x0,T);


% D-Region (continuous-time)

[K_CT_Dregion_centralized, rho_CT_Dregion_centralized, feas_CT_Dregion_centralized] = LMI_CT_DeDicont_Dregion(Atot,Bdec,Cdec,N,CentrStruct);
sysCTtraj = ss((Atot+Btot*K_CT_Dregion_centralized),zeros(5,2),Ctot,[]);
x_traj_CT_Dregion = initial(sysCTtraj,x0,T);


% H2 (continuous-time)

[K_CT_H2_centralized, rho_CT_H2_centralized, feas_CT_H2_centralized] = LMI_CT_DeDicont_H2(Atot,Bdec,Cdec,N,CentrStruct);
sysCTtraj = ss((Atot+Btot*K_CT_H2_centralized),zeros(5,2),Ctot,[]);
x_traj_CT_H2 = initial(sysCTtraj,x0,T);

% Disk_Center (discrete-time)

[K_DT_Disk_Center_centralized, rho_DT_Disk_Center_centralized, feas_DT_Disk_Center_centralized] = LMI_DT_DeDicont_Disk_Center(Ftot,Gdec,Hdec,N,CentrStruct,0.6);
sysDTtraj = ss((Ftot+Gtot*K_DT_Disk_Center_centralized),zeros(5,2),Htot,[],h);
x_traj_DT_Disk_Center = initial(sysDTtraj,x0,[0:h:T(end)]);


figure(1)
hold on
title('Centralized controllers (continuous time)')
grid on
plot(T,x_traj_CT_stabilizing(:,3))
plot(T,x_traj_CT_disk(:,3))
% plot(T,x_traj_CT_DCE(:,3))
plot(T,x_traj_CT_Dregion(:,3))
plot(T,x_traj_CT_H2(:,3))
legend('Stabilizing','Disk','D-Region','H2')
ylabel('x3')
xlabel('Time [s]')
hold off

figure(2)
hold on
title('Centralized controllers (discrete time)')
plot([0:h:Tmax],x_traj_DT_stabilizing(:,3),'.-')
plot([0:h:Tmax],x_traj_DT_Disk_Center(:,3),'.-')
plot([0:h:Tmax],x_traj_DT_disk(:,3),'.-')
grid on
legend('Stabilizing','Disk Centered','Disk')
ylabel('x3')
xlabel('Time [s]')
hold off

%% Decentralized controllers

% Decentralized control structure
CentrStruct=diag(ones(N,1));

% CT and DT fixed modes
[fm_CT_decentralized] = di_fixed_modes(Atot,Bdec,Cdec,N,CentrStruct,3);
[fm_DT_decentralized] = di_fixed_modes(Ftot,Gdec,Hdec,N,CentrStruct,3);


% Stabilizing (Continuous-time and discrete-time)

[K_CT_stabilizing_decentralized, rho_CT_stabilizing_decentralized, feas_CT_stabilizing_decentralized] = LMI_CT_DeDicont(Atot,Bdec,Cdec,N,CentrStruct);
[K_DT_stabilizing_decentralized, rho_DT_stabilizing_decentralized, feas_DT_stabilizing_decentralized] = LMI_DT_DeDicont(Ftot,Gdec,Hdec,N,CentrStruct);

sysCTtraj = ss((Atot+Btot*K_CT_stabilizing_decentralized),zeros(5,2),Ctot,[]);
sysDTtraj = ss((Ftot+Gtot*K_DT_stabilizing_decentralized),zeros(5,2),Htot,[],h);
x_traj_CT_stabilizing_dec = initial(sysCTtraj,x0,T);
x_traj_DT_stabilizing_dec = initial(sysDTtraj,x0,[0:h:T(end)]);


% Disk (discrete-time)

[K_DT_disk_decentralized, rho_DT_disk_decentralized, feas_DT_disk_decentralized] = LMI_DT_DeDicont_Disk(Ftot,Gdec,Hdec,N,CentrStruct);
sysDTtraj = ss((Ftot+Gtot*K_DT_disk_decentralized),zeros(5,2),Htot,[],h);
x_traj_DT_disk_dec = initial(sysDTtraj,x0,[0:h:T(end)]);


% Disk (continuous-time)

[K_CT_disk_decentralized, rho_CT_disk_decentralized, feas_CT_disk_decentralized] = LMI_CT_DeDicont_Disk(Atot,Bdec,Cdec,N,CentrStruct);
sysCTtraj = ss((Atot+Btot*K_CT_disk_decentralized),zeros(5,2),Ctot,[]);
x_traj_CT_disk_dec = initial(sysCTtraj,x0,T);


% D-Region (continuous-time)

[K_CT_Dregion_decentralized, rho_CT_Dregion_decentralized, feas_CT_Dregion_decentralized] = LMI_CT_DeDicont_Dregion(Atot,Bdec,Cdec,N,CentrStruct);
sysCTtraj = ss((Atot+Btot*K_CT_Dregion_decentralized),zeros(5,2),Ctot,[]);
x_traj_CT_Dregion_dec = initial(sysCTtraj,x0,T);


% H2 (continuous-time)

[K_CT_H2_decentralized, rho_CT_H2_decentralized, feas_CT_H2_decentralized] = LMI_CT_DeDicont_H2(Atot,Bdec,Cdec,N,CentrStruct);
sysCTtraj = ss((Atot+Btot*K_CT_H2_decentralized),zeros(5,2),Ctot,[]);
x_traj_CT_H2_dec = initial(sysCTtraj,x0,T);


% Disk_Center (discrete-time)

[K_DT_Disk_Center_decentralized, rho_DT_Disk_Center_decentralized, feas_DT_Disk_Center_decentralized] = LMI_DT_DeDicont_Disk_Center(Ftot,Gdec,Hdec,N,CentrStruct,0.45);
sysDTtraj = ss((Ftot+Gtot*K_DT_Disk_Center_decentralized),zeros(5,2),Htot,[],h);
x_traj_DT_Disk_Center_dec = initial(sysDTtraj,x0,[0:h:T(end)]);


figure(3)
hold on
title('Decentralized controllers (continuous time)')
grid on
plot(T,x_traj_CT_stabilizing_dec(:,3))
% plot(T,x_traj_CT_disk_dec(:,3)) --> UNFEASIBLE
plot(T,x_traj_CT_Dregion_dec(:,3))
plot(T,x_traj_CT_H2_dec(:,3))
legend('Stabilizing','D-Region','H2')
ylabel('x3')
xlabel('Time [s]')
hold off

figure(4) % ALL UNFEASIBLE IN DT SINCE THE SIMPLER STABILIZING CONTROL LAW IS UNFEASIBLE (DOESN'T EXIST A CONTROLLER)
hold on
title('Decentralized controllers (discrete time)')
plot([0:h:Tmax],x_traj_DT_stabilizing_dec(:,3),'.-')
plot([0:h:Tmax],x_traj_DT_Disk_Center_dec(:,3),'.-')
plot([0:h:Tmax],x_traj_DT_disk_dec(:,3),'.-')
grid on
legend('Stabilizing','Disk Centered','Disk')
ylabel('x3')
xlabel('Time [s]')
hold off

%% Distributed controller (Unidirectional String: 1-->2)

% %String
% CentrStruct = [1 0; 1 1];
% 
% % CT and DT fixed modes
% [fm_CT_distributed] = di_fixed_modes(Atot,Bdec,Cdec,N,CentrStruct,3);
% [fm_DT_distributed] = di_fixed_modes(Ftot,Gdec,Hdec,N,CentrStruct,3);
% 
% 
% % Stabilizing (Continuous-time and discrete-time)
% 
% [K_CT_stabilizing_distributed, rho_CT_stabilizing_distributed, feas_CT_stabilizing_distributed] = LMI_CT_DeDicont(Atot,Bdec,Cdec,N,CentrStruct);
% [K_DT_stabilizing_distributed, rho_DT_stabilizing_distributed, feas_DT_stabilizing_distributed] = LMI_DT_DeDicont(Ftot,Gdec,Hdec,N,CentrStruct);
% 
% sysCTtraj = ss((Atot+Btot*K_CT_stabilizing_distributed),zeros(5,2),Ctot,[]);
% sysDTtraj = ss((Ftot+Gtot*K_DT_stabilizing_distributed),zeros(5,2),Htot,[],h);
% x_traj_CT_stabilizing_distr = initial(sysCTtraj,x0,T);
% x_traj_DT_stabilizing_distr = initial(sysDTtraj,x0,[0:h:T(end)]);
% 
% 
% % Disk (discrete-time)
% 
% [K_DT_disk_distributed, rho_DT_disk_distributed, feas_DT_disk_distributed] = LMI_DT_DeDicont_Disk(Ftot,Gdec,Hdec,N,CentrStruct);
% sysDTtraj = ss((Ftot+Gtot*K_DT_disk_distributed),zeros(5,2),Htot,[],h);
% x_traj_DT_disk_distr = initial(sysDTtraj,x0,[0:h:T(end)]);
% 
% 
% % Disk (continuous-time)
% 
% [K_CT_disk_distributed, rho_CT_disk_distributed, feas_CT_disk_distributed] = LMI_CT_DeDicont_Disk(Atot,Bdec,Cdec,N,CentrStruct);
% sysCTtraj = ss((Atot+Btot*K_CT_disk_distributed),zeros(5,2),Ctot,[]);
% x_traj_CT_disk_distr = initial(sysCTtraj,x0,T);
% 
% 
% % D-Region (continuous-time)
% 
% [K_CT_Dregion_distributed, rho_CT_Dregion_distributed, feas_CT_Dregion_distributed] = LMI_CT_DeDicont_Dregion(Atot,Bdec,Cdec,N,CentrStruct);
% sysCTtraj = ss((Atot+Btot*K_CT_Dregion_distributed),zeros(5,2),Ctot,[]);
% x_traj_CT_Dregion_distr = initial(sysCTtraj,x0,T);
% 
% 
% % H2 (continuous-time)
% 
% [K_CT_H2_distributed, rho_CT_H2_distributed, feas_CT_H2_distributed] = LMI_CT_DeDicont_H2(Atot,Bdec,Cdec,N,CentrStruct);
% sysCTtraj = ss((Atot+Btot*K_CT_H2_distributed),zeros(5,2),Ctot,[]);
% x_traj_CT_H2_distr = initial(sysCTtraj,x0,T);
% 
% 
% % Disk_Center (discrete-time)
% 
% [K_DT_Disk_Center_distributed, rho_DT_Disk_Center_distributed, feas_DT_Disk_Center_distributed] = LMI_DT_DeDicont_Disk_Center(Ftot,Gdec,Hdec,N,CentrStruct,0.45);
% sysDTtraj = ss((Ftot+Gtot*K_DT_Disk_Center_distributed),zeros(5,2),Htot,[],h);
% x_traj_DT_Disk_Center_distr = initial(sysDTtraj,x0,[0:h:T(end)]);
% 
% 
% figure(5)
% hold on
% title('Distributed (unidirectional string) controllers (continuous time)')
% grid on
% plot(T,x_traj_CT_stabilizing_distr(:,3))
% % plot(T,x_traj_CT_disk_distr(:,3)) --> UNFEASIBLE
% plot(T,x_traj_CT_Dregion_distr(:,3))
% plot(T,x_traj_CT_H2_distr(:,3))
% legend('Stabilizing','D-Region','H2')
% ylabel('x3')
% xlabel('Time [s]')
% hold off
% 
% figure(6) % ALL UNFEASIBLE IN DT SINCE THE SIMPLER STABILIZING CONTROL LAW IS UNFEASIBLE (DOESN'T EXIST A CONTROLLER)
% hold on
% title('Distributed (unidirectional string) controllers (discrete time)')
% plot([0:h:Tmax],x_traj_DT_stabilizing_distr(:,3),'.-')
% plot([0:h:Tmax],x_traj_DT_Disk_Center_distr(:,3),'.-')
% plot([0:h:Tmax],x_traj_DT_disk_distr(:,3),'.-')
% grid on
% legend('Stabilizing','Disk Centered','Disk')
% ylabel('x3')
% xlabel('Time [s]')
% hold off

%% Distributed controller (Unidirectional String: 2-->1)

%String
CentrStruct = [1 1; 0 1];

% CT and DT fixed modes
[fm_CT_distributed] = di_fixed_modes(Atot,Bdec,Cdec,N,CentrStruct,3);
[fm_DT_distributed] = di_fixed_modes(Ftot,Gdec,Hdec,N,CentrStruct,3);


% Stabilizing (Continuous-time and discrete-time)

[K_CT_stabilizing_distributed, rho_CT_stabilizing_distributed, feas_CT_stabilizing_distributed] = LMI_CT_DeDicont(Atot,Bdec,Cdec,N,CentrStruct);
[K_DT_stabilizing_distributed, rho_DT_stabilizing_distributed, feas_DT_stabilizing_distributed] = LMI_DT_DeDicont(Ftot,Gdec,Hdec,N,CentrStruct);

sysCTtraj = ss((Atot+Btot*K_CT_stabilizing_distributed),zeros(5,2),Ctot,[]);
sysDTtraj = ss((Ftot+Gtot*K_DT_stabilizing_distributed),zeros(5,2),Htot,[],h);
x_traj_CT_stabilizing_distr = initial(sysCTtraj,x0,T);
x_traj_DT_stabilizing_distr = initial(sysDTtraj,x0,[0:h:T(end)]);


% Disk (discrete-time)

[K_DT_disk_distributed, rho_DT_disk_distributed, feas_DT_disk_distributed] = LMI_DT_DeDicont_Disk(Ftot,Gdec,Hdec,N,CentrStruct);
sysDTtraj = ss((Ftot+Gtot*K_DT_disk_distributed),zeros(5,2),Htot,[],h);
x_traj_DT_disk_distr = initial(sysDTtraj,x0,[0:h:T(end)]);


% Disk (continuous-time)

[K_CT_disk_distributed, rho_CT_disk_distributed, feas_CT_disk_distributed] = LMI_CT_DeDicont_Disk(Atot,Bdec,Cdec,N,CentrStruct);
sysCTtraj = ss((Atot+Btot*K_CT_disk_distributed),zeros(5,2),Ctot,[]);
x_traj_CT_disk_distr = initial(sysCTtraj,x0,T);


% D-Region (continuous-time)

[K_CT_Dregion_distributed, rho_CT_Dregion_distributed, feas_CT_Dregion_distributed] = LMI_CT_DeDicont_Dregion(Atot,Bdec,Cdec,N,CentrStruct);
sysCTtraj = ss((Atot+Btot*K_CT_Dregion_distributed),zeros(5,2),Ctot,[]);
x_traj_CT_Dregion_distr = initial(sysCTtraj,x0,T);


% H2 (continuous-time)

[K_CT_H2_distributed, rho_CT_H2_distributed, feas_CT_H2_distributed] = LMI_CT_DeDicont_H2(Atot,Bdec,Cdec,N,CentrStruct);
sysCTtraj = ss((Atot+Btot*K_CT_H2_distributed),zeros(5,2),Ctot,[]);
x_traj_CT_H2_distr = initial(sysCTtraj,x0,T);


% Disk_Center (discrete-time)

[K_DT_Disk_Center_distributed, rho_DT_Disk_Center_distributed, feas_DT_Disk_Center_distributed] = LMI_DT_DeDicont_Disk_Center(Ftot,Gdec,Hdec,N,CentrStruct,0.45);
sysDTtraj = ss((Ftot+Gtot*K_DT_Disk_Center_distributed),zeros(5,2),Htot,[],h);
x_traj_DT_Disk_Center_distr = initial(sysDTtraj,x0,[0:h:T(end)]);


figure(7)
hold on
title('Distributed (unidirectional string) controllers (continuous time)')
grid on
plot(T,x_traj_CT_stabilizing_distr(:,3))
% plot(T,x_traj_CT_disk_distr(:,3)) --> UNFEASIBLE
plot(T,x_traj_CT_Dregion_distr(:,3))
plot(T,x_traj_CT_H2_distr(:,3))
legend('Stabilizing','D-Region','H2')
ylabel('x3')
xlabel('Time [s]')
hold off

figure(8) % ALL UNFEASIBLE IN DT SINCE THE SIMPLER STABILIZING CONTROL LAW IS UNFEASIBLE (DOESN'T EXIST A CONTROLLER)
hold on
title('Distributed (unidirectional string) controllers (discrete time)')
plot([0:h:Tmax],x_traj_DT_stabilizing_distr(:,3),'.-')
plot([0:h:Tmax],x_traj_DT_Disk_Center_distr(:,3),'.-')
plot([0:h:Tmax],x_traj_DT_disk_distr(:,3),'.-')
grid on
legend('Stabilizing','Disk Centered','Disk')
ylabel('x3')
xlabel('Time [s]')
hold off








