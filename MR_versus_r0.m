% MR and SR versus Ku in the case Kr=Ku=6
clear all;
figure(1);
addpath([pwd,'\func']);
MC_times=10; % Monte Carlo times
%% Default Parameters
fc=60e9;
lambda=3e8/fc;
maxite=4;
Nbs=100;
N=240;
distance_bs_ris=10;
distance_bs_ue=50;
radius_ue=10;
h_UE=1.5;
h_BS=10;
radius_move=0.1;
SNR_dB=150;
%% positions
POS_BS=[0,0,h_BS];
[POS_RIS_all]=getRISPOS(POS_BS,distance_bs_ris,sqrt(Nbs),sqrt(Nbs));
% UEs distributed within a circle region
maxKu=6;
center=[0,-distance_bs_ue,0];
[POS_UE_mc]=getUEPOS(MC_times,maxKu,radius_ue,center,h_UE);
%% Calculate path loss (PL), compensate PL using SNR
PathLoss=(lambda^2/(16*pi^2*distance_bs_ris*distance_bs_ue))^2;
compensate_PL=10^(SNR_dB/10);
SNR_dB_virtual=-10*log10(compensate_PL)+SNR_dB;
Pn=10^(-SNR_dB_virtual/10);
%% Monte Carlo (MC) simulation
radius_move_List=[0.01,0.02,0.04,0.08,0.16,0.32];
minR_move=zeros(length(radius_move_List),1);
minR_noMove=zeros(length(radius_move_List),1);
minR_BenchMark=zeros(length(radius_move_List),1);
minR_UpBound=zeros(length(radius_move_List),1);
Kr=6;
Ku=Kr;
for r_idx=1:length(radius_move_List)
    radius_move=radius_move_List(r_idx);
    POS_RIS=POS_RIS_all(1:Ku,:);
    Ps=1e0;
    SNR=Ps/Pn;
    [Nr,Nrx,Nry]=RIS_UPAMapping(Kr); % get Nr,Nrx,Nry for given N and Kr
    for ite=1:MC_times
        POS_UE=POS_UE_mc{ite}(1:Ku,:);
        % generate channel with random RISs' phases, where
        % "A_b, A_r, A_r_tilde, gamma_avg_square" are independent of RISS' phases.
        [~,A_r,A_r_tilde,~,h_ris2ue,H_ap2ris,h_ap2ue,gamma_avg_square]=GenChannel(POS_BS,POS_RIS,POS_UE,Kr,Ku,Nr,Nrx,Nry,Nbs,lambda,compensate_PL);
        %% RCG AO algorithm
        [MR_PhaseOnly,phi_PhaseOnly,~,MR_Move,phi_Move_opt,v_Move_opt,~]=RCG_Joint(maxite,radius_move,A_r,A_r_tilde,POS_BS,POS_RIS,POS_UE,Kr,Ku,Nr,Ps,Pn,Nrx,Nry,Nbs,lambda,compensate_PL);
        %% the upper bound of MR, and the MR in Case 1
        upbound=log2(gamma_avg_square*SNR*Kr*Nr^2*Nbs/Ku+1);
        Case1=log2(gamma_avg_square*SNR*Nr^2*Nbs/Ku+1);
        minR_move(r_idx) = minR_move(r_idx) + MR_Move/MC_times;
    end
    stop=1;
end
figure;
plot(radius_move_List,minR_move,'linewidth',1.5,'Marker','o');hold on;
xlabel('r_{0}');
ylabel('MR');
title('Pn=-100 dBm, Ps=50 dBm, N=240, Kr=Ku=6');
legend('Move');
grid on;