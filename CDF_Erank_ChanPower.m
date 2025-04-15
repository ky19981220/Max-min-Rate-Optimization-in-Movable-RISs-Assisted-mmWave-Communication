% CDF of configured channel's effective rank and power in the case Kr=Ku=6
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
Kr_List=[6];
ER_BenchMark=zeros(MC_times,1);
ER_orig=zeros(MC_times,1);
ER_move=zeros(MC_times,1);
trH_BenchMark=zeros(MC_times,1);
trH_orig=zeros(MC_times,1);
trH_move=zeros(MC_times,1);
trH_UB=zeros(MC_times,1);
trH_LB=zeros(MC_times,1);
SR_noMove=zeros(MC_times,1);
SR_move=zeros(MC_times,1);
MR_upbound=zeros(MC_times,1);
MR_case1=zeros(MC_times,1);
MR_BenchMark=zeros(MC_times,1);
SR_BenchMark=zeros(MC_times,1);
for Kr_idx=1:length(Kr_List)
    Kr=Kr_List(Kr_idx);
    Ku=Kr;
    POS_RIS=POS_RIS_all(1:Ku,:);
    Ps=1e0;
    SNR=Ps/Pn;
    [Nr,Nrx,Nry]=RIS_UPAMapping(Kr); % get Nr,Nrx,Nry for given N and Kr
    for ite=1:MC_times
        POS_UE=POS_UE_mc{ite}(1:Ku,:);
        % generate channel with random RISs' phases, where
        % "A_b, A_r, A_r_tilde, gamma_avg_square" are independent of RISS' phases.
        [A_b,A_r,A_r_tilde,~,h_ris2ue,H_ap2ris,h_ap2ue,gamma_avg_square]=GenChannel(POS_BS,POS_RIS,POS_UE,Kr,Ku,Nr,Nrx,Nry,Nbs,lambda,compensate_PL);
        %% RCG AO algorithm
        [MR_PhaseOnly,phi_PhaseOnly,H_PhaseOnly,MR_Move,phi_Move_opt,v_Move_opt,H_Move]=RCG_Joint(maxite,radius_move,A_r,A_r_tilde,POS_BS,POS_RIS,POS_UE,Kr,Ku,Nr,Ps,Pn,Nrx,Nry,Nbs,lambda,compensate_PL);
        [ER_orig(ite),trH_orig(ite)]=Effec_Rank(H_PhaseOnly);
        [ER_move(ite),trH_move(ite)]=Effec_Rank(H_Move);
        %% the upper bound of MR, and the MR in Case 1
        trH_UB(ite)=Nbs*N^2*gamma_avg_square;
        trH_LB(ite)=Nbs*N^2*gamma_avg_square/Ku;
        %% BenchMark, note that the benchmark is slowest in the simulated schemes
        [MR_BenchMark,SR_BenchMark,F_BF,angles_phi,H_BenchMark]=BenchMark(H_ap2ris, h_ris2ue, h_ap2ue, Nbs, Ku, N, Ps, Pn*ones(Ku,1), angle(phi_PhaseOnly));
        [ER_BenchMark(ite),trH_BenchMark(ite)]=Effec_Rank(H_BenchMark);
    end
    stop=1;
end

x_list=1:0.2:Kr;
f_ER_BenchMark=zeros(length(x_list),1);
f_ER_orig=zeros(length(x_list),1);
f_ER_move=zeros(length(x_list),1);
for x_idx=1:length(x_list)
    x=x_list(x_idx);
    f_ER_BenchMark(x_idx)=sum(ER_BenchMark<x)/MC_times;
    f_ER_orig(x_idx)=sum(ER_orig<x)/MC_times;
    f_ER_move(x_idx)=sum(ER_move<x)/MC_times;
end
figure;
plot(x_list,f_ER_move,'linewidth',1.5);hold on;
plot(x_list,f_ER_orig,'linewidth',1.5);hold on;
plot(x_list,f_ER_BenchMark,'linewidth',1.5);hold on;
xlabel('Channel effective rank');
ylabel('Empirical CDF');
title('Pn=-100 dBm, Ps=50 dBm, Kr=Ku=6, N=240');
legend('Move','Phase-only','BenchMark');
grid on;
stop=1;

x_list=0:20:500;
f_ER_UB=zeros(length(x_list),1);
f_ER_LB=zeros(length(x_list),1);
for x_idx=1:length(x_list)
    x=x_list(x_idx);
    f_ER_BenchMark(x_idx)=sum(trH_BenchMark<x)/MC_times;
    f_ER_orig(x_idx)=sum(trH_orig<x)/MC_times;
    f_ER_move(x_idx)=sum(trH_move<x)/MC_times;
    f_ER_UB(x_idx)=sum(trH_UB<x)/MC_times;
    f_ER_LB(x_idx)=sum(trH_LB<x)/MC_times;
end
figure;
plot(x_list,f_ER_move,'linewidth',1.5);hold on;
plot(x_list,f_ER_orig,'linewidth',1.5);hold on;
plot(x_list,f_ER_BenchMark,'linewidth',1.5);hold on;
plot(x_list,f_ER_UB,'linewidth',1.5);hold on;
plot(x_list,f_ER_LB,'linewidth',1.5);hold on;
xlabel('Channel power');
ylabel('Empirical CDF');
title('Pn=-100 dBm, Ps=50 dBm, Kr=Ku=6, N=240');
legend('Move','Phase-only','BenchMark','UB','Case 1');
grid on;

function [ER,trH]=Effec_Rank(H)
    trH=trace(H'*H);
    [~,S,~]=svd(H);
    S=abs(diag(S));
    nor_S=S./sum(S);
    ER=exp(-sum(nor_S.*log(nor_S)));
end