% MR and SR versus Ku in the case Kr=Ku
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
Kr_List=[2,3,4,5,6].';
minR_move=zeros(length(Kr_List),1);
minR_PhaseOnly=zeros(length(Kr_List),1);
minR_BenchMark=zeros(length(Kr_List),1);
minR_UpBound=zeros(length(Kr_List),1);
minR_Case1=zeros(length(Kr_List),1);
sumR_BenchMark=zeros(length(Kr_List),1);
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
        [~,A_r,A_r_tilde,~,h_ris2ue,H_ap2ris,h_ap2ue,gamma_avg_square]=GenChannel(POS_BS,POS_RIS,POS_UE,Kr,Ku,Nr,Nrx,Nry,Nbs,lambda,compensate_PL);
        %% RCG AO algorithm
        [MR_PhaseOnly,phi_PhaseOnly,~,MR_Move,phi_Move_opt,v_Move_opt,~]=RCG_Joint(maxite,radius_move,A_r,A_r_tilde,POS_BS,POS_RIS,POS_UE,Kr,Ku,Nr,Ps,Pn,Nrx,Nry,Nbs,lambda,compensate_PL);
        %% the upper bound of MR, and the MR in Case 1
        upbound=log2(gamma_avg_square*SNR*Kr*Nr^2*Nbs/Ku+1);
        Case1=log2(gamma_avg_square*SNR*Nr^2*Nbs/Ku+1);
%         [C_MIMO]=getC_MIMO(A_r_tilde, Nr, Kr, Ku, Nbs, SNR);
        %% display results
        fprintf('MR_PhaseOnly = %f, MR_Move = %f, UB = %f, Case 1 = %f\n',MR_PhaseOnly,MR_Move,upbound,Case1);
        %% BenchMark, note that the benchmark is slowest in the simulated schemes
        [MR_BenchMark,SR_BenchMark,~,~,~]=BenchMark(H_ap2ris, h_ris2ue, h_ap2ue, Nbs, Ku, N, Ps, Pn*ones(Ku,1), angle(phi_PhaseOnly));
        minR_move(Kr_idx) = minR_move(Kr_idx) + MR_Move/MC_times;
        minR_PhaseOnly(Kr_idx) = minR_PhaseOnly(Kr_idx) + MR_PhaseOnly/MC_times;
        minR_BenchMark(Kr_idx) = minR_BenchMark(Kr_idx) + MR_BenchMark/MC_times;
        minR_UpBound(Kr_idx) = minR_UpBound(Kr_idx) + upbound/MC_times;
        minR_Case1(Kr_idx) = minR_Case1(Kr_idx) + Case1/MC_times;
        sumR_BenchMark(Kr_idx)= sumR_BenchMark(Kr_idx) + SR_BenchMark/MC_times;
    end
    stop=1;
end
figure;
plot(Kr_List,minR_move,'linewidth',1.5,'Marker','o');hold on;
plot(Kr_List,minR_PhaseOnly,'linewidth',1.5,'Marker','v');hold on;
plot(Kr_List,minR_BenchMark,'linewidth',1.5,'Marker','s');hold on;
plot(Kr_List,minR_UpBound,'linewidth',1.5,'Marker','^');hold on;
plot(Kr_List,minR_Case1,'linewidth',1.5,'Marker','>');hold on;
xlabel('Kr=Ku');
ylabel('MR (bps/Hz)');
legend('Move','Phase-only','Benchmark','UB','Case 1');
stop=1;

figure;
plot(Kr_List,Kr_List.*minR_move,'linewidth',1.5,'Marker','o');hold on;
plot(Kr_List,Kr_List.*minR_PhaseOnly,'linewidth',1.5,'Marker','v');hold on;
plot(Kr_List,sumR_BenchMark,'linewidth',1.5,'Marker','s');hold on;
plot(Kr_List,Kr_List.*minR_UpBound,'linewidth',1.5,'Marker','^');hold on;
plot(Kr_List,Kr_List.*minR_Case1,'linewidth',1.5,'Marker','>');hold on;
xlabel('Kr=Ku');
ylabel('SR (bps/Hz)');
legend('Move','Phase-only','Benchmark','UB','Case 1');
stop=1;