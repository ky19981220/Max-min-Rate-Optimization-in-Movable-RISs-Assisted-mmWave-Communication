function [A_b,A_r,A_r_tilde,A_k,h_ris2ue,H_ap2ris,h_ap2ue,gamma_avg_square]=GenChannel(POS_BS,POS_RIS,POS_UE,Kr,Ku,Nr,Nrx,Nry,Nbs,lambda,compensate_PL)
    N=Kr*Nr;
    a_b=cell(Kr,1);
    for r=1:Kr
        a_b{r}=SteerVectorPlanar(sqrt(Nbs),sqrt(Nbs),POS_BS,POS_RIS(r,:),lambda,1);
    end
    A_b=[];
    for r=1:Kr
        A_b=[A_b,a_b{r}];
    end
    look=A_b'*A_b;
    a_r=cell(Kr,1);
    for r=1:Kr
        a_r{r}=SteerVectorPlanar(Nrx,Nry,POS_BS,POS_RIS(r,:),lambda,1);
    end
    a_rk=cell(Kr,Ku);
    for r=1:Kr
        for u=1:Ku
            a_rk{r,u}=SteerVectorPlanar(Nrx,Nry,POS_RIS(r,:),POS_UE(u,:),lambda,1);
        end
    end
    A_ri=cell(Kr,1);
    A_r=[];
    for r=1:Kr
        for u=1:Ku
            A_ri{r}=[A_ri{r},a_rk{r,u}];
        end
        A_r=[A_r;A_ri{r}];
    end
    % Path Loss
    D_gi=cell(Kr,1);
    for r=1:Kr
        D_gi{r}=norm(POS_RIS(r,:)-POS_BS);
    end
    D_qik=cell(Kr,Ku);
    V_ik=zeros(Kr,Ku)+1j*zeros(Kr,Ku);
    for r=1:Kr
        for u=1:Ku
            D_qik{r,u}=norm(POS_RIS(r,:)-POS_UE(u,:));
            V_ik(r,u)=sqrt(compensate_PL)*lambda^2*exp(1j*2*pi/lambda*(D_gi{r}+D_qik{r,u}) )/(16*pi^2*D_gi{r}*D_qik{r,u});
        end
    end
    gamma_avg_square=sum(sum(abs(V_ik).^2) )/Kr/Ku;
    Q=V_ik;
    A_r_tilde=kron(Q,ones(Nr,1)).*A_r;
    A_ri_tilde=cell(Kr,1);
    for r=1:Kr
        A_ri_tilde{r}=A_r_tilde( (r-1)*Nr+1:r*Nr,: );
    end
    h_ris2ue=zeros(N,Ku)+1j*zeros(N,Ku);
    for u=1:Ku
        h_ris2ue_k=[];
        for r=1:Kr
            h_ris2ue_k=[h_ris2ue_k;sqrt(compensate_PL)*lambda^2/(16*pi^2*D_gi{r}*D_qik{r,u})*exp(1j*2*pi/lambda*D_qik{r,u})*a_rk{r,u}];
        end
        h_ris2ue(:,u)=h_ris2ue_k;
    end
    H_ap2ris=[];
    for r=1:Kr
        H_ap2ris=[H_ap2ris;exp(-1j*2*pi/lambda*D_gi{r})*(ones(Nr,1)*a_b{r}')];
    end
    %%
    A_k=cell(Ku,1);
    for u=1:Ku
        A_k{u}=[];
        for r=1:Kr
            A_k{u}=[A_k{u};V_ik(r,u)*(conj(a_rk{r,u}).*ones(Nr,1))*a_b{r}'];
        end
    end
    h_ap2ue=zeros(Nbs,Ku);
end