function [minUserRate,SumRate,wk,theta,hk]=BenchMark(H_ap2ris, h_ris2ue, h_ap2ue, Nbs, Ku, N, Ps, Pn, theta_ini)
    %% Outer Iteration
    wk=zeros(Nbs,Ku)+1j*zeros(Nbs,Ku);
    hk=zeros(Nbs,Ku)+1j*zeros(Nbs,Ku);
    %% Initialization
%     theta=2*pi*rand(N,1);
%     for k=1:Ku
%         hk(:,k)=(exp(1j*theta.')*diag(h_ris2ue(:,k)')*H_ap2ris + h_ap2ue(:,k)')';
%         wk(:,k)=sqrt(Ps/Ku)*(hk(:,k)./norm(hk(:,k)));
%     end
    %% Initialization
    theta=theta_ini;
    for k=1:Ku
        hk(:,k)=(exp(1j*theta.')*diag(h_ris2ue(:,k)')*H_ap2ris + h_ap2ue(:,k)')';
    end
    wk=hk/(hk'*hk);
    wk=sqrt(Ps/trace(wk'*wk))*wk;
    %%
    [Ak,Bk,Ck,r_kl,xi_k,z_km]=calcul_AkBk(hk,wk,Ku,Pn);
    [UserRate]=calcul_Rate(hk,wk,Pn,Ku);
    minUserRate=min(UserRate);
    SumRate=sum(UserRate);
    rho2=20;
    idx=0;
    while(1)
        idx=idx+1;
        eps2=1e-3;
        xi_k=min(r_kl)*ones(Ku,1);
        xid_k=zeros(Ku,1);
        zd_km=zeros(Ku,Ku);
        [wk_new,hk_new,theta_new,gamma_out]=Inner_Iteration(...
            hk,wk,Ak,Bk,Ck,theta,xi_k,z_km,xid_k,zd_km,rho2,Nbs,Ku,N,Ps,h_ris2ue,h_ap2ue,H_ap2ris);

        [Ak,Bk,Ck,r_kl,xi_k,z_km]=calcul_AkBk(hk_new,wk_new,Ku,Pn);
        [UserRate_new]=calcul_Rate(hk_new,wk_new,Pn,Ku);
        minUserRate_new=min(UserRate_new);
        SumRate_new=sum(UserRate_new);
        if (minUserRate_new-minUserRate) < eps2
            break;
        end
        minUserRate=minUserRate_new;
        SumRate=SumRate_new;
        %% Updation
        hk=hk_new;
        wk=wk_new;
        theta=theta_new;
    end

end

function [UserRate]=calcul_Rate(hk,wk,Pn,Ku)
    B=hk'*wk;
    UserRate=zeros(Ku,1);
    for k=1:Ku
        UserRate(k)=log2( 1+abs(B(k,k))^2/( norm(B(k,:))^2-abs(B(k,k))^2+Pn(k) ) );
    end
end
function [Ak,Bk,Ck,r_kl,xi_k,z_km]=calcul_AkBk(hk,wk,Ku,Pn)
    p=[1 0];
    Ak=cell(Ku,1);
    Bk=cell(Ku,1);
    Ck=zeros(Ku,1);
    r_kl=zeros(Ku,1);
    z_km=zeros(Ku,Ku)+1j*zeros(Ku,Ku);
    sumk=zeros(Ku,1);
    for k=1:Ku
        for m=1:Ku
            sumk(k)=sumk(k)+abs(hk(:,k)'*wk(:,m))^2;
        end
    end
    for k=1:Ku
        Ak{k}=[1 wk(:,k)'*hk(:,k);hk(:,k)'*wk(:,k) sumk(k)+Pn(k)];
        Bk{k}=inv(Ak{k}) * p.' /(p/Ak{k}*p.')*p/Ak{k};
        for m=1:Ku
            z_km=hk(:,k)'*wk(:,m);
        end
    end
    for k=1:Ku
        Ck(k)= real(  log(p/Ak{k}*p.')+trace(Bk{k}*Ak{k})-Bk{k}(1,1)-Pn(k)*Bk{k}(2,2)  );
        r_kl(k)=(1/log(2))*( Ck(k)-2*real( Bk{k}(1,2)*hk(:,k)'*wk(:,k) )-Bk{k}(2,2)*sumk(k) );
    end
    xi_k=min(r_kl)*ones(Ku,1);
end
function [wk,hk,theta,gamma]=Inner_Iteration(...
    hk,wk,Ak,Bk,Ck,theta,xi_k,z_km,xid_k,zd_km,rho2,Nbs,Ku,N,Ps,h_ris2ue,h_ap2ue,H_ap2ris)

    eps1=1e-4;
    gamma=xi_k(1);
    while(1)
        % Step 1: Update gamma
        gamma_new=( sum(xi_k+xid_k)*rho2+1 )/(rho2*Ku);
        % Step 2: Update z_km,xi_k
        f_kb=zeros(Ku,1);
        xi_k_new=zeros(Ku,1);
        z_km_new=zeros(Ku,Ku);
        beta_lower=zeros(Ku,1);
        beta_upper=10*ones(Ku,1);
        beta_mid=zeros(Ku,1);
        const_beta=zeros(Ku,1);
        for k=1:Ku
            f_kb(k)=Ck(k)-2*real( Bk{k}(1,2)*( hk(:,k)'*wk(:,k)-zd_km(k,k) ) ) - Bk{k}(2,2)* norm(hk(:,k)'*wk-zd_km(k,:) )^2 ...
                - log(2)*(gamma_new-xid_k(k));
            if f_kb(k) >= 0
                xi_k_new(k) = gamma_new - xid_k(k);
                for m=1:Ku
                    z_km_new(k,m) = hk(:,k)'*wk(:,m) - zd_km(k,m);
                end
            else
                idx=0;
                while(1) % Bisection Search
                    idx=idx+1;
                    beta_mid(k)=(beta_upper(k)+beta_lower(k))/2;
                    xi_k_new(k) = gamma_new - xid_k(k) - log(2)*beta_mid(k)/2;
                    for m=1:Ku
                        if m==k
                            z_km_new(k,m) = (hk(:,k)'*wk(:,m) - zd_km(k,m) - beta_mid(k)*conj(Bk{k}(1,2))) / (1 + beta_mid(k)*Bk{k}(2,2));
                        else
                            z_km_new(k,m) = (hk(:,k)'*wk(:,m) - zd_km(k,m)) / (1 + beta_mid(k)*Bk{k}(2,2));
                        end
                    end
                    const_beta(k) = Ck(k)-2*real(Bk{k}(1,2)*z_km_new(k,k))-Bk{k}(2,2)*norm(z_km_new(k,:))^2 - log(2)*xi_k_new(k);
                    if const_beta(k) < 0
                        beta_lower(k)=beta_mid(k);
                    else
                        beta_upper(k)=beta_mid(k);
                    end
                    if (abs(const_beta(k))<1e-5) || (idx>1000)
                        if idx>1000
                            fprintf('Warning: Unexpected Iteration Times = 1000!\n');
                        end
                        break;
                    end
                end

            end
        end
        % Step 3: Update wk
        wk_new=zeros(Nbs,Ku)+1j*zeros(Nbs,Ku);
        alpha_lower=0;
        alpha_upper=1e4;
        idx=0;
        [U,D]=eig(hk'*hk);
        D=real(diag(D));
        while(1)
            idx=idx+1;
            if idx==1
                wk_new=hk/(hk'*hk)*(z_km_new+zd_km);
                tr=trace(wk_new'*wk_new); 
                if tr <= Ps
                    break;
                end
            end
            alpha_mid=(alpha_upper+alpha_lower)/2;
            for k=1:Ku
                D_tmp=diag( 1./(D+alpha_mid*ones(Ku,1)) );
                inv_tmp=(eye(Nbs)-hk*(U*D_tmp*U')*hk')/alpha_mid;
                wk_new(:,k)=inv_tmp *( hk*(z_km_new(:,k)+zd_km(:,k)) );
            end
            tr=trace(wk_new'*wk_new); 
            if tr < Ps
                alpha_upper = alpha_mid;
            else
                alpha_lower = alpha_mid;
            end
            if (abs(tr-Ps)<1e-5) || (idx>1000)
                if idx>1000
                    fprintf('Warning: Unexpected Iteration Times = 1000!\n');
                end
                break;
            end
        end
        % Step 4: Update Theta
        H_rk=cell(Ku,1);
        for k=1:Ku
            H_rk{k}=diag(h_ris2ue(:,k)')*H_ap2ris;
        end
        Tp1_km=z_km_new+zd_km-h_ap2ue'*wk_new;
        Tp2_km_n=cell(Ku,Ku);
        tmp2=cell(Ku,Ku);
        for k=1:Ku
            for m=1:Ku
                tmp2{k,m} = H_rk{k}*wk_new(:,m);
                vec_tmp = exp(1j*theta).* ( H_rk{k}*wk_new(:,m) );
                sum_tmp = sum(vec_tmp);
                Tp2_km_n{k,m}=zeros(N,1);
                for n=1:N
                    Tp2_km_n{k,m}(n)=sum_tmp-vec_tmp(n);
                end
            end
        end
        Tpn=zeros(N,1)+1j*zeros(N,1);
        for n=1:N
            for k=1:Ku
                for m=1:Ku
                     Tpn(n)=Tpn(n)+tmp2{k,m}(n)*conj( Tp1_km(k,m) - Tp2_km_n{k,m}(n) );
                end
            end
        end
        theta_new=-angle(Tpn);
        % Step 5: Update zd_km, xid_k
        zd_km_new = z_km_new - hk'*wk + zd_km;
        xid_k_new =xi_k_new-gamma_new+xid_k;
        %% Stop if not improve
        if (gamma_new - gamma) < eps1
            if gamma_new < gamma
                stop=1;
            end
            break;
        end
        %% Update Inner Iteration parameters
        gamma = gamma_new;
        xi_k = xi_k_new;
        z_km = z_km_new;
        wk = wk_new;
        theta = theta_new;
        xid_k = xid_k_new;
        zd_km = zd_km_new;
        hk_new = zeros(Nbs,Ku) + 1j*zeros(Nbs,Ku);
        Theta=diag(exp(1j*theta));
        for k=1:Ku
            hk_new(:,k)=(exp(1j*theta.')*diag(h_ris2ue(:,k)')*H_ap2ris + h_ap2ue(:,k)')';
        end
        hk = hk_new;

    end
    
end