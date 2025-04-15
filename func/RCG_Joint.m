function [MR_PhaseOnly,phi_PhaseOnly,H_PhaseOnly,MR_Move,phi_Move_opt,v_Move_opt,H_Move]=RCG_Joint(maxite,radius_move,A_r,A_r_tilde,POS_BS,POS_RIS,POS_UE,Kr,Ku,Nr,Ps,Pn,Nrx,Nry,Nbs,lambda,compensate_PL)
    % Joint phase and displacement optimization
    %% Inputs:
    % maxite:           maximum iteration times of algorithm 1
    % radius_move:      r_0
    % ...
    % compensate_PL:    compensate factor for path loss
    %% Outputs: 
    % MR_PhaseOnly:     optimal MR of Phase-only scheme  
    % phi_PhaseOnly:    optimized phase of Phase-only scheme  
    % H_PhaseOnly:      resulting channel of Phase-only scheme  
    % MR_Move:          optimal MR of "Move" scheme  
    % phi_Move_opt:     optimized phase of "Move" scheme  
    % v_Move_opt:       optimized displacement of "Move" scheme  
    % H_Move:           resulting channel of "Move" scheme 
    SNR=Ps/Pn;
    f_min=inf;
    phi_ini=exp(1j*2*pi*rand(Kr*Nr,1));
    phi_cur=phi_ini;
    v_ini=zeros(Kr,3);
    v_cur=v_ini;
    A_r_cur=A_r;
    A_r_tilde_cur=A_r_tilde;
    for iter=1:maxite
        % Phase Optimization subproblem: P2_1
        [phi_new,f_P2_1]=ManOptPhiZF(Nr,Kr,Ku,A_r_tilde_cur,SNR,Nbs,phi_cur);
        phi_cur=phi_new;
        if iter==1
            phi_PhaseOnly=phi_cur;
        end
        if f_P2_1<f_min
            f_min=f_P2_1;
            phi_Move_opt=phi_cur;v_Move_opt=v_cur;
        end
        % Displacement Optimization subproblem: P2_2
        [x_new,f_P2_2]=ManOptPosZF(v_cur,radius_move,A_r_cur,compensate_PL,Kr,Ku,Nr,Nbs,phi_cur,Ps,Pn,lambda,POS_RIS,POS_BS,POS_UE);
        v_cur=x_new;
        POS_RIS_moved = zeros(Kr,3);
        for r=1:Kr
            POS_RIS_moved(r,:) = POS_RIS(r,:) + x_new(r, : );
        end
        [~,A_r_cur,A_r_tilde_cur,~,~,~,~,~]=GenChannel(...
            POS_BS,POS_RIS_moved,POS_UE,Kr,Ku,Nr,Nrx,Nry,Nbs,lambda,compensate_PL);
        if f_P2_2<f_min
            f_min=f_P2_2;
            phi_Move_opt=phi_cur;v_Move_opt=v_cur;
        end
    end
    
    %% "Phase-only" scheme
    % generate channel with RISs' phases optimized "Phase-only" scheme
    [A_b,~,A_r_tilde,~,~,~,~,~]=GenChannel(POS_BS,POS_RIS,POS_UE,Kr,Ku,Nr,Nrx,Nry,Nbs,lambda,compensate_PL);
    % calculate the MR of "Phase-only" scheme
    PHASE_i_bar=zeros(Nr*Kr,Kr);
    for r=1:Kr
        PHASE_i_bar((r-1)*Nr+1:r*Nr,r)=phi_PhaseOnly( (r-1)*Nr+1:r*Nr );
    end
    R_matrix=PHASE_i_bar'*A_r_tilde;
    f_P2_1=real( trace(inv(R_matrix' * R_matrix)) );
    MR_PhaseOnly=log2(1+SNR*Nbs/f_P2_1);
    H_PhaseOnly = A_b * R_matrix;
    %% "Move (joint phase and dispalcement optimization)" scheme
    % generate channel with RISs' phases and displacements optimized "Move (joint phase and dispalcement optimization)" scheme
    POS_RIS_moved=POS_RIS+v_Move_opt;
    [A_b_moved,~,A_r_tilde_moved,~,~,~,~,~]=GenChannel(POS_BS,POS_RIS_moved,POS_UE,Kr,Ku,Nr,Nrx,Nry,Nbs,lambda,compensate_PL);
    % calculate the MR of "Move" scheme
    PHASE_i_bar=zeros(Nr*Kr,Kr);
    for r=1:Kr
        PHASE_i_bar((r-1)*Nr+1:r*Nr,r)=phi_Move_opt( (r-1)*Nr+1:r*Nr );
    end
    R_matrix=PHASE_i_bar'*A_r_tilde_moved;
    H_Move = A_b_moved * R_matrix;
    f_P2_2=real(trace( inv(H_Move'*H_Move) ));
    MR_Move=log2(1+SNR/f_P2_2);
end

function [v_Move,f_P2_2]=ManOptPosZF(v_ini,radius_move,A_r,compensate_PL,Kr,Ku,Nr,Nbs,vec_phi,Ps,Pn,lambda,POS_RIS,POS_BS,POS_UE)
    D_gi=cell(Kr,1);
    for r=1:Kr
        D_gi{r}=norm(POS_RIS(r,:)-POS_BS);
    end
    manifold_moving = poincareballfactory(3,Kr);
    problem_moving.M = manifold_moving;
    problem_moving.cost = @(x) mycost_moving(x,radius_move,A_r,compensate_PL,Kr,Ku,Nr,Nbs,vec_phi,Ps,Pn,lambda,POS_RIS,POS_BS,POS_UE,D_gi);
    problem_moving.egrad= @(x) myegrad_moving(x,radius_move,A_r,compensate_PL,Kr,Ku,Nr,Nbs,vec_phi,Ps,Pn,lambda,POS_RIS,POS_BS,POS_UE,D_gi);
    options.stopfun = @mystopfun;
    options.tolgradnorm=1e-6;
    options.maxiter = 500;
    checkgradient(problem_moving);
    [v_normalized, f_P2_2, ~, ~] = conjugategradient(problem_moving,v_ini.'./radius_move,options);
    v_Move=radius_move*v_normalized.';
end
function [obj]=mycost_moving(x0,radius_move,A_r,compensate_PL,Kr,Ku,Nr,Nbs,phi,Ps,Pn,lambda,POS_RIS,POS_BS,POS_UE,D_gi)
    x=radius_move*x0.';
    POS_RIS_new = zeros(Kr,3);
    for r=1:Kr
        POS_RIS_new(r,:) = POS_RIS(r,:) + x(r, : );
    end
    D_qik=zeros(Kr,Ku);
    D_0ik=zeros(Kr,Ku);
    V_ik=zeros(Kr,Ku)+1j*zeros(Kr,Ku);
    vx=zeros(Kr,Ku);
    vy=zeros(Kr,Ku);
    vz=zeros(Kr,Ku);
    for r=1:Kr
        for u=1:Ku
            D_qik(r,u)=norm(POS_RIS_new(r,:)-POS_UE(u,:));
            D_0ik(r,u)=norm(POS_RIS(r,:)-POS_UE(u,:));
            vx(r,u)=(POS_RIS(r,1) + x(r,1) - POS_UE(u,1))/D_qik(r,u);
            vy(r,u)=(POS_RIS(r,2) + x(r,2) - POS_UE(u,2))/D_qik(r,u);
            vz(r,u)=(POS_RIS(r,3) + x(r,3) - POS_UE(u,3))/D_qik(r,u);
            V_ik(r,u)=sqrt(compensate_PL)*lambda^2*exp(1j*2*pi/lambda*D_qik(r,u) )/(16*pi^2*D_gi{r}*D_0ik(r,u));
        end
    end
    PHASE_i_bar=zeros(Nr*Kr,Kr);
    for r=1:Kr
        PHASE_i_bar((r-1)*Nr+1:r*Nr,r)=phi( (r-1)*Nr+1:r*Nr );
    end
    R_tilde=V_ik.*(PHASE_i_bar'*A_r);
    tmp0=R_tilde'*R_tilde;
    obj=real( trace( inv(tmp0) ) );
end
function [egrad]=myegrad_moving(x0,radius_move,A_r,compensate_PL,Kr,Ku,Nr,Nbs,phi,Ps,Pn,lambda,POS_RIS,POS_BS,POS_UE,D_gi)
    x=radius_move*x0.';
    POS_RIS_new = zeros(Kr,3);
    for r=1:Kr
        POS_RIS_new(r,:) = POS_RIS(r,:) + x(r, : );
    end
    D_qik=zeros(Kr,Ku);
    D_0ik=zeros(Kr,Ku);
    V_ik=zeros(Kr,Ku)+1j*zeros(Kr,Ku);
    vx=zeros(Kr,Ku);
    vy=zeros(Kr,Ku);
    vz=zeros(Kr,Ku);
    for r=1:Kr
        for u=1:Ku
            D_qik(r,u)=norm(POS_RIS_new(r,:)-POS_UE(u,:));
            D_0ik(r,u)=norm(POS_RIS(r,:)-POS_UE(u,:));
            vx(r,u)=(POS_RIS(r,1) + x(r,1) - POS_UE(u,1))/D_qik(r,u);
            vy(r,u)=(POS_RIS(r,2) + x(r,2) - POS_UE(u,2))/D_qik(r,u);
            vz(r,u)=(POS_RIS(r,3) + x(r,3) - POS_UE(u,3))/D_qik(r,u);
            V_ik(r,u)=sqrt(compensate_PL)*lambda^2*exp(1j*2*pi/lambda*D_qik(r,u) )/(16*pi^2*D_gi{r}*D_0ik(r,u));
        end
    end
    PHASE_i_bar=zeros(Nr*Kr,Kr);
    for r=1:Kr
        PHASE_i_bar((r-1)*Nr+1:r*Nr,r)=phi( (r-1)*Nr+1:r*Nr );
    end
    R_tilde=V_ik.*(PHASE_i_bar'*A_r);
    tmp0=R_tilde'*R_tilde;
    B=-(inv(tmp0)*inv(tmp0)*R_tilde').';
    egrad=zeros(Kr,3);
    for r=1:Kr
        bx = 1j*2*pi/lambda * (R_tilde(r,:).' .* vx(r,:).');
        by = 1j*2*pi/lambda * (R_tilde(r,:).' .* vy(r,:).');
        bz = 1j*2*pi/lambda * (R_tilde(r,:).' .* vz(r,:).');
        egrad(r,:) =  2 * real(B(r,:) * [bx,by,bz]); 
    end
    egrad = radius_move*egrad.';
end
function [vec_phi,f_P2_1]=ManOptPhiZF(Nr,Kr,Ku,A_r_tilde,SNR,Nbs,phase_ini)
    manifold = complexcirclefactory(Kr * Nr);
    problem.M = manifold;
    problem.cost = @(phi) mycost_ZF(A_r_tilde, Nr, Kr, phi);
    problem.egrad = @(phi) myegrad_ZF(A_r_tilde, Nr, Kr, phi);
    options.stopfun = @mystopfun;
    options.tolgradnorm=1e-6;
    options.maxiter = 1000;
    checkgradient(problem);
    [vec_phi, f_P2_1, info2, options2] = conjugategradient(problem,[phase_ini],options);
end
function [obj] = mycost_ZF(A_r_tilde, Nr, Kr, phi)
    PHASE_i_bar=zeros(Nr*Kr,Kr);
    for r=1:Kr
        PHASE_i_bar((r-1)*Nr+1:r*Nr,r)=phi( (r-1)*Nr+1:r*Nr );
    end
    tmp0=A_r_tilde'*PHASE_i_bar*PHASE_i_bar'*A_r_tilde;
    obj=real(trace( inv(tmp0) ));
end
function [egrad] = myegrad_ZF(A_r_tilde, Nr, Kr, phi)
    PHASE_i_bar=zeros(Nr*Kr,Kr);
    for r=1:Kr
        PHASE_i_bar((r-1)*Nr+1:r*Nr,r)=phi( (r-1)*Nr+1:r*Nr );
    end
    tmp0=A_r_tilde'*PHASE_i_bar*PHASE_i_bar'*A_r_tilde;
    egrad=zeros(Kr*Nr,1)+1j*zeros(Kr*Nr,1);
    tmp=-A_r_tilde/(tmp0)/(tmp0)*A_r_tilde'*PHASE_i_bar;
    for r=1:Kr
        egrad( (r-1)*Nr+1:r*Nr ) = 2*tmp( (r-1)*Nr+1:r*Nr, r );
    end
end