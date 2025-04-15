function [POS_RIS]=getRISPOS(POS_BS,radius,Nbx,Nbz)
% Requirements: Nbx = Nbz > 2
    D_BS2RIS=5;
    Nidx=floor(sqrt(2)*Nbz/4);
    cos_theta=zeros(Nidx,1);
    sin_theta=zeros(Nidx,1);
    cos_phi=zeros(Nidx,1);
    sin_phi=zeros(Nidx,1);
    theta_1=pi/2;
    phi_1=pi/2;
    for idx_theta=1:Nidx
        cos_theta(idx_theta)=idx_theta*2/Nbz;
        sin_theta(idx_theta)=sqrt(1-cos_theta(idx_theta)^2);
        cos_phi(idx_theta)=idx_theta*2/Nbz/sin_theta(idx_theta);
        sin_phi(idx_theta)=sqrt(1-cos_phi(idx_theta)^2);
    end
    cot_phi=cos_phi./sin_phi;
    cot_theta=cos_theta./sin_theta;
    x=zeros(4*Nidx,1);
    y=zeros(4*Nidx,1);
    z=zeros(4*Nidx,1);
    for idx=1:Nidx
        vec=[cot_phi(idx)*D_BS2RIS, D_BS2RIS, norm([cot_phi(idx)*D_BS2RIS,D_BS2RIS])*cot_theta(idx)];
        vec=radius*( vec./norm(vec) );
        x0=vec(1);y0=vec(2);z0=vec(3);
        x(4*(idx-1)+1:4*idx)=[x0;-x0;x0;-x0];
        y(4*(idx-1)+1:4*idx)=[y0;y0;y0;y0];
        z(4*(idx-1)+1:4*idx)=[z0;-z0;-z0;z0];
    end
    POS_RIS=POS_BS+[x,y,z];
end