function a_N=SteerVectorPlanar(Nx,Nz,POS_tx,POS_rx,lambda,In_Out)
% RIS UPA placed in x-z plane, orient toward the negative y-axis direction
N=Nx*Nz;
d=lambda/2;
vec=POS_rx-POS_tx;
cos_phi=vec(1)/sqrt( vec(1)^2+vec(2)^2 );
cos_theta=vec(3)/sqrt( vec(1)^2+vec(2)^2+vec(3)^2 );
sin_theta=sqrt( vec(1)^2+vec(2)^2 )/sqrt( vec(1)^2+vec(2)^2+vec(3)^2 );
vec_x=exp(-1j*In_Out*2*pi/lambda*[0:Nx-1]*d*cos_phi*sin_theta);
vec_z=exp(-1j*In_Out*2*pi/lambda*[0:Nz-1]*d*cos_theta);
a_N=kron(vec_z,vec_x).';
end