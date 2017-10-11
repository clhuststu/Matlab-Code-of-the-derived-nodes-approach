%**************************************************************************
% This matlab code will illustrate the application of the derived nodes
% approach to the temporal finite element method. After running the code,
% you wll see the stability lobes prediction of the slot milling process 
% (a/D=1), as shown in the lower right of table 2 in our manuscript.
%**************************************************************************
close all;
clear all;
clc;
tic
% the modal parameters are as follows, which can also be found in Ref.[14]
N=2;      % number of teeth
Kt=6e8;   % tangential cutting force coefficient (N/m^2)
Kn=2e8;   % normal cutting force coefficient (N/m^2)
w0=922*2*pi;   % angular natural frequency (rad/s)
zeta=0.011;    % relative damping (1)
m_t=0.03993;   % mass (kg)
aD=1;          % radial immersion ratio a/D
up_or_down=-1;      % 1: up-milling,-1: down-milling
if up_or_down==1    % up_milling
    fist=0;         % start angle
    fiex=acos(1-2*aD);     % exit angle
elseif up_or_down==-1      % down-milling
    fist=acos(2*aD-1);     % start angle
    fiex=pi;               % exit angle
end
stx=200;                   % steps of spindle speed
sty=100;                   % steps of depth of cut
w_st=0e-3;                 % starting depth of cut
w_fi=4e-3;                 % final depth of cut
o_st=5e3;                  % starting spindle speed
o_fi=7.5e3;                % final spindle speed

% the computional parameters
m=40;
% the recursive matrix
P=eye(m);           % matrix P in Eq.(7)
p=ones(m-1,1);
P=P+diag(p,1);
Q=zeros(m,m+2);     % matrix Q in Eq.(7)
q=ones(m-1,1);
Q(1:m-1,1:m+1)=full(spdiags([1/4*q ...
    3/2*q 1/4*q],0:2,m-1,m+1));
Q(m,m+2)=1;
H=inv(P)*Q;         % matrix H in Eq.(8)
HkronI=kron(H,I);
% Dsicretization of the specific cutting force coefficient h(t)
for i=1:m+1
    dtr=(fiex-fist)/m;     % the discretization step
    h(i)=0;              
    for j=1:N              % loop for tooth j
        fi=fist+(i-1)*dtr+(j-1)*2*pi/N;
        if (fi>=fist)*(fi<=fiex)
            g=1;           % tooth is in the cut
        else
            g=0;           % tooth is out of cut
        end
        h(i)=h(i)+g*(Kt*cos(fi)+Kn*sin(fi))*sin(fi);
    end
end
%------------Begin of the proposed 2nd-TFEA-----------------
A=[-zeta*w0,1/m_t;...
    m_t*((zeta*w0)^2-w0^2),-zeta*w0];      % A in Eq.(26)
I=eye(size(A));
for x=2:stx+1                              % sweeping spindle speeds
    o=o_st+(x-1)*(o_fi-o_st)/stx;          % the spindle speed ¦¸
    T=60/o/N;                              % time delay
    tf=T*(1-(fiex-fist)*N/(2*pi));
    F0=expm(A*tf);
    tau=(T-tf)/m;                          % time step
    %----the constant coefficients for calclating the Matrix P and N----
    C12=1;C11=0;C10=-1;
    C22=-tau/6;C21=-2*tau/3;C20=-tau/6;
    C32=0;C31=-tau/3;C30=-tau/6;
    C42=-tau/6;C41=-tau/3;C40=0;
    %--------------------------------------------------------------------
    for y=1:sty+1                          % sweeping depth of cut
        w=w_st+(y-1)*(w_fi-w_st)/sty;      % the depth of cut
        %------------Construction of the mapping matrix---------------
        Gamma=zeros(2*(m+2),2*(m+2));
        Sigma=zeros(2*(m+2),2*(m+2));
        Gamma(2*(m+1)-1:2*(m+1),2*(m+1)-1:2*(m+1))=I;
        Sigma(2*(m+1)-1:2*(m+1),1:2)=F0;
        for i=1:m
            B0i=[0 0;-w*h(m-i+1) 0];       % Bi in Eq.(28)
            B1i=[0 0;-w*h(m-i+2) 0];       % Bi+1 in Eq.(28)
            P2i=C32*B0i+C42*B1i;           % P2m in Eq.(21)
            P1i=C31*B0i+C41*B1i;           % P1m in Eq.(21)
            P0i=C30*B0i+C40*B1i;           % P0m in Eq.(21)
            N2i=C12*I+C22*A+P2i;           % N2m in Eq.(21)
            N1i=C11*I+C21*A+P1i;           % N1m in Eq.(21)
            N0i=C10*I+C20*A+P0i;           % N0m in Eq.(21)
            Gamma(2*i-1:2*i,2*i-1:2*i)=N2i+N1i*H(i,i);         
            Gamma(2*i-1:2*i,2*i+1:2*(i+1))=N0i+N1i*H(i,i+1);   
            Gamma(2*i-1:2*i,2*i+3:2*(m+2))=N1i*HkronI(2*i-1:2*i,2*i+3:2*(m+2));
            Sigma(2*i-1:2*i,2*i-1:2*i)=P2i+P1i*H(i,i);         
            Sigma(2*i-1:2*i,2*i+1:2*i+2)=P0i+P1i*H(i,i+1);     
            Sigma(2*i-1:2*i,2*i+3:2*(m+2))=P1i*HkronI(2*i-1:2*i,2*i+3:2*(m+2));
        end
        %-----Begin of the extra equations as shown in Eq.(25)-------- 
        B0=[0 0;-w*h(1) 0];
        B1=[0 0;-w*h(2) 0];
        P2d=tau/30*B0+2*tau/15*B1;
        P1d=-tau/15*B0+tau/15*B1;
        P0d=-2*tau/15*B0-tau/30*B1;
        N2d=-2/3*I+tau/6*A+P2d;
        N1d=4/3*I+P2d;
        N0d=-2/3*I-tau/6*A+P0d;
        Gamma(2*(m+2)-1:2*(m+2),2*m-1:2*m)=N2d;
        Gamma(2*(m+2)-1:2*(m+2),2*(m+1)-1:2*(m+1))=N0d;
        Gamma(2*(m+2)-1:2*(m+2),2*(m+2)-1:2*(m+2))=N1d;
        Sigma(2*(m+2)-1:2*(m+2),2*m-1:2*m)=P2d;
        Sigma(2*(m+2)-1:2*(m+2),2*(m+1)-1:2*(m+1))=P0d;
        Sigma(2*(m+2)-1:2*(m+2),2*(m+2)-1:2*(m+2))=P1d;
        %-----The end of the extra equations as shown in Eq.(25)-------
        Phi=inv(Gamma)*Sigma;
        ss(x,y)=o;                          % matrix of spindle speeds
        dc(x,y)=w;                          % matrix of depth of cut
        ei(x,y)=max(abs(eig(Phi)));         % matrix of eigenvalues
    end
    stx+1-x
end
%-----------The end of the proposed 2nd-TFEA-------------
toc
contour(ss,dc,ei,[1,1],'r-','linewidth',1.5);
xlabel('\Omega(rpm)');ylabel('w(m)');
xlim([5000,7500]);
title('m=40');
