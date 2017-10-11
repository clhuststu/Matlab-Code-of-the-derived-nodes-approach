%**************************************************************************
% This matlab code will illustrate the application of the derived nodes 
% approach to the full discretization method. After running the code, 
% you wll see the stability lobes prediction of the slot milling  process 
% (a/D=1), as shown in the lower right of table 1 in our manuscript.
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

% the computational parameters
m=40;               % number of discretization interval
D=zeros(m+3,m+3);   % matrix D  
d=ones(m+2,1);
d(1:2)=0;
D=D+diag(d,-1);
D(3,1)=1;
D(m+3,m:m+3)=[1/4 3/2 1/4 -1];
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
% Discretization of the specific cutting force coefficient h(t) in Eq.(25)
for i=1:m+1
    dtr=2*pi/N/m;   % ¦¤t
    h(i)=0;
    for j=1:N       % loop for tooth j
        fi=i*dtr+(j-1)*2*pi/N;
        if (fi>=fist)*(fi<=fiex)
            g=1;    % tooth is in the cut
        else
            g=0;    % tooth is out of cut
        end
        h(i)=h(i)+g*(Kt*cos(fi)+Kn*sin(fi))*sin(fi);
    end
end
%------------Begin of the proposed 2nd-FDM------------
A=[-zeta*w0,1/m_t;...
    m_t*((zeta*w0)^2-w0^2),-zeta*w0];   % A in Eq.(26)
I=eye(size(A));
invA=inv(A);
% start of computation
for x=2:stx+1       % sweeping spindle speeds
    o=o_st+(x-1)*(o_fi-o_st)/stx;  % the spindle speed
    T=60/o/N;                      % time delay
    tau=T/m;                       % time step
    %---------Calculation of F0,F1,F2,F3,F4----------
    F0=expm(A*tau);                % F0 in Eq.(11)
    F1=invA*(F0-I);
    F2=invA*(F0*tau-F1);
    F3=invA*(F0*tau^2-2*F2);
    F4=invA*(F0*tau^3-3*F3);
    f0=expm(A*tau/2);
    %----The end of calculation of F0,F1,F2,F3,F4----
    G10=(2*F4-tau*F3)/tau^3;       % G10 in Eq.(11)
    G11=(2*F3-tau*F2)/tau^2-G10;   % G11 in Eq.(11)
    G20=4*(tau*F3-F4)/tau^3;       % G20 in Eq.(11)
    G21=4*(tau*F2-F3)/tau^2-G20;   % G21 in Eq.(11)
    G30=(2*F4-3*tau*F3+tau^2*F2)/tau^3;  % G30 in Eq.(11)
    G31=(2*F3-3*tau*F2+tau^2*F1)/tau^2 ...
    -G30;                                % G40 in Eq.(11)
    for y=1:sty+1                 % sweeping depth of cut
        w=w_st+(y-1)*(w_fi-w_st)/sty;  % the depth of cut
        Fi=eye(m+3,m+3); % construct transition matrix Fi
        for i=1:m
            B0i=[0 0;-w*h(i) 0];       % Bi in Eq.(28)
            B1i=[0 0;-w*h(i+1) 0];     % Bi+1 in Eq.(28)
            %----the coefficients matrix for calculating N1i,N2i,N3i-----
            Cons1=G11*B1i+G10*B0i;
            Cons2=G21*B1i+G20*B0i;
            Cons3=G31*B1i+G30*B0i;
            Li=inv(I-Cons3-Cons2/4);   % Li in Appendix A
            %------------------------------------------------------------
            N2i=Li*Cons2;              % N2i in Eq.(12)
            %--Block multiplication of Di+Ri x kron(H,I)--
            D(1:2,1:2)=Li*(F0+Cons1)+(3/2-H(1,1))*N2i;
            D(1:2,3)=N2i(1:2,1:1)*(1/4-H(1,2));
            D(1:2,4:m)=-N2i(1:2,1:1)*H(1,3:m-1);
            D(1:2,m+1)=-Li*Cons3(1:2,1:1)-N2i(1:2,1:1)...
            *H(1,m);
            D(1:2,m+2)=-Li*Cons1(1:2,1:1)-N2i(1:2,1:1)...
            *H(1,m+1);
            D(1:2,m+3)=-N2i(1:2,1:1)-N2i(1:2,1:1)...
            *H(1,m+2);
            %---end of calculation of Di+Ri x kron(H,I)--
            Fi=D*Fi;
        end
        %-------------Construction of Phi---------------
        Phi=zeros(m+4,m+4);
        Phi(1:m+2,1:m+3)=Fi(1:m+2,1:m+3);
        B0=[0 0;-w*h(1) 0];
        Bd=[0 0;-w*(h(1)+h(2))/2 0];
        invImBd=inv(I-Bd*tau/4);
        Phi(m+3:m+4,1:2)=invImBd*f0*(I+B0*tau/4);
        Phi(m+3:m+4,m+2)=-invImBd*f0*B0(1:2,1:1)*tau/4;
        Phi(m+3:m+4,m+3:m+4)=-invImBd*Bd*tau/4;
        %--------The end of Construction of Phi---------
        ss(x,y)=o;            % matrix of spindle speeds
        dc(x,y)=w;            % matrix of depth of cut
        ei(x,y)=max(abs(eig(Phi)));       % eigenvalues
    end
    stx+1-x
end
%-----------The end of the proposed 2nd-FDM-------------
toc
contour(ss,dc,ei,[1,1],'r-','linewidth',1.5);
xlabel('\Omega(rpm)');ylabel('w(m)');
xlim([5000,7500]);
title('m=40');