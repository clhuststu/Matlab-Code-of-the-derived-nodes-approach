%**************************************************************************
% This matlab code will illustrate the realization of the previous second-
% order temporal finite element method. After running the code, you wll see
% the stability lobes prediction of the slot milling process (a/D=1), as 
% shown in the lower left of table 2 in our manuscript.
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
m=40;     % number of discretization interval of the forced vibration stage
% Dsicretization of the specific cutting force coefficient h(t)
for i=1:m+1
    dtr=(fiex-fist)/m;     % ¦¤t
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
%------------Begin of the previous 2nd-TFEA-----------------
A=[-zeta*w0,1/m_t;...
    m_t*((zeta*w0)^2-w0^2),-zeta*w0];      % A in Eq.(26)
I=eye(size(A));
for x=2:stx+1                              % sweeping spindle speeds
    o=o_st+(x-1)*(o_fi-o_st)/stx;          % the spindle speed ¦¸
    T=60/o/N;                              % time delay
    tf=T*(1-(fiex-fist)*N/(2*pi));% time length of the free vibration stage
    F0=expm(A*tf);
    tau=(T-tf)/m;                          % time step
    %----the constant coefficients for calclating the Matrix P and N----
    C1_12=1;C1_11=0;C1_10=-1;
    C1_22=-2/3;C1_21=4/3;C1_20=-2/3;
    C2_12=-tau/3;C2_11=-4*tau/3;C2_10=-tau/3;
    C2_22=tau/3;C2_21=0;C2_20=-tau/3;
    C3_12=tau/24;C3_11=-tau/4;C3_10=-7*tau/24;
    C3_22=tau/40;C3_21=-2*tau/15;C3_20=-9*tau/40;
    C4_12=-tau/12;C4_11=-5*tau/6;C4_10=-tau/12;
    C4_22=tau/12;C4_21=0;C4_20=-tau/12;
    C5_12=-7*tau/24;C5_11=-tau/4;C5_10=tau/24;
    C5_22=9*tau/40;C5_21=2*tau/15;C5_20=-tau/40;
    %--------------------------------------------------------------------
    for y=1:sty+1                           % sweeping depth of cut
        w=w_st+(y-1)*(w_fi-w_st)/sty;       % the depth of cut
        %------------Construction of the mapping matrix---------------
        Gamma=zeros(2*(m+1),2*(m+1));
        Sigma=zeros(2*(m+1),2*(m+1));
        Gamma(2*(m+1)-1:2*(m+1),2*(m+1)-1:2*(m+1))=I;
        Sigma(2*(m+1)-1:2*(m+1),1:2)=F0;
        for i=1:m/2
            B0i=[0 0;-w*h(m-2*i+1) 0];      % Bi
            B1i=[0 0;-w*h(m-2*i+2) 0];      % Bi+1
            B2i=[0 0;-w*h(m-2*i+3) 0];      % Bi+2
            %----calculation of the coefficient matrix P and N-------
            P12=C3_12*B0i+C4_12*B1i+C5_12*B2i;
            P11=C3_11*B0i+C4_11*B1i+C5_11*B2i;
            P10=C3_10*B0i+C4_10*B1i+C5_10*B2i;
            P22=C3_22*B0i+C4_22*B1i+C5_22*B2i;
            P21=C3_21*B0i+C4_21*B1i+C5_21*B2i;
            P20=C3_20*B0i+C4_20*B1i+C5_20*B2i;
            N12=C1_12*I+C2_12*A+P12;
            N11=C1_11*I+C2_11*A+P11;
            N10=C1_10*I+C2_10*A+P10;
            N22=C1_22*I+C2_22*A+P22;
            N21=C1_21*I+C2_21*A+P21;
            N20=C1_20*I+C2_20*A+P20;
            %--The end of calculation of the coefficient matrix P and n----
            Gamma(2*(2*i-2)+1:2*(2*i-1),2*(2*i-2)+1:2*(2*i-1))=N12;
            Gamma(2*(2*i-2)+1:2*(2*i-1),2*(2*i-1)+1:2*(2*i))=N11;
            Gamma(2*(2*i-2)+1:2*(2*i-1),2*(2*i)+1:2*(2*i+1))=N10;
            Gamma(2*(2*i-1)+1:2*(2*i),2*(2*i-2)+1:2*(2*i-1))=N22;
            Gamma(2*(2*i-1)+1:2*(2*i),2*(2*i-1)+1:2*(2*i))=N21;
            Gamma(2*(2*i-1)+1:2*(2*i),2*(2*i)+1:2*(2*i+1))=N20;
            Sigma(2*(2*i-2)+1:2*(2*i-1),2*(2*i-2)+1:2*(2*i-1))=P12;
            Sigma(2*(2*i-2)+1:2*(2*i-1),2*(2*i-1)+1:2*(2*i))=P11;
            Sigma(2*(2*i-2)+1:2*(2*i-1),2*(2*i)+1:2*(2*i+1))=P10;
            Sigma(2*(2*i-1)+1:2*(2*i),2*(2*i-2)+1:2*(2*i-1))=P22;
            Sigma(2*(2*i-1)+1:2*(2*i),2*(2*i-1)+1:2*(2*i))=P21;
            Sigma(2*(2*i-1)+1:2*(2*i),2*(2*i)+1:2*(2*i+1))=P20;
        end
        %-------The end of construction of the mapping matrix----------
        Phi=inv(Gamma)*Sigma;
        ss(x,y)=o;                         % matrix of spindle speeds
        dc(x,y)=w;                         % matrix of depth of cut
        ei(x,y)=max(abs(eig(Phi)));        % matrix of eigenvalues
    end
    stx+1-x
end
%-----------The end of the previous 2nd-TFEA-------------
toc;
contour(ss,dc,ei,[1,1],'r-','linewidth',1.5);
xlabel('\Omega(rpm)');ylabel('w(m)');
xlim([5000,7500]);
title('m=40');