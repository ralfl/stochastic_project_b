clear all
close all

Ftsize=18;

om_g=18;  %[]
zeta =0.02;
zeta_g=0.3;
xi=100;    % [mm]
S0=1;        % [m^2/s]
m=400000;   % [kg]


Fs = 100;                   % Sampling frequency
T = 1/Fs;                   % Sample time
L = 2500;                   % Length of signal
t = 0:T:L*T-T;     % Time vector
dt=T;
f_cut=10;
om=0:(f_cut/L*2*pi):f_cut*2*pi-(f_cut/L*2*pi);

e=4*(exp(-.25*t)-exp(-.5*t));

% for ii=1:length(L)
% e(ii)=4*(exp(-.25*t(ii)-exp(-.5*t(ii))));
% end

figure 
subplot(2,2,2)
plot(t,e)
title('Modulating function e(t)','Fontsize',Ftsize+2,'FontWeight','n','FontName','Times New Roman')
xlabel('Time [s]','Fontsize',Ftsize,'FontWeight','n','FontName','Times New Roman')
ylabel('Envelope [-]','Fontsize',Ftsize,'FontWeight','n','FontName','Times New Roman')


for i=1:length(om)
Sbb(i)=S0*(4*zeta_g^2*om_g^2*om(i)^2+om_g^4)/((om_g^2-om(i)^2)^2+4*zeta_g^2*om_g^2*om(i)^2);
end

subplot(2,2,1)
plot(om,Sbb);
title('PSD  Kanai-Tajimi Filter','Fontsize',Ftsize+2,'FontWeight','n','FontName','Times New Roman')
xlabel('Frequency \omega [rad/s]','Fontsize',Ftsize,'FontWeight','n','FontName','Times New Roman')
ylabel('[m^2/s^4/Hz]','Fontsize',Ftsize,'FontWeight','n','FontName','Times New Roman')

SQPS=(Sbb*pi/T).^.5;
SQPS_re=randn(1,L).*SQPS;
SQPS_im=randn(1,L).*SQPS;
B=complex(SQPS_re,SQPS_im);




b=real(ifft(B,[],L));
subplot(2,2,3)
plot(t,b)
title('Inverse FFT from Spectrum  b(t)','Fontsize',Ftsize+2,'FontWeight','n','FontName','Times New Roman')
xlabel('Time [s]','Fontsize',Ftsize,'FontWeight','n','FontName','Times New Roman')
ylabel('Acceleration [m/s^2]','Fontsize',Ftsize,'FontWeight','n','FontName','Times New Roman')

a_F=e.*b;
F=a_F.*m;

subplot(2,2,4) 
plot(t,F)
title('Earthquake Force  F=m\cdote(t)\cdot b(t)','Fontsize',Ftsize+2,'FontWeight','n','FontName','Times New Roman')
xlabel('Time [s]','Fontsize',Ftsize,'FontWeight','n','FontName','Times New Roman')
ylabel('Force [N]','Fontsize',Ftsize,'FontWeight','n','FontName','Times New Roman')

%% Central diffrence method SDOF

E=210e9;   % [N/m^2)
b=.7;        % [m]
H=6;          % [m]
I=b^4/12;   % [m^4]    
k=4*12*E*I/H^3;
c=2*zeta*sqrt(k*m);
a=zeros(L,1);
x=zeros(L,1);
x(1)=0;
a(1)=0;%-(c*v0+k*x(1))/m;

for i=2:L-1
    a(i)=F(i)-(k-2*m/dt^2)*x(i)-(m/dt^2-c/2/dt)*x(i-1);
    x(i+1)=a(i)'/(m/(dt^2)+c/2/dt);
end

figure
plot(t,x)
title('Time response SDOF','Fontsize',Ftsize+2,'FontWeight','n','FontName','Times New Roman')
xlabel('Time [s]','Fontsize',Ftsize,'FontWeight','n','FontName','Times New Roman')
ylabel('Displacement [m]','Fontsize',Ftsize,'FontWeight','n','FontName','Times New Roman')

%%