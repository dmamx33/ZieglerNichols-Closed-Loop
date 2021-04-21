%%
clear all
close all 
clc
month = 12; %%%
pole = 2;
%%%%%%%%%%%%%%%% Symbolic Calculations %%%%%%%%%%%%%%%%
syms Kp w real
syms s G 
G1(s)=1/s;
G2(s)=1/(s+month);
G3(s)=1/(s+pole);
G(s)=G1(s) * G2(s) * G3(s);
G(s)= expand(G(s));
%G(s)=1/((6*s^3)+(11*s^2)+(6*s)+(6))%%%%%%QUITAR
[numg,deng] = numden(G(s));
num = sym2poly(numg);
den = sym2poly(deng);
disp('Plant:')
pretty(G(s))
%%G(s)=1/((6*s^3)+(11*s^)+(6*s)+(6))
%F(s)=expand(Kp*G(s)/((Kp*G(s))+1));
F(s)=simplify(Kp*G(s)/((Kp*G(s))+1));
disp('Feedback plant with Proportional Kp:')
pretty(F(s));
[numF,denF] = numden(F(s));
CE = subs(denF,s,i*w);
CE_r = real (CE);
CE_i = imag (CE);
Critical =  solve([CE_r==0,CE_i==0],[Kp w] );
ArrayKsim = Critical.Kp;
Arraywsim = Critical.w;
ArrayK = double(ArrayKsim);
ArrayW = double(Arraywsim);
Kcr=ArrayK(find(ArrayK > 0,1,'first'));
wcr=ArrayW(find(ArrayW > 0,1, 'first'));
%Kcr = double(Critical.Kp(1,1));
%wcr = double(Critical.w(1,1));

% if Kcr < 0
%     Kcr = abs(Kcr);
% end
% if wcr < 0
%     wcr = abs(wcr);
% end
 Tcr = 2*pi/wcr;
% return
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 disp('Critial Values ZN closed loop test:')
% Kcr = (month+pole)*2*month
% Tcr=2*pi/(sqrt(month*2))
Kcr
Tcr
%%%%%%%%%%%%%%%%% ZIEGLER NICHOLS %%%%%%%%%%%%%%%%%%
disp('PID Constants Values ZN closed loop test:')
PZNprop = 0.5 * Kcr;        %% P 
PZNpi = 0.45*Kcr;           %% PI
P= 0.6 * Kcr            %% PID
tau_i = 0.5 *Tcr
tau_d = 0.125*Tcr
%%%%%%% ZIEGLER NICHOLS MODIFIED SOME OVERSHOOT - PID
% disp('PID Constants Values ZN Modified PID Some overshoot:')
% P_MZN_SO= 0.33 * Kcr
% TI_MZN_SO = 0.5 *Tcr
% TD_MZN_SO = Tcr / 3
% PIDMZN = tf([P_MZN_SO*tau_d2*tau_i2 P_MZN_SO*tau_i2 P_MZN_SO],[0 tau_i2 0]);
% compensatedPLANTMZN = feedback (series(PID2,PLANT),1);
% YcompMZN = lsim(compensatedPLANT2,Ustep,t);
% %%%%%%% ZIEGLER NICHOLS MODIFIED NO OVERSHOOT
% disp('PID Constants Values ZN Modified PID no overshoot:')
% P_MZN_NO= 0.2 * Kcr
% TI_MZN_NO = 0.5 *Tcr
% TD_MZN_NO = Tcr / 3
%%%%%%% TYREUS LUYBEN PI
% disp('PID Constants Values TL for PI:')
% P_TL_PI=  Kcr / 3.2
% TI_TL_PI = 2.2 *Tcr
% %%%%%%% TYREUS LUYBEN PID
% disp('PID Constants Values TL for PID:')
% P_TL_PID=  Kcr / 3.2
% TI_TL_PID = 2.2 *Tcr
% TD_TL_PID = Tcr/6.3
%%%%%%%%%%%%%%%%%% TIME VECTOR %%%%%%%%%%%%%%%%%%%%%
Ts=0.001;
It=0;
Ft=100;
t=It:Ts:Ft-Ts;
%%%%%%%%%%%%%%%%%% INPUT U VECTOR %%%%%%%%%%%%%%%%%%%%%
Ustep = ones(length(t),1);
%%%%%%%%%%%%%%%%%% PLANTS  %%%%%%%%%%%%%%%%%%%%%
PLANT = tf(num,den);
fPLANT=feedback(PLANT,1);
fcritPLANT=feedback( series(Kcr,PLANT) , 1 );
PID = tf([P*tau_d*tau_i P*tau_i P],[0 tau_i 0]);
compensatedPLANT = feedback (series(PID,PLANT),1);
%%%%%%%%%%%%%%%%%% SIMULATIONS  %%%%%%%%%%%%%%%%%%%%%
Y = lsim(PLANT,Ustep,t);
Yf = lsim(fPLANT,Ustep,t);
Ycrit = lsim(fcritPLANT,Ustep,t);
Ycomp = lsim(compensatedPLANT,Ustep,t);
%%%%%%%%%%%%%%%%%%%%   Manual readjust   %%%%%%%%%%%%%%%%%%%%%% 
disp('PID Constants Values ZN closed loop test with manual modifications:')
%%%%%lla no mhober
P_2 = P * 0.5 %%%%%% MODIFY TO ACHIEVE YOUR DESIRED BEHAVIOR 0.4 %%%% 0.5 
tau_i2 = tau_i * 20%%%%%% MODIFY TO ACHIEVE YOUR DESIRED BEHAVIOR 10  %%% 30
tau_d2 = tau_d * 3%%%%%% MODIFY TO ACHIEVE YOUR DESIRED BEHAVIOR 3  %%%% 3
PID2 = tf([P_2*tau_d2*tau_i2 P_2*tau_i2 P_2],[0 tau_i2 0]);
compensatedPLANT2 = feedback (series(PID2,PLANT),1);
Ycomp2 = lsim(compensatedPLANT2,Ustep,t);
%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(t,Y)
hold on
plot(t,Yf)
plot(t,Ycrit)
plot(t,Ycomp)
plot(t,Ycomp2,'r')
plot(t,Ustep,'--')
axis([0 25 0 2])
grid 
legend('Output Uncompensated','Output feedbacked','Output with crictical gain','Output with PID ZN','Output with PID ZN manual modification','Input')
xlabel('Time')
ylabel('Output')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%plot(t,Y)
hold on
%plot(t,Yf)
plot(t,Ycrit)
%plot(t,Ycomp)
%plot(t,Ycomp2,'r')
plot(t,Ustep)
axis([0 10 0 2])
grid 
title('Testing plant with critic gain')
legend('Output','Input')
xlabel('Time')
ylabel('Amplitude')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%plot(t,Y)
hold on
%plot(t,Yf)
%plot(t,Ycrit)
plot(t,Ycomp)
%plot(t,Ycomp2,'r')
plot(t,Ustep)
axis([0 10 0 2])
grid 
title('Plant with PID and ZN tunning')
legend('Output','Input')
xlabel('Time')
ylabel('Amplitude')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(t,Y,'LineWidth',2)
hold on
plot(t,Yf,'LineWidth',2)
%plot(t,Ycrit)
plot(t,Ycomp,'LineWidth',2)
plot(t,Ycomp2,'LineWidth',2)
plot(t,Ustep,'--')
axis([0 10 0 2])
grid 
title('Plant with PID and ZN tunning and manual tunning')
legend('Plant without PID open loop','Plant without PID closed loop',' Plant closed loop with ZN PID',' Plant closed loop with ZN PID and manual modifications','Input')
xlabel('Time')
ylabel('Amplitude')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
