function combined_smc_plots

tspan = [0 10];
r0 = [0.3; 0.5; 0.2];

[t1,r_conv] = ode45(@convsmcfun,tspan,r0);
K=6;
eta = 4;
c1=10;
c2=5;
Zeta_Conv = c1*r_conv(:,1) + c2*r_conv(:,2) + r_conv(:,3);
U_conv = -K*sign(Zeta_Conv) + eta*Zeta_Conv;

[t2,r_cont] = ode45(@contsmcfun,tspan,r0);
k= 1.0/56.7;
Zeta_cont = c1*r_cont(:,1) + c2*r_cont(:,2) + r_cont(:,3);
U_cont = -15*sign(r_cont(:,1)).*(abs(r_cont(:,1))).^(1/4)-23*sign(r_cont(:,2)).*(abs(r_cont(:,2))).^(1/3)-9*sign(r_cont(:,3)).*(abs(r_cont(:,3))).^(1/2);
%f = -k*sin(acos(-r_cont(:,3)/k)).*r_cont(:,3).*tan(asin(-r_cont(:,2)));
%g = k*sin(acos(-r_cont(:,3)/k));

r0 = [0.3; 0.5; 0.2; 0];
[t3,r_ismc] = ode45(@ismcfun,tspan,r0);
Zeta_ismc = r_ismc(:,4);
a1 = 15;
a2 = 23;
a3 = 9;
b1 = 1/4;
b2 = 1/3;
b3 = 1/2;
a4 = 10;
U_nom_ismc = -a1*sign(r_ismc(:,1)).*(abs(r_ismc(:,1))).^(b1)-a2*sign(r_ismc(:,2)).*(abs(r_ismc(:,2))).^(b2)-a3*sign(r_ismc(:,3)).*(abs(r_ismc(:,3))).^(b3);
U_ismc = U_nom_ismc -a4*sign(Zeta_ismc);


r0 = [0.3; 0.5; 0.2; 0; 0];
[t4,r_sta] = ode45(@stafun,tspan,r0);
Zeta_sta = r_sta(:,4);
Theta_sta = r_sta(:, 5);
a1 = 15;
a2 = 23;
a3 = 9;
b1 = 1/4;
b2 = 1/3;
b3 = 1/2;
a5 = 0.13;

U_nom_sta = -a1*sign(r_sta(:,1)).*(abs(r_sta(:,1))).^(b1)-a2*sign(r_sta(:,2)).*(abs(r_sta(:,2))).^(b2)-a3*sign(r_sta(:,3)).*(abs(r_sta(:,3))).^(b3);
U_stc_sta = -a5*sign(Zeta_sta).*(abs(Zeta_sta)).^(1/2) + Theta_sta;
U_sta = U_nom_sta + U_stc_sta;

plot(t1, r_conv(:,1), t2, r_cont(:,1), t3, r_ismc(:,1), t4, r_sta(:,1))
title('Comparison')
xlabel('Time (t) (secs.)')
ylabel('Distance away from plane (r1 = y) (cms.)')
legend('Conventional SMC', 'Continuous SMC', 'ISMC', 'ISMC+STA')
figure
plot(t1, r_conv(:,3), t2, r_cont(:,3), t3, r_ismc(:,3), t4, r_sta(:,3))
title('Comparison')
xlabel('Time (t) (secs.)')
ylabel('r3 State (X-Rotation)')
legend('Conventional SMC', 'Continuous SMC', 'ISMC', 'ISMC+STA')
figure
plot(t1, U_conv, t2, U_cont, t3, U_ismc, t4, U_sta)
title('Comparison')
xlabel('Time (t) (secs.)')
ylabel('Control Input')
legend('Conventional SMC', 'Continuous SMC', 'ISMC', 'ISMC+STA')
hold on
axis image

function xdot = convsmcfun(t,x)
K=6;
eta = 4;
c1=10;
c2=5;
k= 1.0/56.7;
f = -k*sin(acos(-x(3)/k))*x(3)*tan(asin(-x(2)));
g = k*sin(acos(-x(3)/k));
d = 2*sin(t)+3;
zeta = c1*x(1) + c2*x(2) + x(3);
i2 = -1*(c1*x(2) + c2*x(3) + f + K*sign(zeta) + eta*zeta)/g;
xdot(1,1) = x(2);
xdot(2,1) = x(3);
% Without Disturbance
%xdot(3,1) = f + g*i2;
% With Disturbance d=2sint + 3
xdot(3,1) = f + g*i2 + d;

function xdot = contsmcfun(t,x)
k= 1.0/56.7;
f = -k*sin(acos(-x(3)/k))*x(3)*tan(asin(-x(2)));
g = k*sin(acos(-x(3)/k));
d = 2*sin(t)+3;
u_nom = -15*sign(x(1))*(abs(x(1)))^(1/4)-23*sign(x(2))*(abs(x(2)))^(1/3)-9*sign(x(3))*(abs(x(3)))^(1/2);
u = u_nom;
i2 = (u-f)/g;
xdot(1,1) = x(2);
xdot(2,1) = x(3);
% Without Disturbance
%xdot(3,1) = f + g*i2;
% With Disturbance d=2sint + 3
xdot(3,1) = f + g*i2 + d;

function xdot = ismcfun(t,x)
k= 1.0/56.7;
f = -k*sin(acos(-x(3)/k))*x(3)*tan(asin(-x(2)));
g = k*sin(acos(-x(3)/k));
d = 2*sin(t)+3;
a1 = 15;
a2 = 23;
a3 = 9;
b1 = 1/4;
b2 = 1/3;
b3 = 1/2;

u_nom = -a1*sign(x(1))*(abs(x(1)))^(b1)-a2*sign(x(2))*(abs(x(2)))^(b2)-a3*sign(x(3))*(abs(x(3)))^(b3);
zeta = x(4);
a4 = 10;
u_discon = -a4*sign(zeta);
u = u_nom + u_discon;
i2 = (u-f)/g;
xdot(1,1) = x(2);
xdot(2,1) = x(3);
% Without Disturbance
%xdot(3,1) = f + g*i2;
% With Disturbance d=2sint + 3
xdot(3,1) = f + g*i2 + d;
xdot(4,1) = u + d - u_nom;

function xdot = stafun(t,x)
k= 1.0/56.7;
f = -k*sin(acos(-x(3)/k))*x(3)*tan(asin(-x(2)));
g = k*sin(acos(-x(3)/k));
d = 2*sin(t)+3;
a1 = 15;
a2 = 23;
a3 = 9;
b1 = 1/4;
b2 = 1/3;
b3 = 1/2;
UM=6;

u_nom = -a1*sign(x(1))*(abs(x(1)))^(b1)-a2*sign(x(2))*(abs(x(2)))^(b2)-a3*sign(x(3))*(abs(x(3)))^(b3);
zeta = x(4);
theta = x(5);
a5 = 0.13;
a6 = 300;
u_stc = -a5*(zeta^(1/2))*sign(zeta) + theta;
u = u_nom + u_stc;
i2 = (u-f)/g;
xdot(1,1) = x(2);
xdot(2,1) = x(3);
% Without Disturbance
%xdot(3,1) = f + g*i2;
% With Disturbance d=2sint + 3
xdot(3,1) = f + g*i2 + d;
xdot(4,1) = u + d - u_nom;
xdot(5,1) = -a6*sign(zeta);
% If conditioned:
% If you take big enough UM (say 10 or above), it will be same as unconditioned case
%if abs(u_stc) > UM
%   xdot(5,1) = -u_stc;
%else
%   xdot(5,1) = -a6*sign(zeta);
%end