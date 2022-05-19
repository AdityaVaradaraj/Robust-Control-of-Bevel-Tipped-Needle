function ismc_sta

tspan = [0 10];
r0 = [0.3; 0.5; 0.2; 0; 0];

[t,r] = ode45(@myodefun,tspan,r0);
Zeta = r(:,4);
Theta = r(:, 5);
a1 = 15;
a2 = 23;
a3 = 9;
b1 = 1/4;
b2 = 1/3;
b3 = 1/2;
a5 = 0.13;

U_nom = -a1*sign(r(:,1)).*(abs(r(:,1))).^(b1)-a2*sign(r(:,2)).*(abs(r(:,2))).^(b2)-a3*sign(r(:,3)).*(abs(r(:,3))).^(b3);
U_stc = -a5*sign(Zeta).*(abs(Zeta)).^(1/2) + Theta;
U = U_nom + U_stc;
plot(t, r(:,1), 'b')
title('ISMC + STA With Disturbance')
xlabel('Time (t) (secs.)')
ylabel('Distance away from plane (r1 = y) (cms.)')
figure
plot(t, r(:,2), 'r', t, r(:,3), 'b')
title('ISMC + STA With Disturbance')
xlabel('Time (t) (secs.)')
ylabel('States denoting rotations (r2 and r3)')
legend('r2 state', 'r3 state')
figure
plot3(r(:,1), r(:,2), r(:,3))
xlabel('r1')
ylabel('r2')
zlabel('r3')
figure
plot(t, Zeta, 'b')
title('ISMC + STA With Disturbance')
xlabel('Time (t) (secs.)')
ylabel('Sliding Surface')
figure
plot(t, U, 'b')
title('ISMC + STA With Disturbance')
xlabel('Time (t) (secs.)')
ylabel('Control Input')
figure
plot(t, U_nom, 'b')
title('ISMC + STA With Disturbance')
xlabel('Time (t) (secs.)')
ylabel('Nominal Control Input')
hold on
axis image

function xdot = myodefun(t,x)
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
UM = 6;

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
