function continuous_smc

tspan = [0 10];
r0 = [0.3; 0.5; 0.2];

[t,r] = ode45(@myodefun,tspan,r0);
K=6;
eta = 4;
c1=10;
c2=5;
k= 1.0/56.7;
Zeta = c1*r(:,1) + c2*r(:,2) + r(:,3);
U = -15*sign(r(:,1)).*(abs(r(:,1))).^(1/4)-23*sign(r(:,2)).*(abs(r(:,2))).^(1/3)-9*sign(r(:,3)).*(abs(r(:,3))).^(1/2);
f = -k*sin(acos(-r(:,3)/k)).*r(:,3).*tan(asin(-r(:,2)));
g = k*sin(acos(-r(:,3)/k));
i2 = g\(U-f);
plot(t, r(:,1), 'b')
title('Continuous SMC With Disturbance')
xlabel('Time (t) (secs.)')
ylabel('Distance away from plane (r1 = y) (cms.)')
figure
plot(t, r(:,2), 'r', t, r(:,3), 'b')
title('Continuous SMC With Disturbance')
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
title('Continuous SMC With Disturbance')
xlabel('Time (t) (secs.)')
ylabel('Sliding Surface')
figure
plot(t, U, 'b')
title('Continuous SMC With Disturbance')
xlabel('Time (t) (secs.)')
ylabel('Control Input')
hold on
axis image

function xdot = myodefun(t,x)
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