% declare symbolic variables.
syms theta_1(t) theta_2(t) L_1 L_2 m_1 m_2 g;

% declare the equations of motion.
% (rearranged for cleaner code)
w_1 = diff(theta_1);
w_2 = diff(theta_2);
alpha_1 = (L_2/L_1)*(m_2/(m_1 + m_2))*cos(theta_1 - theta_2);
alpha_2 = (L_1/L_2)*cos(theta_1-theta_2);
f_1 = -(L_2/L_1)*(m_2/(m_1 + m_2))*diff(theta_2)^2*sin(theta_1-theta_2)-(g/L_1)*sin(theta_1);
f_2 = (L_1/L_2)*diff(theta_1)^2*sin(theta_1-theta_2)-(g/L_2)*sin(theta_2);

% actual equations of motion,
% built from previously declared equations.
symsEqn_1 = diff(w_1) + alpha_1*diff(w_2) == f_1;
symsEqn_2 = diff(w_2) + alpha_2*diff(w_1) == f_2;

% substitute symbolic variables with real values.
eqn_1 = subs(symsEqn_1,[L_1 L_2 m_1 m_2 g],[1 1 1 1 9.8]);
eqn_2 = subs(symsEqn_2,[L_1 L_2 m_1 m_2 g],[1 1 1 1 9.8]);

% reduce to first-order differential equations,
% returned as a symbolic vector.
[V,S] = odeToVectorField(eqn_1,eqn_2);

% convert to MATLAB function
M = matlabFunction(V,'vars',{'t','Y'});

% replace theta1 and theta2 with real values.
% WARNING: Code won't work unless real values
% are substituted in.
% eg: initCond = [0 0 pi/12 0];
initCond = [theta2 0 theta1 0];

% solves the differential equation, t spanning
% from 0 to 20.
sols = ode113(M,[0 20],initCond);

% extracting theta1 (column 3) and theta2 (column 1)
% from the solution, to plot only theta1 and theta2
solstheta = sols.y([3 1],:);

% plotting the solution
plot(sols.x,solstheta)
legend('\theta_1','\theta_2')
title('Solutions of State Variables')
xlabel('Time (s)')
ylabel('Solutions (rad)')