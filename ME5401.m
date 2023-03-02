% initialize parameters
a = 3;
b = 1;
c = 1;
d = 9;
g = 9.8;
M_f = 2.14 + c / 20;
M_r = 5.91 - b/10;
M_c = 1.74;
L_Ff = 0.05;
L_r = 0.128;
L_c = 0.259;
J_x = 0.5+(c-d)/100;
alpha = 15.5-a/3+b/2;
gamma = 11.5+(a-c)/(b+d+3);
H_f = 0.18;
H_r = 0.161;
H_c = 0.098;
L_F = 0.133;
L_R = 0.308+(a-d)/100;
mu_x = 3.33-b/20+a*c/60;
beta = 27.5-d/2;
delta = 60+(a-b)*c/10;
den = M_f*H_f*H_f+M_r*H_r*H_r+M_c*H_c*H_c+J_x;
a_51 = -(M_c * g)/den;
a_52 = (M_f * H_f + M_r*H_r + M_c*H_c)*g/den;
a_53 = ((M_r*L_r*L_F + M_c*L_c*L_F + M_f*L_Ff*L_R)*g)/((L_R + L_F)*den);
a_54 = -M_c*H_c*alpha/den;
a_55 = -mu_x/den;
a_56 = M_f*H_f*L_Ff*gamma/den;
b_51 = M_c*H_c*beta/den;
b_52 = -M_f*H_f*L_Ff*delta/den;
A = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 6.5 -10 -alpha 0 0; a_51 a_52 a_53 a_54 a_55 a_56; 5 -3.6 0 0 0 -gamma];
B = [0 0; 0 0; 0 0; beta 11.2; b_51 b_52; 40 delta];
C = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0];
% 1:
K_1 = unity_place(A,B,[-1 -2 -3 -4 -5 -6]);
% 2:
Q = diag([1 10 20 1 2 3]);
R = 0.1*eye(2);
K_2 = my_lqr(A,B,Q,R);
% 3:
% observability matrix:
O = [C; C*A; C*A^2; C*A^3; C*A^4; C*A^5];
L = place_pole(A',C',[-5 -10 -15 -20 -25 -30])';
% 4:
C_2 = [1 0 0 0 0 0; 0 0 1 0 0 0];
% G(s) is the plant transfer function:
syms s;
G_4 = C_2 * (s * eye(6) - A) ^ (-1) * B;
B_star = [C_2(1,:)*A*B; C_2(2,:)*A*B];
C_star = [C_2(1,:)*(A+2*eye(6))^2; C_2(2,:)*(A+4*eye(6))^2];
K_4 = B_star^(-1) * C_star;
F = B_star^(-1);
s=1;
H_4 = C_2/(s*eye(6) - (A - B*K_4))*B*F;
% 5:
% initialize parameters
a = 3;
b = 1;
c = 1;
d = 9;
Y_sp = -0.1*C*A^(-1)*B*[-0.5 + (a-b)/20; 0.1+(b-c)/(a+d+10)];
% when x is stable and x_dot equals zero, we can get this Y_sp set point.
K_5 = sym('k', [2 9]);
A_5 = [A zeros(6,3); -C zeros(3,3)];
B_5 = [B;zeros(3,2)];
A_5_cl = A_5 - B_5*K_5;
Q=diag([10 9 8 7 6 5 4 3 2]);
Q=3*eye(9);
R=diag([0.1 0.1]);
K_5 = real(my_lqr(A_5, B_5, Q, R));
%6:
% numerical example:
y_sp=[1;2;3];
u = sym('u',[2 1]);
right_equ = -C*(A-B*K_2)^(-1)*B*u;
phi=[-0.0165 0.0774; -0.0580 0.0781; 0.0572 0.0407];
(phi'*phi)^(-1)*phi'*y_sp;