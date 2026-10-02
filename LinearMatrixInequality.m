%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Erfan Rahmani %%%%%%%%%%%
%%%%%%%%%%%% 401123132 %%%%%%%%%%%%%
%%%%%%%%%%%%%%% IIOT %%%%%%%%%%%%%%%
%%%%%%%% Take Home Midterm %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

%% Initials

J = 0.03;
B = 0.1;
kt = 0.5;
n = 2;
N = 4;
A = [0 1;
    0 -B/J];
Ad_tilde = kron(eye(N-1),A);
B = [0;
    kt/J];
Bd = kron(eye(N),B);

%% Computing Kd_tilde

L = [1 0 0 -1;
    0 1 0 -1;
    0 0 1 -1;
    0 0 0 0];

Ec = [1 0 0 0;
    0 1 0 0;
    0 0 1 0];

X = kron(Ec,eye(n));
Lc = kron(L,eye(n));
Lc_tilde = X*Lc;
Pseudo = [eye(6);
    zeros(2,6)];
syms k11 k12 k21 k22 k31 k32 k41 k42

Kd = [k11 k12 0 0 0 0 0 0;
    0 0 k21 k22 0 0 0 0;
    0 0 0 0 k31 k32 0 0;
    0 0 0 0 0 0 k41 k42];

Kd_tilde = Kd*Lc*Pseudo;

%% LMI Solving for P and Kd_tilde

lambda = 1;
e = eig(Ad_tilde - lambda*eye(length(Ad_tilde)));
for i = 1:length(Ad_tilde)
    if e(i) < 0
        % System is stable.
        Stability = 1;
    else
        % System is not stable.
        Stability = 0;
        break
    end
end

cvx_begin sdp
variable P_lambda(6,6) symmetric semidefinite
(Ad_tilde - lambda*eye(6))'*P_lambda + P_lambda*(Ad_tilde - lambda*eye(6)) < -.1*eye(6);
cvx_end

cvx_begin
variables K1(1,2) K2(1,2) K3(1,2) K4(1,2)
KD_Tilde = [K1 0 0 0 0;
    0 0 K2 0 0;
    0 0 0 0 K3;
    0 0 0 0 0 0];
(Ad_tilde + Lc_tilde*Bd*KD_Tilde)'*P_lambda + P_lambda*(Ad_tilde + Lc_tilde*Bd*KD_Tilde) < -.1*eye(6);
cvx_end


disp('Eigenvalues of A + B*K1:')
eig(A + B*K1)
disp('Eigenvalues of A + B*K2:')
eig(A + B*K2)
disp('Eigenvalues of A + B*K3:')
eig(A + B*K3)
disp('Eigenvalues of A + B*K4:')
eig(A + B*K4)
disp('Eigenvalues of Ad_tilde + Lc_tilde*Bd*Kd_tilde:')
eig(Ad_tilde + Lc_tilde*Bd*KD_Tilde)

P_lambda
K1,K2,K3,K4