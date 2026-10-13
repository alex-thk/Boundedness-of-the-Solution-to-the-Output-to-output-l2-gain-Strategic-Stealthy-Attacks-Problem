clear all; 
close all; 
clc;

eps = 1e-3; 

%% system matrices
Ap = [ ...
   -0.1068,   0,       0.0275,    0; 
    0,       -0.0903,  0,         0.0258;
    0,        0,      -0.0275,    0;
    0,        0,       0,        -0.0258 ];
Bp = [ ...
    0.0802,     0; 
    0,         0.0807;
    0,         0.1345;
    0.1337,    0 ];
Cm = [ ...
    0.2000,  0,      0,      0;
    0,       0.2000, 0,      0 ];
Cp = [ ...
    6.3246,  0,      0,      0;
    0,       6.3246, 0,      0;
    0,       0,      4.4721, 0;
    0,       0,      0,      4.4721;
    0,       0,      0,      0;
    0,       0,      0,      0 ];
Dp = [ zeros(4,2); 3.1623*eye(2) ];
Gamma_u = eye(2);  % from springer chap. 

%% LQE KALMAN FILTER AS OBSERVER 
% Kf -> kalman gain matrix 
% eigs_kalman_filter -> eigs of (A-Kf*C)
Qn = 1e-2 * eye(4);   % process‐noise covariance
Rn = 1e-1 * eye(2);   % measurement‐noise covariance
[Kf, ~, eigs_kalman_filter] = lqe(Ap, eye(4), Cm, Qn, Rn);

%% LQ REGULATOR WITH STATE ESTIMATE FROM LQE 
Q = Cp'*Cp; 
R = Dp'*Dp; 
sys = ss(Ap, Bp, Cp, Dp); 
[L, ~, ~] = lqr(sys, Q, R); 

%% Closed loop system matrices 
% Case where E_p = 0.
Gamma_u = eye(2); 

A = [Ap-Bp*L,  Bp*L;
           zeros(4),  Ap-Kf*Cm]; 
B = [Bp*Gamma_u; 
     Bp*Gamma_u];

Cp_cl = [Cp-Dp*L,  Dp*L]; 
Cr = [zeros(2,4),  Cm]; 

Dp_cl = Dp*Gamma_u; 
Dr = zeros(2,2); 


sysp = ss(A,B,Cp_cl,Dp_cl); 
sysr = ss(A,B,Cr,Dr);

sysr_discrete = c2d (sysr, 1);
tzero(sysr_discrete);

zeros_sysp = tzero(sysp); 
zeros_sysr = tzero(sysr); 

%% simulations in DISCRETE TIME 
% Convert to discrete time
Ts = 1;
sysp_d = c2d(sysp, Ts);
sysr_d = c2d(sysr, Ts);

% Get unstable zero of detector system (discrete)
zeros_sysr = tzero(sysr_d);
idx = find(abs(zeros_sysr) > 1 + 1e-3, 1); % unstable discrete zero
z = zeros_sysr(idx);

% Compute attack direction in discrete time
A_d = sysr_d.A;
B_d = sysr_d.B;
C_d = sysr_d.C;
D_d = sysr_d.D;

Rosenbrock = [z*eye(size(A_d)) - A_d, -B_d;
              C_d,                   D_d];
nv = null(Rosenbrock); % Null space of Rosenbrock matrix
x_non_normalized =  nv(1:end-2,1); 
x = x_non_normalized / norm(x_non_normalized); 
v_non_normalized = nv(end-1:end,1);   % Input direction
v = v_non_normalized / norm(v_non_normalized);       % Normalize

%% certificate Y from Theorem 3, Eq. (17)
Y_blk3 = nv * nv'; 
Y_blk2 = 0; 
Y_blk1 = x_non_normalized * (abs(z)^2 -1) * x_non_normalized'; 

Y = blkdiag(Y_blk1, Y_blk2, Y_blk3); 

%% check our system theoretical Farkas Lemma conditions on the Y like in Lemma
Cp_cl = sysp_d.C;
Dp_cl = sysp_d.D;
Cr = sysr_d.C; 
Dr = sysr_d.D; 
n = 8;
u = 2;
H = [eye(n), zeros(n,1), [eye(n), zeros(n,u)]];  
G = [zeros(n), zeros(n,1), [A_d, B_d]];
W = [zeros(u,n), zeros(u,1), [Cr, Dr]];
S = [zeros(n-u,n), zeros(n-u,1), [Cp_cl, Dp_cl]];
I = blkdiag(eye(n), 0, eye(n+u)); 

norm(H*Y*H' - G*Y*G', 'fro') % iva 
trace((W'*W + eye(n+1+n+u) - I)*Y) % ivb 
trace((S'*S-eps*I)*Y) % v

%% Simulation
Ts = 1;
Tsim = 1000;
Ns = Tsim / Ts;
t = (0:Ns-1) * Ts; % Discrete time vector

% Attack signal in discrete time: a[k] = z^k * v
a_t = (z .^ t)' * v.';
             
% Simulate using discrete systems
[y_c,~,~] = lsim(sysp_d, a_t, t);
[y_r,~,~] = lsim(sysr_d, a_t, t);

% Compute norms
yc_norm = vecnorm(y_c, 2, 2);
yr_norm = vecnorm(y_r, 2, 2);
a_norm  = vecnorm(a_t,  2, 2);

% Cumulative output energies (Fig. 2): sum_{t<=k} ||y(t)||^2
yc_cum_energy = cumsum(sum(y_c.^2,2));
yr_cum_energy = cumsum(sum(y_r.^2,2));


% saving signals for pgfplots 
T = table(t', ...
    y_c(:,1), y_c(:,2), y_c(:,3), y_c(:,4), y_c(:,5), y_c(:,6), ...
    y_r(:,1), y_r(:,2), ...
    a_t(:,1), a_t(:,2), ...
    yc_norm, yr_norm, a_norm, ...
    yc_cum_energy, yr_cum_energy, ...
    'VariableNames', {'Time', ...
                      'y_c_1', 'y_c_2', 'y_c_3', 'y_c_4', 'y_c_5', 'y_c_6', ...
                      'y_r_1', 'y_r_2', ...
                      'a_1', 'a_2', ...
                      'yc_norm', 'yr_norm', 'a_norm', ...
                      'yc_cum_energy', 'yr_cum_energy'});

% Save to CSV (data/ folder)
data_dir = fullfile(fileparts(mfilename('fullpath')), 'data');
if ~exist(data_dir, 'dir'), mkdir(data_dir); end
writetable(T, fullfile(data_dir, 'full_signal_output_NEW.csv'));


%% FARKAS PROBLEM (FIND Y)
yalmip('clear');
n = size(A,1); 
u = size(B,2); 

F_matrices = cell((psi_map(n,n,n)+2),1); % last referred to gamma and F0

u = size(B,2); 
for i = 1:n 
    for j = i:n 
        idx = psi_map(i,j,n); 
        Uij = get_Uij(n,i,j);
        F_dagger_ij = get_F3_block(i,j,n,A_d,B_d); 

        F_matrices{idx} = sparse([Uij, zeros(n,1), zeros(n,n+u); ...
                                   zeros(1,n), 0, zeros(1,n+u); ...
                                   zeros(n+u,n), zeros(n+u,1), F_dagger_ij]); 
    end 
end 

 F_matrices{end-1} = sparse([zeros(n,n), zeros(n,1), zeros(n,n+u); ...
                    zeros(1,n), 1, zeros(1,n+u); ...
                    zeros(n+u,n), zeros(n+u,1), ([Cr Dr]'*[Cr Dr])]);
 
 F_matrices{end} = sparse([-eps*eye(n), zeros(n,1), zeros(n,n+u); ...
                    zeros(1,n), 0, zeros(1,n+u); ...
                    zeros(n+u,n), zeros(n+u,1), ([Cp_cl Dp_cl]'*[Cp_cl Dp_cl])-(eps*eye(n+u))]); %F0 

%% YALMIP setup
Y_var = sdpvar(n+n+u+1, n+n+u+1, 'symmetric');
warmstart(Y_var, Y); 
Constraints = [Y_var >= 0, Y_var == Y_var', trace(F_matrices{end}*Y_var) >= 0];
m = (n*(n+1)/2)+1;
for i = 1:m 
    Constraints = [Constraints, trace(F_matrices{i} * Y_var) == 0];
end
options = sdpsettings('solver', 'mosek','verbose',0);
obj = norm(Y_var - Y,'fro');
sol_dual = optimize(Constraints, obj, options);
% check(Constraints); 

if sol_dual.problem == 0
    disp('Solver converged successfully!');
    disp('Optimal Y:'), disp(value(Y_var));
    disp(sol_dual.info);
    Yv = value(Y_var); 
else
    warning('Solver did not converge. Info: %s', sol_dual.info);
end

% 1) Absolute Frobenius‐norm error
resid = norm(Yv - Y, 'fro');
fprintf('Absolute Frobenius‐norm ||Yv − Y||_F = %.3e\n', resid);

%% PRIMAL PROBLEM (ANC. FORMULATION) 
P = sdpvar(n, n, 'symmetric'); % Symmetric matrix P
gamma = sdpvar(1);             % Scalar variable gamma
upper_triangle_mask = triu(true(n)); % n-by-n logical mask selecting the upper triangle
zeta = [P(upper_triangle_mask); gamma];  

Fsum = 0;
for k = 1:(m-1)
    Fsum = Fsum + zeta(k)*F_matrices{k};
    % disp('Added index to contraint')
    % disp(k)
end
% zeta(m) is gamma -> multiply F_matrices{end}
Fsum = Fsum + zeta(m)*F_matrices{end-1}; 

Constraints = Fsum - F_matrices{end} >= 0;

Objective = zeta(m);  % gamma (o2o-gain)
options = sdpsettings('solver','mosek','verbose',0);

sol_primal = optimize(Constraints, Objective, options);

if sol_primal.problem == 0
    disp('Solver converged successfully!');
    disp('Optimal gamma:'), disp(value(zeta(m)));
    disp('Optimal P:'), disp(value(zeta(1:end-1))); 
    disp(sol_primal.info);
else
    warning('Solver did not converge. Info: %s', sol_primal.info);
end

%% check our lemma conditions with Yv and verify that they are of course satisfied
H = [eye(n), zeros(n,1), [eye(n), zeros(n,u)]];  
G = [zeros(n), zeros(n,1), [A_d, B_d]];
W = [zeros(u,n), zeros(u,1), [Cr, Dr]];
S = [zeros(n-u,n), zeros(n-u,1), [Cp_cl, Dp_cl]];
I = blkdiag(eye(n), 0, eye(n+u)); 

norm(H*Yv*H' - G*Yv*G', 'fro') % iva 
trace((W'*W + eye(n+1+n+u) - I)*Yv) % ivb 
trace((S'*S-eps*I)*Yv) % v
