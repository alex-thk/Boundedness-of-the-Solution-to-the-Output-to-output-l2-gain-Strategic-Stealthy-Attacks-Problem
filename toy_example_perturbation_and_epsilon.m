clearvars; close all; clc
rng(1);   % fixed seed so that Fig. 3 is reproducible

% Requires get_Uij.m, get_F3_block.m, psi_map.m on path

%% parameters
eps_val = 1e-3;   % small epsilon used throughout Part 1 (variability test)
N       = 100;    % number of random perturbations per matrix
options = sdpsettings('solver','mosek','verbose',0);

% generate 100 random perturbations in [-1,1] (same draw reused for all three
% perturbation sweeps)
variability = -1 + 2*rand(1,N);

%% fixed (nominal) system data -- shared-zero example from Section V
A  = [2 0; -1 1];
B  = [1; 0];
Cr = [1 1];
Dr = 1;
Cp = [2 2];
Dp = 0;

%% PART 1: variability test (Cp, A, B perturbations), log-ratio boxplot

Gamma_Cp_perturbed = nan(N,2);  % col1 = farkas (Problem 3), col2 = original (Problem 2)
Gamma_A_perturbed  = nan(N,2);
Gamma_B_perturbed  = nan(N,2);

%% perturb Cp
for k = 1:N
    d      = variability(k);
    Cp_pert = [2+d, 2+d];   % perturb *both* entries equally

    [gF, feasF] = solve_problem3_farkas(A, B, Cr, Dr, Cp_pert, Dp, eps_val, options);
    [gO, feasO] = solve_problem2_original(A, B, Cr, Dr, Cp_pert, Dp, options);

    if feasF, Gamma_Cp_perturbed(k,1) = gF; else, warning('Farkas LMI failed at k=%d (Cp perturbation)', k); end
    if feasO, Gamma_Cp_perturbed(k,2) = gO; else, warning('Original LMI failed at k=%d (Cp perturbation)', k); end
end

%% perturb A
for k = 1:N
    d     = variability(k);
    A_pert = [2+d, 0; -1, 1+d];   % perturb *both* entries equally

    [gF, feasF] = solve_problem3_farkas(A_pert, B, Cr, Dr, Cp, Dp, eps_val, options);
    [gO, feasO] = solve_problem2_original(A_pert, B, Cr, Dr, Cp, Dp, options);

    if feasF, Gamma_A_perturbed(k,1) = gF; else, warning('Farkas LMI failed at k=%d (A perturbation)', k); end
    if feasO, Gamma_A_perturbed(k,2) = gO; else, warning('Original LMI failed at k=%d (A perturbation)', k); end
end

%% perturb B
for k = 1:N
    d      = variability(k);
    B_pert = [1+d; 0];

    [gF, feasF] = solve_problem3_farkas(A, B_pert, Cr, Dr, Cp, Dp, eps_val, options);
    [gO, feasO] = solve_problem2_original(A, B_pert, Cr, Dr, Cp, Dp, options);

    if feasF, Gamma_B_perturbed(k,1) = gF; else, warning('Farkas LMI failed at k=%d (B perturbation)', k); end
    if feasO, Gamma_B_perturbed(k,2) = gO; else, warning('Original LMI failed at k=%d (B perturbation)', k); end
end

%% Compute log-ratio
logratio_Cp = log( Gamma_Cp_perturbed(:,1) ./ Gamma_Cp_perturbed(:,2) );
logratio_A  = log( Gamma_A_perturbed(:,1)  ./ Gamma_A_perturbed(:,2) );
logratio_B  = log( Gamma_B_perturbed(:,1)  ./ Gamma_B_perturbed(:,2) );

T = table(logratio_A, logratio_B, logratio_Cp, ...
    'VariableNames', {'logratio_A', 'logratio_B', 'logratio_Cp'});
data_dir = fullfile(fileparts(mfilename('fullpath')), 'data');
if ~exist(data_dir, 'dir'), mkdir(data_dir); end
writetable(T, fullfile(data_dir, 'logratio_values.csv'));

vars = T.Properties.VariableNames;
for k = 1:numel(vars)
    x = T.(vars{k});
    x = x(~isnan(x));

    pr      = prctile(x, [25 50 75]);
    Q1 = pr(1); Med = pr(2); Q3 = pr(3);
    IQR     = Q3 - Q1;
    whiskLo = min(x(x >= Q1 - 1.5*IQR));
    whiskHi = max(x(x <= Q3 + 1.5*IQR));

    fprintf('%s: %.6g, %.6g, %.6g, %.6g, %.6g\n', vars{k}, Q1, Med, Q3, whiskLo, whiskHi);
end

media = mean([logratio_Cp; logratio_A; logratio_B]);
stdev = std([logratio_Cp; logratio_A; logratio_B]);

%% decreasing-epsilon sequence test
% Restate the nominal system explicitly (A, B, Cp were overwritten by the
% perturbation loops above).
A  = [2 0; -1 1];
B  = [1; 0];
Cr = [1 1];
Dr = 1;
Cp = [2 2];
Dp = 0;

eps_seq = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6];

gamma_farkas_seq = nan(size(eps_seq));
feasible_seq     = false(size(eps_seq));

for kk = 1:numel(eps_seq)
    [gamma_farkas_seq(kk), feasible_seq(kk)] = ...
        solve_problem3_farkas(A, B, Cr, Dr, Cp, Dp, eps_seq(kk), options);
    if ~feasible_seq(kk)
        warning('Problem 3 infeasible/solver failed at eps = %.g', eps_seq(kk));
    end
end

% reference value: original (unrelaxed) Problem 2, eps-free
[gamma_ref, feasible_ref] = solve_problem2_original(A, B, Cr, Dr, Cp, Dp, options);
if ~feasible_ref
    warning('Original Problem 2 solve failed on the nominal system.');
end

abs_error_seq = abs(gamma_farkas_seq - gamma_ref);

T_eps = table(eps_seq(:), gamma_farkas_seq(:), feasible_seq(:), abs_error_seq(:), ...
    'VariableNames', {'eps_k','gamma_star_farkas','feasible','abs_error_vs_original'});

disp(T_eps)
writetable(T_eps, fullfile(data_dir, 'epsilon_sequence_test.csv'));


%% Local functions
function [gamma_val, feasible] = solve_problem3_farkas(A, B, Cr, Dr, Cp, Dp, eps_val, options)
% solve prob 3 (Farkas/SDP reformulation) at a given eps_val, using the
% same F_h construction 
    n = size(A,1);
    u = size(B,2);

    yalmip('clear');
    P     = sdpvar(n,n,'symmetric');
    gamma = sdpvar(1);
    upper_mask = triu(true(n));
    zeta = [P(upper_mask); gamma];

    F_matrices = cell(psi_map(n,n,n)+2,1);
    for i = 1:n
        for j = i:n
            idx = psi_map(i,j,n);
            Uij = get_Uij(n,i,j);
            F3_ij = get_F3_block(i,j,n,A,B);
            F_matrices{idx} = sparse([Uij,          zeros(n,1),   zeros(n,n+u); ...
                                       zeros(1,n),   0,            zeros(1,n+u); ...
                                       zeros(n+u,n), zeros(n+u,1), F3_ij]);
        end
    end

    % F_{psi(n,n)+1}: coefficient of gamma
    F_matrices{end-1} = sparse([zeros(n,n),   zeros(n,1),   zeros(n,n+u); ...
                                 zeros(1,n),   1,            zeros(1,n+u); ...
                                 zeros(n+u,n), zeros(n+u,1), ([Cr Dr]'*[Cr Dr])]);

    % F0(eps_val) = blkdiag(-eps*I_n, 0, Hp'Hp - eps*I_{n+u})
    F_matrices{end} = sparse([-eps_val*eye(n), zeros(n,1),   zeros(n,n+u); ...
                               zeros(1,n),      0,            zeros(1,n+u); ...
                               zeros(n+u,n),    zeros(n+u,1), ([Cp Dp]'*[Cp Dp]) - eps_val*eye(n+u)]);

    Fsum = 0;
    m = length(zeta);
    for i = 1:m-1
        Fsum = Fsum + zeta(i)*F_matrices{i};
    end
    Fsum = Fsum + zeta(end)*F_matrices{end-1};

    Constraints = (Fsum - F_matrices{end}) >= 0;
    sol = optimize(Constraints, gamma, options);

    if sol.problem == 0
        gamma_val = double(gamma);
        feasible  = true;
    else
        gamma_val = NaN;
        feasible  = false;
    end
end

function [gamma_val, feasible] = solve_problem2_original(A, B, Cr, Dr, Cp, Dp, options)
% solve original(unrelaxed) Problem 2 LMI 
    n = size(A,1);

    yalmip('clear');
    P2     = sdpvar(n,n,'symmetric');
    gamma2 = sdpvar(1);
    R = [A'*P2*A - P2, A'*P2*B; B'*P2*A, B'*P2*B];

    LMI2 = -R + gamma2*([Cr,Dr]'*[Cr,Dr]) - ([Cp,Dp]'*[Cp,Dp]) >= 0;
    Constraints2 = [LMI2, P2 >= 0, gamma2 >= 0];

    sol = optimize(Constraints2, gamma2, options);

    if sol.problem == 0
        gamma_val = double(gamma2);
        feasible  = true;
    else
        gamma_val = NaN;
        feasible  = false;
    end
end
