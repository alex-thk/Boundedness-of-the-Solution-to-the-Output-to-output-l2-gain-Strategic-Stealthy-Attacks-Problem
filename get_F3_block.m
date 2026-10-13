function F_dagger_ij = get_F3_block(i,j,n,A,B)
% Builds [F_{psi(i,j)}]_3 = -Omega' (J_2 kron U_ij) Omega + blkdiag(U_ij, 0_{u x u}), Omega = blkdiag(A,B).
    Omega = blkdiag(A,B);
    last_term = blkdiag(get_Uij(n,i,j), zeros(size(B,2), size(B,2)));
    
    F_dagger_ij = -Omega'*kron(ones(2,2), get_Uij(n,i,j))*Omega + last_term; 
    F_dagger_ij = sparse(F_dagger_ij); 
end

