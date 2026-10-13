function Uij = get_Uij(n, i, j)
    if i == j
        % e_i e_i^T
        Uij = sparse(i, i, 1, n, n);
    else
        % e_i e_j^T + e_j e_i^T
        Uij = sparse([i j], [j i], [1 1], n, n);
    end
end
