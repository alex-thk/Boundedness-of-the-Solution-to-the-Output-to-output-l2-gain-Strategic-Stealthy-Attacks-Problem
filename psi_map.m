function [vec_idx] = psi_map(i,j,n)
    vec_idx = n*(i-1) - (((i-2)*(i-1))/2) + j - i + 1; 
end

