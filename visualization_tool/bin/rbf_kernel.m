function k = rbf_kernel(gamma, x, y)

    k = exp(-gamma * norm(x-y)^2);

end

