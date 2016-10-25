function block = jordanblock(lambda, n)
    block = zeros(n,n);
    block(1:n+1:n*n) = lambda;
    block(n+1:n+1:n*n) = 1;
end