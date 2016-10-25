function [ev,mult] = heltalsev(A,tol)

if nargin < 2
    tol = eps*1e10;
end
    
[V,D] = eig(A) %V eigenvectors, D eigenvalues

for i=1:size(D,1)
    if isreal(D(i,i))
        if calculate_rest(D(i,i)) > tol
            disp('Not integer')
            return
        end
    else
        if calculate_rest(real(D(i,i))) > tol & calculate_rest(imag(D(i,i))) > tol 
            disp('Imaginary and not integer')            
        end
    end
end

disp('Integer eigenvalues')
D = round(D)
ev = unique(diag(D));
mult = [];
for i=1:length(ev)
    mult(i,1) = length(find(diag(D) == ev(i)));
end
end
