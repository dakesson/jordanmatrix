function J = jordanmatris(A,tol)

if nargin < 2, tol = eps*1e7; end
[ev,mult] = heltalsev(A,tol)

jordan = zeros(size(A));
progress = 1; %How far the construction has come
for i=1:length(ev) %Loop over eigenvalues
    j=1;
    
    %Initialize
    ri = []; pi = []; bi = []; ni = [];
    while length(ri) < 2 | pi(end) ~= pi(end-1) %Compute table
        ri(j) = rank((A - ev(i)*eye(length(A)))^j);
        pi(j) = length(A) - ri(j);
        
        if j == 1
            bi(j) = pi(j);
        else
            bi(j) = pi(j) - pi(j-1);
        end

        if j > 1
            ni(j-1) = bi(j-1) - bi(j); %Fix, for last pos
        end
        
        j = j+1;
    end
    ni(j) = 0;
    
    positions = find(ni == 1)
    
    for k=1:length(positions)
        n = positions(k);
        jordan(progress:progress+n-1,progress:progress+n-1) = jordanblock(ev(i), n)
        progress = progress + n;
    end
    
    J = jordan
end