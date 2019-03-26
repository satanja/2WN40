iter = 1000;
tol = 10^-6;
%% 10a

A = [5 4 1 1; 4 5 1 1; 1 1 4 2; 1 1 2 4];

x0 = ones(4, 1);
x0 = x0 / norm(x0);
[x, lambda, conv] = rayleighquotients(A, x0, tol, iter)

%% 10 b


for n = 1 : 10
   A1 = diag(2 * ones(n + 1, 1));
   A2 = diag(ones(n, 1), 1);
   A3 = diag(ones(n, 1), -1);
   A = A1 + A2 + A3;
   x0 = ones(n + 1, 1);
   [x, lambda, conv] = rayleighquotients(A, x0, tol, iter)
end


function [x, lambda, conv] = rayleighquotients(A, x0, tol, maxint)
x = x0;
lambda = 0;
conv = 0;

lambda1 = 0; % \mu_{k - 1}
lambda2 = 0; % \mu_{k - 2}

    for k = 0 : maxint
        y = A * x;
        
        %shift last computed approximations and compute new approximation
        lambda2 = lambda1; 
        lambda1 = lambda;   
        lambda = dot(y, x);
        
        x = y / norm(y);  
        
        delta = (lambda - lambda1) / (lambda1 - lambda2); % \Delta_k
        ltilde = lambda + delta / (1 - delta) * (lambda - lambda1); % \tilde{\mu_k}}
        cond = abs((ltilde - lambda) / lambda);
        if (k >= 2 && cond < tol)
            conv = 1;
            break;
        end
    end
end