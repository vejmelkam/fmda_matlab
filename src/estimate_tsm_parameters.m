
%
% A trend surface model estimator of regression coefficients and the
% variance of the microscale variability.  An iterative Least squares
% and Fay-Herriot estimator of the variance.
%
%  synopsis: [beta,sigma2]=estimate_tsm_parameters(X,Z,gammas2)
%  
%
%

function [beta,sigma2]=estimate_tsm_parameters(X,Z,gammas2)
    [N,k] = size(X);

    sigma2 = 0.0;
    sigma2_old = -10;
    
    while(abs(sigma2 - sigma2_old)/max(sigma2_old,1e-8) > 1e-3)
        
        sigma2_old = sigma2;

        % transform into model with white unit variance errors and solve
        % for beta
        SigSqrt = diag(gammas2.^0.5) + sqrt(sigma2)*eye(N);
        X2 = SigSqrt * X;
        Z2 = SigSqrt * Z;
        [Q, R] = qr(X2);
        beta = R \ Q'*Z2;

        % compute residuals and re-estimate sigma2
        res2 = (Z - X*beta).^2;
        sigma2 = bisect_solve(res2,gammas2,k);
        
    end
    
