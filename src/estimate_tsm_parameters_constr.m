
%
% A trend surface model estimator of regression coefficients and the
% variance of the microscale variability.  An iterative constrained
% Least squares and Fay-Herriot estimator of the variance.
%
%  synopsis: [beta,sigma2]=estimate_tsm_parameters_constr(X,Z,gammas2,Xe,MaxFM)
%  
%
%

function [beta,sigma2]=estimate_tsm_parameters_constr(X,Z,gammas2,Xe,MaxFM)
    [N,k] = size(X);

    sigma2 = 0.0;
    sigma2_old = -10;
    sigma2_old2 = -10;
    iter = 0;
    
    while((abs(sigma2 - sigma2_old)/max(sigma2_old,1e-8) > 1e-3) && ...
          (abs(sigma2 - sigma2_old2)/max(sigma2_old2,1e-8) > 1e-3) && ...
          iter < 20)
        
        sigma2_old2 = sigma2_old;
        sigma2_old = sigma2;

        % transform into model with white unit variance errors and solve
        % for beta
        
        SigSqrt = diag(gammas2.^0.5) + sqrt(sigma2)*eye(N);
        X2 = SigSqrt * X;
        Z2 = SigSqrt * Z;
        beta = cls(X2,Z2,MaxFM,Xe);   

        % compute residuals and re-estimate sigma2
        res2 = (Z - X*beta).^2;
        sigma2 = bisect_solve(res2,gammas2,k);
        
 %       fprintf('beta = [%s] sigma2=%g\n', num2str(beta'), sigma2);
        iter = iter + 1;
    end
    
    if(iter==20)
        beta=nan;
        sigma2=nan;
        warning('tsm_constr: failure to converge in 20 iters');
    end
    
    
