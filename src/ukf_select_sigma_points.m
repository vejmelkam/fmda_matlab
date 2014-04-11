
% This function selects a sigma ensemble that conforms to the requirements
% for deterministic sampling in the UKF filter.  This method automatically
% generates 2*N points (in addition to the state vector itself), where N is
% the dimension of the state vector.
%
% Exception: the space for the observation noise is accoutned for but
% observation noise does not affect the development of the forecast.  This
% function relies on ukf_update to fill in the observation noise sigma_f
% part.
%
%  synopsis: X = ukf_select_sigma_points(x,Pf,Q,nv,kappa)
%
%   x  - the current state vector (of original system)
%   Pf - the current covariance matrix of the state
%   Q  - the process noise covariance (must be diagonal)
%   nv - dimension of observation noise
%

function X = ukf_select_sigma_points(x,Pf,Q,nv,kappa)

    if(nargin < 5)
        kappa = 0;
    end
    
    nx = size(x,1);
    n = 2*nx + nv;
    
    % augmented state & space for sigma points
    xa = [x;zeros(nx+nv,1)];
    X = zeros(n, 2*n+1);
    
    % construct Cholesky factor of augmented covariance
    SF = zeros(n);
    SF(1:nx,1:nx) = chol((n+kappa)*Pf);
    SF(nx+1:2*nx,nx+1:2*nx) = sqrt(n+kappa)*diag(sqrt(diag(Q)));
    
    % matlab decomposition is A'*A, so we use the rows
    for i=1:n
       X(:,i) = xa + SF(i,:)';
       X(:,i+n) = xa - SF(i,:)';
    end
    
    X(:, 2*n+1) = xa;
