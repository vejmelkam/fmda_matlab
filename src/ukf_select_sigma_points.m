
% This function selects a sigma ensemble that conforms to the requirements
% for deterministic sampling in the UKF filter.  This method automatically
% generates 2*N points (in addition to the state vector itself), where N is
% the dimension of the state vector.
%
% Sigma point selection according to Simon, 2010, page: 444.
%
%  synopsis: X = ukf_select_sigma_points(x, Sigma, W0)
%
%   x - the current state vector
%   Sigma - the current covariance matrix of the state
%

function X = ukf_select_sigma_points(x, Sigma, kappa)

    if(nargin < 3)
        kappa = 0;
    end
    
    N = size(x,1);
    X = zeros(N, 2*N+1);
    
    % matlab decomposition is A'*A, so we use the rows
    SF = chol((N+kappa)*Sigma);
    for i=1:N
       X(:,i) = x + SF(i,:)';
       X(:,i+N) = x - SF(i,:)';
    end
    
    X(:, 2*N+1) = x;

