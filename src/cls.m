

function beta = cls(X,z,s_max,Xe)

    [n,k] = size(X);
    [nc,~] = size(Xe);
    m = 2*nc;
    
    eta = 0.95;
    
    G = X'*X;
    g = -X' * z;
    
    % initial conditions
    x = ones(k,1) * 0.01;
    l = ones(m,1) * 0.01;
    s = ones(m,1) * 0.01;
    
    % check if any column is the constant factor
    for i=1:k
        if(all(X(:,i)==X(1,i)))
            x(i) = mean(z)/X(1,i);
            break;
        end
    end
    
    % numerical tolerance
    tol = max(1024, sqrt(n)) * eps;
    
    % initial residuals & complementarity gap
    [rx,rl,rsl,mu] = residuals(G,g,Xe,s_max,x,l,s);
        
    iter = 1;
    while(iter <= 200 && norm(rx) >= tol && norm(rl) >= tol && abs(mu) >= tol)
        
        % we have L*L'=G_bar
        G0 = G + Xe'*diag(l(1:nc)./s(1:nc))*Xe + Xe'*diag(l(nc+1:m)./s(nc+1:m))*Xe;
        L = chol(G0, 'lower');
        [~,dl_a,ds_a] = solve_system(Xe,L,1./s,l,rx,rl,rsl);

        % find affine step size
        alpha_a = min([1.0, min(-l(dl_a<0) ./ dl_a(dl_a < 0)), min(-s(ds_a<0) ./ ds_a(ds_a < 0))]);
        
        % compute affine complementarity measure
        mu_a = (s+alpha_a*ds_a)'*(l+alpha_a*dl_a) / m;
        
        % centering par
        sigma = (mu_a/mu)^3;
        
        % solve modified PC system
        [dx,dl,ds] = solve_system(Xe,L,1./s,l,rx,rl,rsl+ds_a.*dl_a-sigma*mu*ones(m,1));
        
        % compute step size
        alpha = min([1.0, min(-l(dl<0) ./ dl(dl<0)), min(-s(ds<0) ./ ds(ds<0))]);
        
        % update all vars
        x = x + eta*alpha*dx;
        l = l + eta*alpha*dl;
        s = s + eta*alpha*ds;
        
        % compute new residuals & complementarity measure
        [rx,rl,rsl,mu] = residuals(G,g,Xe,s_max,x,l,s);
        
%         fprintf('i = %d n(res)=%g n(rx)=%g n(rl)=%g mu=%g sig=%g cond(G0)=%g\n', ...
%             iter, norm(X*x-z), norm(rx), norm(rl), mu, sigma, cond(G0));
        
%         figure(1);
%         plot(z);
%         hold on
%         plot(X*x, 'r-', 'linewidth', 2);
%         plot([0,n],[2.5,2.5], 'k-');
%         plot([0,n],[0,0], 'k-');
%         hold off
%         legend('z', 'Constr. sol');
%         title(sprintf('Iteration %d', iter));
%         pause;
    
        iter = iter + 1;
    end
    
    beta = x;
    
    
    function [rx,rl,rsl,mu] = residuals(G,g,Xe,s_max,x,l,s)
        rx = G*x + g + Xe'*l(1:nc) - Xe'*l(nc+1:m);
        rl = s - [-Xe*x+s_max;Xe*x];
        rsl = l.*s;
        mu = s'*l/m;
    end
    
    
    function [dx,dl,ds] = solve_system(Xe,L,s_1,l,rx,rl,rsl)
        
        % build the right hand side
        r1 = rsl - l .* rl;
        r0 = -Xe'*(s_1(1:nc).*r1(1:nc))+Xe'*(s_1(nc+1:m).*r1(nc+1:m));

        dx = L'\(L\(-rx - r0));
        Atdx = [-Xe*dx;Xe*dx];
        ds = Atdx - rl;
        dl = -s_1.*(rsl + l.*ds);
    end
        
end
