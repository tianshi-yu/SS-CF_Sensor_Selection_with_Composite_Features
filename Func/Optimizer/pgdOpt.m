%-----------------------------------------------%
% Function: solve penalised regression optimisation problem with proximal gradient descent
% Input: 
%       X: between class scatter, n x d
%       Y: label n x 1
%       A: relationship matrix, d x p
%       alpha: scalar parameter balances sensor and feature sparsity
%       lambda: scalar tuning parameter
%       wi: initialised projection vector for warm start
%       varagin: options for the oprimizer
% Output: 
%       sol: struct of solutions
%-----------------------------------------------%
function out = pgdOpt(X,Y,A,alpha,lambda,wi,varargin)
    %---Process options---%
    [maxit,tol,stepInit,verbose,gamma] = ...
    processOptions(varargin{:},'maxit',1000,'tol',1e-3,'step',1,...
    'verbose',false,'gamma',0.1);
    
    %---Constans---%
    d = size(A,1);
    p = size(A,2);
    n = size(X,1);
    
    %---Solve---%
    rng('default')
    J = [];
    if(isempty(wi))
        w = rand(d,1);
    else
        w = wi;
    end
    t = 1;
    wOld1 = w;
    wOld2 = w;
    wNew = w;

    while true % Main optimisation
        
        step = stepInit;
        while true % Optimise step size
            %-Try proximal gradient-%
            for i = 1:p
                I = logical(A(:,i)); % feature index of current sensor group
                %-Gradient update-%
                gd = (1/n).*X(:,I)'* (X*wOld1 - Y);
                r = wOld1(I) - step .* gd;
                
                %-Proximal update-% 
                ps = sqrt(sum(I));
                st =  softTh(r,alpha*step*lambda);
                wNew(I) = plusSgn( 1 - ( ((1-alpha)*step*lambda*ps) / norm(st,2) ) ) * st;
            end
            
            %-Backtracking line search to adjust step size-%
            temp1 = norm(Y-X*wNew,2)^2/(2*n);
            temp2 = norm(Y-X*wOld1,2)^2/(2*n) - ((1/n)*X'*(X*wOld1 - Y))'*(wOld1-wNew) + 0.5*norm(wOld1-wNew,2)^2/step;    
            
            if(temp1>temp2)
                % Shrink the step
                step = gamma*step;
            else
                break;
            end
            
        end
        
        %-Proximal gradient update-%
        w = wNew;
        
        %-Objective function value-%
        J(t,1) = norm(Y-X*w,2)^2/(2*n);
        
        %-Check convergence-%
        temp1 = norm((w-wOld1),'Inf');
        temp2 = norm(w,2);
        
        
        str = ['PGD: step=', num2str(t),', step size=', num2str(step), ', function value=', num2str(J(t,1)), ', rerr=',num2str(temp1/temp2)];
        if(verbose)
           disp(str);
        end
        
        if(temp1<=tol*temp2 || t>maxit)
            break;
        end
        
        %-Increat the step by 1-%
        t = t+1;
        wOld1 = wNew;
        wOld2 = wOld1;
    end
    
    sol.w = w;
    sol.J = J(end,1);
    out = sol;
end







%             G = (wOld - wNew)/step;
%             temp1 = norm(Y-X*(wOld-step*G),2)^2/(2*n);
%             temp2 = norm(Y-X*wOld,2)^2/(2*n) - step*((1/n)*X'*(X*wOld - Y))'*G + 0.5*step*norm(G,2)^2;

%         wOld = w;
%         wNew = w;
%         for i = 1:p
%            I = logical(A(:,i)); % feature index of current sensor group
%            % Gradient update
%            gd = (1/n).*X(:,I)'* (X*wOld - Y);
%            r = wOld(I) - step .* gd;
%            % Proximal update 
%            st =  softTh(r,alpha*step*lambdas);
%            wNew(I) = plusSgn(1 - ( ((1-alpha)*step*lambdas) / norm(st,2)) ) * st;
%         end
%         w = wNew + ((t-2)/(t+1))*(wNew-wOld);
