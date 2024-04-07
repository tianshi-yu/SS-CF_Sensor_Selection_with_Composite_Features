%-----------------------------------------------%
% Function: solve convex form of separability problem with penalty
% Input: 
%       X: between class scatter, n x d
%       Y: label n x 1
%       A: relationship matrix, d x p
%       m: projected subspace dimension
%       alpha: scalar parameter balances sensor and feature sparsity
%       lambda: scalar tuning parameter
%       lamFold: cv fold number for lambda tuning
%       Wi: initialised projection matrix for warm start
%       tol: tolerance of optimal scoring optimisation 
%       maxit: maximum iterations of optimal scoring optimisation
%       varargin: options for the optimiser
% Output: 
%       W: optimal sparse projection matrix
%       Theta: optimal scoring
%       mse: least square objective function value
%-----------------------------------------------%


function [W, Theta, mse] = optScoreReg(X,Y,A,m,alpha,lambda,lamFold,Wi,tol,maxit,varargin)
    rng('default')
    d = size(X,2);
    c = length(unique(Y));

    %% Cross-validation for lambda tuning - if asked
    if(~isempty(lamFold))
        % CV partition
        cv = cvpartition(Y,'KFold',lamFold,'Stratify',true); 
        mse = [];
        for i = 1:lamFold
            % Preparation
            trainIdx = gather(cv.training(i));
            testIdx = gather(cv.test(i));
            Xt = X(trainIdx,:);
            Yt = Y(trainIdx,:); Yt = dummyVarLabel(Yt);
            Xv  = X(testIdx,:);
            Yv  = Y(testIdx,:); Yv = dummyVarLabel(Yv);
            n = size(Xt,1);
            D = (Yt'*Yt)./n;
            
            W = zeros(d,m);
            Theta = eye(c,m);
                % Main loop until converge
                for k = 1:m
                    wNorm = [];
                    t = 0; err = tol + 1;
                    while true

                        % Solve for W using pgd solver
                        if(i==1)
                            fit = pgdOpt(Xt,Yt*Theta(:,k),A,alpha,lambda,Wi(:,k),varargin{:});
                        else
                            fit = pgdOpt(Xt,Yt*Theta(:,k),A,alpha,lambda,Wt(:,k),varargin{:});
                        end
                            
                        
                        w = fit.w(:,1);

                        % Update scoring Theta
                        theta_hat = D\Yt'*Xt*w;
                        if(k~=1)
                            theta_hat = (eye(c,c) - Theta(:,1:k-1)*Theta(:,1:k-1)'*D ) * theta_hat;
                        else
                            theta_hat = (eye(c,c) - ones(c,1)*ones(c,1)'*D ) * theta_hat;
                        end

                        if(sum(abs(theta_hat)) ~= 0)
                            Theta(:,k) = theta_hat./(sqrt( theta_hat'*D*theta_hat ));
                        else
                            Theta(:,k) = theta_hat;
                        end


                        wNorm = [wNorm;norm(w,2)/size(w,1)]; 
                        if(  checkConvergence(wNorm,tol) || t>maxit)
                            break;
                        end

                        wOld = w;
                        t = t + 1;
                    end

                    W(:,k) = w;
                    
                end
                Wt = W;
                
                % Validation
                if(~all(all(W==0)))
                    mse(i) = norm(Yv*Theta-Xv*W,'fro')^2/(2*size(Xv,1));
                else
                    mse(i) = NaN;
                end
            
        end  
    end
    %% Train the model with the Lambda
    
    % Preparations
    n = size(X,1);
    Y = dummyVarLabel(Y);
    W = zeros(d,m);
    w0 = zeros(d,1);
    Theta = eye(c,m);
    D = (Y'*Y)./n;

    % Main loop until converge
    for k = 1:m
        
        wNorm = [];
        t = 0; err = tol + 1;
        while true
            
            % Solve for W using pgd solver
            fit = pgdOpt(X,Y*Theta(:,k),A,alpha,lambda,Wt(:,k),varargin{:});
            w = fit.w(:,1);
               
            % Update scoring Theta
            theta_hat = D\Y'*X*w;
            if(k~=1)
                theta_hat = (eye(c,c) - Theta(:,1:k-1)*Theta(:,1:k-1)'*D ) * theta_hat;
            else
                theta_hat = (eye(c,c) - ones(c,1)*ones(c,1)'*D ) * theta_hat;
            end
            
            if(sum(abs(theta_hat)) ~= 0)
                Theta(:,k) = theta_hat./(sqrt( theta_hat'*D*theta_hat ));
            else
                Theta(:,k) = theta_hat;
            end
            

            wNorm = [wNorm;norm(w,2)/size(w,1)]; 
            if(  checkConvergence(wNorm,tol) || t>maxit)
                break;
            end
            
            wOld = w;
            t = t + 1;
        end
        
        W(:,k) = w;
       
    end
    
    %% Validation
    if(~isempty(lamFold))
        % For cross validation we already have the mse value
        if(~all(isnan(mse)))
            mse = mse(~isnan(mse));
        end
        mse = mean(mse);
    else
        if(~all(all(W==0)))
            mse = norm(Y*Theta-X*W,'fro')^2/(2*size(X,1));
        else
            mse = NaN;
        end
    end
end

