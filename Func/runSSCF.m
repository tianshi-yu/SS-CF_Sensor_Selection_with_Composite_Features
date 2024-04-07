%-----------------------------------------------%
% Function: Sensor Selection with Composite Features (SS-CF)
% Input: 
%       X: data table, n x d
%       Y: label n x 1
%       G: candidate sensors names p x 1
%       A: relationship matrix, d x p
%       L: feature number constraint
%       q: sensor number constraint
%       m: projected subspace dimension
%       alpha: scalar parameter balances sensor and feature sparsity
%       lamPath: search grid of the tuning parameter lambda 
%       lamFold: cv fold number for lambda tuning
%       tol: tolerance of  optimal scoring optimisation 
%       maxIte: maximum iterations of optimal scoring optimisation
%       varargin: options for the optimiser
% Output: 
%       ssN: sensor selection results return the names
%       ss: sensor selection results return the indices
%       fs: feature selection results return the indices
%       W: sparse projection matrix
%       J: final objective function value
%-----------------------------------------------%
function [ssN,ss,fs,W,mse] = runSSCF(X,Y,G,A,L,q,m,alpha,lamPath,lamFold,tol,maxit,verbose,dispPlot,varargin)
    rng('default')
    %---Initialisation---%
    ssN = []; ss = []; fs = []; W = [];
    tooSparse = true;
    Wt = cell(length(lamPath),1); Wtc = cell(length(lamPath),1);
    Theta = cell(length(lamPath),1); 
    mse = NaN(length(lamPath),1);
    ssNum = zeros(length(lamPath),1);

    %PGD and solution path
    for i = 1:length(lamPath)
        if(verbose)
            disp(strcat('*---Solution path: Lambda ', num2str(i),' ---%'));
        end
        optOption = varargin{:};

        Xin = table2array(X); 
        [n1,m1] = find(Y ==unique(Y)');  Yin(n1,1) = m1; 
        Ain = A;
        if(i==1)
            Wi = rand(size(Xin,2),m);
            [Wt{i,1}, Theta{i,1}, mse(i,1)] = optScoreReg(Xin,Yin,Ain,m,alpha,lamPath(i),lamFold,Wi,tol,maxit,optOption);
        else
            Wi = Wt{i-1,1};
            [Wt{i,1}, Theta{i,1}, mse(i,1)] = optScoreReg(Xin,Yin,Ain,m,alpha,lamPath(i),lamFold,Wi,tol,maxit,optOption);
        end

        Wtc{i,1} = Wt{i,1};
        ssNum(i,1) = sum( A'*sum(Wtc{i,1}.^2,2)~=0);
        if(ssNum(i,1) == 0) 
            ssNum(i,1) = nan; 
        end
    end
   
    if(dispPlot)
        figure 
        plot(lamPath,mse,'-*')
        ylabel('MSE');
        xlabel('\lambda')
        drawnow
    end

        
    % Project the results to a feasible region
    for i = 1:max(size(L))
        temp = ssNum-i; % Difference on sensor sparsity
        dt = sort(unique(temp(temp>=0))); dt = dt(~isnan(dt));
        if(isempty(dt))
            warning(strcat("Obtained projection is too sparse to select", num2str(i)," sensors. Maybe consider reduce maximum lambda!"));
            continue;
        end
        t1 = 1;
        t2 = 1;
        while true
            mMse = sort( mse(temp==dt(t1),1), 'ascend' );
            idx = mse==mMse(t2);
            [ss{i,1},fs{i,1},W{i,1},tooSparse] = feasibleSelection(Wtc{idx,1},A,L(i),q(i));
            if(tooSparse && t1<=length(dt))
                if(t2 <length(mMse))
                    t2 = t2 + 1;
                    idx = mse==mMse(t2);
                else
                    t1 = t1 +1;
                    t2 = 1;
                end
            else
                break;
            end
        end
        if(~isempty(G))
             ssN{i,1} = G(ss{i,1});
        end
    end
    
    
        
end
