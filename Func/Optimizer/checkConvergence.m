%-----------------------------------------------%
% Function: check convergence of a sequence
%-----------------------------------------------%
function flag = checkConvergence(x,tol)
    n = size(x,1);
    if(n==1)
        flag = false;
        return
    else
        e = abs(diff(x));
        N = 2;
        if(size(e,1)<=N)
            flag = false; %Not converged
            return
        else
            if(max(e(end-N:end)) < tol )
                flag = true; % Converged
            else
                flag = false; %Not converged
            end
        end
    end
   
    
end

