%-----------------------------------------------%
% Function: project the projection matrix W to a feasible region under sensor and feature number constraints
% W: d x m projection matrix
% A: feature sensors relationship matrix
% L: feature number constraint
% q: sensors number constraint
%-----------------------------------------------%
function [ss,fs,Wf,tooSparse] = feasibleSelection(W,A,L,q)
    % Initialisation
    d = size(A,1);
    p = size(A,2);
    tooSparse = false;
    fs = zeros(d,1);
    ss = zeros(p,1);
    Wf = zeros(size(W,1),size(W,2));
    
    % Check if we can find feasible solution
    temp = sum(W.^2,2);
    [~, I] = sort(temp,'descend');
    feature = I(1:sum(temp~=0)); % Feature index sorted by the weights in W
    sensor = find( (A'* temp~=0) ~= 0 );
    if( size(sensor,1)<q) %size(feature,1)<L ||
        tooSparse = true; % Too few sensors, need to decrease the parameter
        ss = [];
        fs = [];
        Wf = W;
        return;
    end
    
    % Find the feasible solution
    for t = 1:size(feature,1)
        tempF = fs;
        tempF(feature(t),1) = 1;
        tempS = double((A'* tempF~=0) ~= 0);
        
        if(sum(tempF)<=L && sum(tempS) <= q && all(ismember(find(tempS),sensor)) )
            fs = tempF;
            ss = tempS;
            Wf(feature(t),:) = W(feature(t),:);
        end 
    end
    
    % Check feasibility
    fs = find(fs==1);
    ss = find(ss==1);

    % sensor number equality constraint lifted
    if(size(ss,1) < q)
        tooSparse = true;
    end
   
    

end

