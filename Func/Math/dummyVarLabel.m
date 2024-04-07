%-----------------------------------------------%
% Function: Convert label vector to matrix (dummy variable)
%-----------------------------------------------%

function Y = dummyVarLabel(y)
    n = max(size(y,1));
    c = max(size(unique(y)));
    Y = zeros(n,c);
    for i = 1:c
        Y(y==i,i) = 1;
    end
end

