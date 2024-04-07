%-----------------------------------------------%
% Function: Express the sensor feature relationship in a mapping matrix form
%-----------------------------------------------%
function A = relationsMatrix(S,F,source)
    d = max(size(F)); % dimension
    p = max(size(S));

    A = zeros(d,p);
    for i = 1:d
        for j = 1:max(source.RequireSensor)
            sensor = string(table2array(source(i,2+j)));
            idx = find(strcmp(S,sensor));
            A(i,idx) = 1;
        end
    end
end

