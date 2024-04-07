%-----------------------------------------------%
% Function: soft thresholding
%-----------------------------------------------%
function out = softTh(x,th)
    I1 = x>th;
    I2 = x<-th;
    I3 = (x<=th & x>=-th);
    xt = x;
    
    xt(I1) = xt(I1) - th;
    xt(I2) = xt(I2) + th;
    xt(I3) = 0;
    
    out = xt;
end
