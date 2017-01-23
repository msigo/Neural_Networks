function [out] = forwardPass( patterns, w, v, ndata)
    hin = w * [patterns ; ones(1,ndata)];
    hout = [2 ./ (1+exp(-hin)) - 1 ; ones(1,ndata)];
    
    oin = v * hout;
    out = 2 ./ (1+exp(-oin)) - 1;
end

