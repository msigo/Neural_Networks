function delta_W = delta(W, patterns, target, eta)
    if nargin < 4
        eta = 0.1;
    end
    
    patterns = [patterns; ones(1,size(patterns,2))];
    
    delta_W = eta * (W*patterns - target)*patterns';
end

