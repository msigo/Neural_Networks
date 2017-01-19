function delta_W = delta(W, patterns, target)
    eta = 0.001;
    patterns = [patterns; ones(1,size(patterns,2))];
    
    delta_W = -eta * (W*patterns - target)*patterns';
end

