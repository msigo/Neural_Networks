function [patterns, targets] = nsepdata(nData)
    classA(1,:) = [randn(1,nData/4) .* 0.2 - 1.0, ...
    randn(1,nData/4) .* 0.2 + 1.0];
    classA(2,:) = randn(1,nData/2) .* 0.2 + 0.3;
    classB(1,:) = randn(1,nData/2) .* 0.3 + 0.0;
    classB(2,:) = randn(1,nData/2) .* 0.3 - 0.1;

    patterns = [classA, classB];
    targets = [ones(1,nData/2) -ones(1,nData/2)];
end

