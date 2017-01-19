function [patterns, targets] = sepdata( )

classA(1,:) = randn(1,100) .* 0.5 + 1.0;
classA(2,:) = randn(1,100) .* 0.5 + 0.5;
classB(1,:) = randn(1,100) .* 0.5 - 1.0;
classB(2,:) = randn(1,100) .* 0.5 +0;

patterns = [classA, classB];
targets = [ones(1,100) -ones(1,100)];

end

