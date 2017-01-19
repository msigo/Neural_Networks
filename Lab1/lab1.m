clc;
clf;
clear all;

[patterns, targets] = sepdata();
%[patterns, targets] = nsepdata();


[insize, ndata] = size(patterns);
[outsize, ndata] = size(targets);

permute = randperm(200);
patterns = patterns(:, permute);
targets = targets(:, permute);
epoch = 40;
W = zeros(1,3);

for i = 1:epoch
    W  = W + delta(W,patterns,targets);
    p = W(1,1:2);
    k = -W(1, insize+1) / (p*p');
    l = sqrt(p*p');
    plot (patterns(1, find(targets>0)), ...
        patterns(2, find(targets>0)), '*', ...
        patterns(1, find(targets<0)), ...
        patterns(2, find(targets<0)), '+', ...
        [p(1), p(1)]*k + [-p(2), p(2)]/l, ...
        [p(2), p(2)]*k + [p(1), -p(1)]/l, '-');
    drawnow;
    axis([-2 2 -2 2],'square')
end



%%
plot (patterns(1, find(targets>0)), ...
    patterns(2, find(targets>0)), '*', ...
    patterns(1, find(targets<0)), ...
    patterns(2, find(targets<0)), '+');

%%