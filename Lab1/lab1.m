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


eta_array = [0.001, 0.1, 1];
epoch_array = [3, 6, 20];

figure(1)
for i = 1:epoch_array(3)
    
        if i < epoch_array(1)
           subplot_i = 1;

        else if i < epoch_array(2)
           subplot_i = 2;
        
        else if i < epoch_array(3)
           subplot_i = 3;
            end
            end
        end
    
        plot_title = sprintf('Epoch = %f', i);
              
          
        subplot(1,3,subplot_i)
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
        title(plot_title);
        axis([-2 2 -2 2],'square')
        drawnow;
end


figure(2)
for j = 1:size(eta_array,2)
    
        W = zeros(1,3);
        
        
        for i = 1:epoch_array(3)
            plot_title = sprintf('Epoch = %f, eta = %f', i, eta_array(j));
            
            subplot(1,3,j)
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
            title(plot_title);
            axis([-2 2 -2 2],'square')
            drawnow;
        end
end




%%
plot (patterns(1, find(targets>0)), ...
    patterns(2, find(targets>0)), '*', ...
    patterns(1, find(targets<0)), ...
    patterns(2, find(targets<0)), '+');
