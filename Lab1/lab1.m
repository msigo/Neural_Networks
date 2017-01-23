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
%%

clc;
clf;
clear all;

[patterns, targets] = nsepdata(200);
[insize, ndata] = size(patterns);
[outsize, ndata] = size(targets);



hidden=40;
eta = 0.1;
alpha = 0.9;
backprop_epoch = 100;

W = 2*(rand(hidden,3) - 0.5*ones(hidden,3));
V = 2*(rand(1,hidden+1) - 0.5*ones(1,hidden+1));
dw = 0;
dv = 0;


for i = 1:backprop_epoch
        plot_title = sprintf('Epoch = %f', i);
              
        [W,V,dw,dv,out] = backprop(W,V,dw,dv,patterns,targets,ndata,hidden,eta,alpha);
        error(i) = sum(sum(abs(sign(out)- targets)./2));
end
figure(2)
plot (error);
title('Error');

%% The encoder problem

clc;
clf;
clear all;

patterns = eye(8)*2 -1;
targets = patterns;

ndata = 8;

nInputLayers = 8;
hidden = 3;
nOutputLayers = 8;

W = 2*(rand(hidden,nInputLayers + 1) - 0.5*ones(hidden,nInputLayers +1 ));
V = 2*(rand(nOutputLayers,hidden+1) - 0.5*ones(nOutputLayers,hidden+1));
dw = 0;
dv = 0;

backprop_epoch = 400;
eta = 0.4;
alpha = 0.6;

for i = 1:backprop_epoch
        plot_title = sprintf('Epoch = %f', i);
              
        [W,V,dw,dv,out] = backprop(W,V,dw,dv,patterns,targets,ndata,hidden,eta,alpha);
        error(i) = sum(sum(abs(sign(out)- targets)./2));
end
figure(2)
plot (error);
title('Error');



%% Function Approximation


clc;
clf;
clear all;

x = [-5:1:5]';
y = x;
z = exp(-x.*x*0.1)*exp(-y.*y*0.1)' - 0.5;



ndata = size(x,1)*size(y,1);
gridsize = size(x,1);
targets = reshape(z,1,ndata);
[xx,yy] = meshgrid(x,y);
patterns = [reshape(xx,1,ndata); reshape(yy,1,ndata)];

nInputLayers = 2;
hidden = 13;
nOutputLayers = 1;

W = 2*(rand(hidden,nInputLayers + 1) - 0.5*ones(hidden,nInputLayers +1 ));
V = 2*(rand(nOutputLayers,hidden+1) - 0.5*ones(nOutputLayers,hidden+1));
dw = 0;
dv = 0;

backprop_epoch = 200;
eta = 0.1;
alpha = 0.9;

for i = 1:backprop_epoch
        plot_title = sprintf('Epoch = %f', i);
              
        [W,V,dw,dv,out] = backprop(W,V,dw,dv,patterns,targets,ndata,hidden,eta,alpha);
        
        figure(2)
        subplot(1,2,1)
        zz = reshape(out, gridsize, gridsize);
        mesh(x,y,zz);
        axis([-5 5 -5 5 -0.7 0.7]);
        title('Approximated function')
        
        subplot(1,2,2)
        mesh(x,y,z)
        axis([-5 5 -5 5 -0.7 0.7]); 
        title('Target function')
        
        drawnow        
end











