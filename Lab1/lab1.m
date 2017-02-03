clc;
clf;
clear all;

%[patterns, targets] = sepdata();
[patterns, targets] = nsepdata(200);


[insize, ndata] = size(patterns);
[outsize, ndata] = size(targets);

permute = randperm(200);
patterns = patterns(:, permute);
targets = targets(:, permute);
epoch = 40;
W = zeros(1,3);


eta_array = [0.00001, 0.1, 10];
epoch_array = [3, 6, 40];

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
    
        plot_title = sprintf('Epoch = %f , eta = %f', i, 0.1);
              
        %subplot(1,3,subplot_i)
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



hiddens=[2,4,8,16];
eta = 0.1;
alpha = 0.9;
backprop_epoch = 100;



hold on

N = 50;
error = zeros(size(hiddens,2),backprop_epoch);
for j = 1:size(hiddens,2)
    for k = 1:N
        W = 2*(rand(hiddens(j),3) - 0.5*ones(hiddens(j),3));
        V = 2*(rand(1,hiddens(j)+1) - 0.5*ones(1,hiddens(j)+1));
        dw = 0;
        dv = 0;
        for i = 1:backprop_epoch
                plot_title = sprintf('Epoch = %f', i);

                [W,V,dw,dv,out] = backprop(W,V,dw,dv,patterns,targets,ndata,hiddens(j),eta,alpha);
                error(j,i) = error(j,i) +sum(sum(abs(sign(out)- targets)./2));
        end


    end
    
    plot (error(j,:)/N, 'DisplayName', sprintf('%.f hidden neurons',hiddens(j)));
    title(sprintf('Average (N=%.f)Number of errors',N))
    
end
legend('show')
hold off

%% The encoder problem

%clc;
%clf;
%clear all;


patterns = eye(8)*2 -1;
targets = patterns;
ndata = 8;
nInputLayers = 8;
hidden = 3;
nOutputLayers = 8;
N = 1000;

W_result = zeros(hidden, nInputLayers+1);

for j = 1:N
    W = 2*(rand(hidden,nInputLayers + 1) - 0.5*ones(hidden,nInputLayers +1 ));
    V = 2*(rand(nOutputLayers,hidden+1) - 0.5*ones(nOutputLayers,hidden+1));
    dw = 0;
    dv = 0;
    alpha = 0.6;
    backprop_epoch = 400;
    eta = 0.4;

    
    for i = 1:backprop_epoch
            plot_title = sprintf('Epoch = %f', i);
            [W,V,dw,dv,out] = backprop(W,V,dw,dv,patterns,targets,ndata,hidden,eta,alpha);
            (sign(W) + 1)/2;
            error(i) = sum(sum(abs(sign(out)- targets)./2));
    end
    figure(1)
    plot (error);
    title('Error');
    %W_result = (W_result + W)/j;
end

(sign(W) + 1)/2
%(sign(V) + 1)/2

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



backprop_epoch = 150;
eta = 0.1;
alpha = 0.9;



figure(1)
subplot(2,3,1)
mesh(x,y,z)
axis([-5 5 -5 5 -0.7 0.7]); 
title('Target function')

hiddens = [2,4,6,9,13];


for j = 1:5
    hidden = hiddens(j)
    W = 2*(rand(hidden,nInputLayers + 1) - 0.5*ones(hidden,nInputLayers +1 ));
    V = 2*(rand(nOutputLayers,hidden+1) - 0.5*ones(nOutputLayers,hidden+1));
    dw = 0;
    dv = 0;
    
  
    
    subplot(2,3,j+1)
    for i = 1:backprop_epoch
            plot_title = sprintf('Epoch = %f', i);

            [W,V,dw,dv,out] = backprop(W,V,dw,dv,patterns,targets,ndata,hidden,eta,alpha);


            %subplot(1,2,1)
            zz = reshape(out, gridsize, gridsize);
            mesh(x,y,zz);
            axis([-5 5 -5 5 -0.7 0.7]);
            title(sprintf('%.0f nodes, %.0f epochs',hidden, backprop_epoch))

            %subplot(1,2,2)


            drawnow        
    end
end

%% Generalization 


clc;
clf;
clear all;

n = 10;
hidden= 30;
eta = 0.01;
alpha = 0.9;
backprop_epoch = 200;


[patterns, targets] = nsepdata(200);
permute = randperm(200);
patterns = patterns(:, permute);
targets = targets(:, permute);



trainingPatterns = patterns(:, 1:n);
trainingTargets = targets(:,1:n);
validationPatterns = patterns(:,n:end); 
validationTargets = targets(:,n:end);

[insize, ndataTrain] = size(trainingPatterns);
[outsize, ndataTrain] = size(trainingTargets);

[insize, ndataValidate] = size(validationPatterns);
[outsize, ndataValidate] = size(validationTargets);

hiddens = [1,3, 10, 30, 50,200];

%rng(1);

%set(0,'defaultaxescolororder',[0 0 0; 0.5 0.5 0.5]) %black and gray
%set(0,'defaultaxeslinestyleorder',{'-+','-o','-*','-.','-x','s','d','^','v','>','<','p','h'})
figure(1)


for j = 1:size(hiddens,2);
    hold on
    hidden = hiddens(j)
    W = 2*(rand(hidden,3) - 0.5*ones(hidden,3));
    V = 2*(rand(1,hidden+1) - 0.5*ones(1,hidden+1));
    dw = 0;
    dv = 0;


    for i = 1:backprop_epoch
            plot_title = sprintf('Epoch = %f', i);

            [W,V,dw,dv,out] = backprop(W,V,dw,dv,trainingPatterns,trainingTargets,ndataTrain,hidden,eta,alpha);

            out = forwardPass(validationPatterns, W,V,ndataValidate);
            error(i) = sum(sum(abs(sign(out)- validationTargets)./2));
    end
    
    subplot(1,2,1)
    plot (error, 'DisplayName', sprintf('%.f hidden neurons',hiddens(j)));
    title(sprintf('Number of errors with %.f training sets', n))
   % axis([0 30 0 150],'square')
    
    %subplot(1,2,2)
    %plot (error, 'DisplayName', sprintf('%.f hidden neurons',hiddens(j)));
    %title(sprintf('Number of errors with %.f training sets', n))
    axis([100 200 0 100],'square') 
end
hold off
legend('show')





hiddens = [1,3, 10, 30, 50,200];
n = 25
trainingPatterns = patterns(:, 1:n);
trainingTargets = targets(:,1:n);
validationPatterns = patterns(:,n:end); 
validationTargets = targets(:,n:end);

[insize, ndataTrain] = size(trainingPatterns);
[outsize, ndataTrain] = size(trainingTargets);

[insize, ndataValidate] = size(validationPatterns);
[outsize, ndataValidate] = size(validationTargets);


for j = 1:size(hiddens,2);
    hold on
    hidden = hiddens(j)
    W = 2*(rand(hidden,3) - 0.5*ones(hidden,3));
    V = 2*(rand(1,hidden+1) - 0.5*ones(1,hidden+1));
    dw = 0;
    dv = 0;


    for i = 1:backprop_epoch
            plot_title = sprintf('Epoch = %f', i);

            [W,V,dw,dv,out] = backprop(W,V,dw,dv,trainingPatterns,trainingTargets,ndataTrain,hidden,eta,alpha);

            out = forwardPass(validationPatterns, W,V,ndataValidate);
            error(i) = sum(sum(abs(sign(out)- validationTargets)./2));
    end
    
        
    subplot(1,2,2)
    plot (error, 'DisplayName', sprintf('%.f hidden neurons',hiddens(j)));
    title(sprintf('Number of errors with %.f training sets', n))
    
    %axis([0 30 0 150],'square')
    %subplot(1,2,2)
    %plot (error, 'DisplayName', sprintf('%.f hidden neurons',hiddens(j)));
    %title(sprintf('Number of errors with %.f training sets', n))
    axis([100 200 0 100],'square')

    
    
end
legend('show')
hold off


%%

clc;
clf;
clear all;

n = 10;
hidden= 30;
eta = 0.2;
alpha = 0.9;
backprop_epoch = 200;


[patterns, targets] = nsepdata(200);
permute = randperm(200);
patterns = patterns(:, permute);
targets = targets(:, permute);



trainingPatterns = patterns(:, 1:n);
trainingTargets = targets(:,1:n);
validationPatterns = patterns(:,n:end); 
validationTargets = targets(:,n:end);

[insize, ndataTrain] = size(trainingPatterns);
[outsize, ndataTrain] = size(trainingTargets);

[insize, ndataValidate] = size(validationPatterns);
[outsize, ndataValidate] = size(validationTargets);


hiddens = 1:2:100;

%rng(1);

%set(0,'defaultaxescolororder',[0 0 0; 0.5 0.5 0.5]) %black and gray
%set(0,'defaultaxeslinestyleorder',{'-+','-o','-*','-.','-x','s','d','^','v','>','<','p','h'})
figure(1)


for j = 1:size(hiddens,2);
    hold on
    hidden = hiddens(j);
    W = 2*(rand(hidden,3) - 0.5*ones(hidden,3));
    V = 2*(rand(1,hidden+1) - 0.5*ones(1,hidden+1));
    dw = 0;
    dv = 0;


    for i = 1:backprop_epoch
            plot_title = sprintf('Epoch = %f', i);

            [W,V,dw,dv,out] = backprop(W,V,dw,dv,trainingPatterns,trainingTargets,ndataTrain,hidden,eta,alpha);

            out = forwardPass(validationPatterns, W,V,ndataValidate);
            %error(i) = sum(sum(abs(sign(out)- validationTargets)./2));
    end
    error(j) = sum(sum(abs(sign(out)- validationTargets)./2));
    
end
size(error)
subplot(1,2,1)
plot (error)
title(sprintf('Number of errors with %.f training sets', n))
%legend('show')
hold off




n = 25
trainingPatterns = patterns(:, 1:n);
trainingTargets = targets(:,1:n);
validationPatterns = patterns(:,n:end); 
validationTargets = targets(:,n:end);

[insize, ndataTrain] = size(trainingPatterns);
[outsize, ndataTrain] = size(trainingTargets);

[insize, ndataValidate] = size(validationPatterns);
[outsize, ndataValidate] = size(validationTargets);


for j = 1:size(hiddens,2);
    hold on
    hidden = hiddens(j);
    W = 2*(rand(hidden,3) - 0.5*ones(hidden,3));
    V = 2*(rand(1,hidden+1) - 0.5*ones(1,hidden+1));
    dw = 0;
    dv = 0;


    for i = 1:backprop_epoch
            plot_title = sprintf('Epoch = %f', i);

            [W,V,dw,dv,out] = backprop(W,V,dw,dv,trainingPatterns,trainingTargets,ndataTrain,hidden,eta,alpha);

            out = forwardPass(validationPatterns, W,V,ndataValidate);    
    end
    error(j) = sum(sum(abs(sign(out)- validationTargets)./2));
    
end
subplot(1,2,2)
plot (error)
title(sprintf('Number of errors with %.f training sets', n))
%legend('show')
hold off



%%

clc;
clf;
clear all;

n = 10;
hidden= 30;
eta = 0.05;
alpha = 0.9;
backprop_epoch = 500;


[patterns, targets] = nsepdata(200);
permute = randperm(200);
patterns = patterns(:, permute);
targets = targets(:, permute);



trainingPatterns = patterns(:, 1:n);
trainingTargets = targets(:,1:n);
validationPatterns = patterns(:,n:end); 
validationTargets = targets(:,n:end);

[insize, ndataTrain] = size(trainingPatterns);
[outsize, ndataTrain] = size(trainingTargets);

[insize, ndataValidate] = size(validationPatterns);
[outsize, ndataValidate] = size(validationTargets);

hiddens = [1,3, 10, 30, 50,200];

%rng(1);

N = 50;

%set(0,'defaultaxescolororder',[0 0 0; 0.5 0.5 0.5]) %black and gray
%set(0,'defaultaxeslinestyleorder',{'-+','-o','-*','-.','-x','s','d','^','v','>','<','p','h'})


figure(1)
error = zeros(size(hiddens,2),backprop_epoch);
for j = 1:size(hiddens,2);

    for k = 1:N
        hold on
        hidden = hiddens(j)
        W = 2*(rand(hidden,3) - 0.5*ones(hidden,3));
        V = 2*(rand(1,hidden+1) - 0.5*ones(1,hidden+1));
        dw = 0;
        dv = 0;


        for i = 1:backprop_epoch
                plot_title = sprintf('Epoch = %f', i);

                [W,V,dw,dv,out] = backprop(W,V,dw,dv,trainingPatterns,trainingTargets,ndataTrain,hidden,eta,alpha);

                out = forwardPass(validationPatterns, W,V,ndataValidate);
                error(j,i) = error(j,i) + sum(sum(abs(sign(out)- validationTargets)./2));
        end
    

    end
    
    %subplot(1,2,1)
    plot (error(j,:)/N, 'DisplayName', sprintf('%.f hidden neurons',hiddens(j)));
    title(sprintf('Average error (N= %.f) with %.f training sets', N,n))
    %axis([0 50 0 150],'square')
    %axis([100 200 0 100],'square') 
    
end
hold off
legend('show')





hiddens = [1,3, 10, 30, 50,200];
n = 25
trainingPatterns = patterns(:, 1:n);
trainingTargets = targets(:,1:n);
validationPatterns = patterns(:,n:end); 
validationTargets = targets(:,n:end);

[insize, ndataTrain] = size(trainingPatterns);
[outsize, ndataTrain] = size(trainingTargets);

[insize, ndataValidate] = size(validationPatterns);
[outsize, ndataValidate] = size(validationTargets);

figure(2)

backprop_epoch = 500;
error = zeros(size(hiddens,2),backprop_epoch);
for j = 1:size(hiddens,2);
    for k = 1:N
        hold on
        hidden = hiddens(j)
        W = 2*(rand(hidden,3) - 0.5*ones(hidden,3));
        V = 2*(rand(1,hidden+1) - 0.5*ones(1,hidden+1));
        dw = 0;
        dv = 0;

        for i = 1:backprop_epoch
                plot_title = sprintf('Epoch = %f', i);

                [W,V,dw,dv,out] = backprop(W,V,dw,dv,trainingPatterns,trainingTargets,ndataTrain,hidden,eta,alpha);

                out = forwardPass(validationPatterns, W,V,ndataValidate);
                error(j,i) = error(j,i) + sum(sum(abs(sign(out)- validationTargets)./2));
        end
    end

        
    %subplot(1,2,2)
    plot (error(j,:)/N, 'DisplayName', sprintf('%.f hidden neurons',hiddens(j)));
    title(sprintf('Average error (N= %.f) with %.f training sets', N,n))
    %axis([0 50 0 150],'square')
    
    %subplot(1,2,2)
    %plot (error, 'DisplayName', sprintf('%.f hidden neurons',hiddens(j)));
    %title(sprintf('Number of errors with %.f training sets', n))
    %axis([100 200 0 100],'square')

    
    
end
legend('show')
hold off





