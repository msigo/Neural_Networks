
%% Batch mode training using least squares, sin(2x)
clear all, close all

% 1)
x = 0:0.1:2*pi;
x = x';

% 2)
f = sin(2*x);

% 3)
units = 7; % 0.1 maximum absolute error.
%units = 25;  % 0.01 maximum absolute error.
%units = 56; % 0.001 maxumum absolute error

%units = 63; n = N;

makerbf
Phi = calcPhi(x,m,var);
% 4)
w = Phi\f;
y = Phi*w;
rbfplot1(x,y,f,units);

%% Bath mode training using least squares, square(2x)

clear all, close all

% 1)
x = 0:0.1:6.3;
x = x';

% 2)
f = square(2*x);

% 3)
%units = 56; % 0.1 maximum absolute error.
%units = 57;  % bigger ?!?!
%units = 63; % 0.06 error
%units = 64; % 0 error
%units = 56; % 0.001 maxumum absolute error

%units = 63; n = N;

makerbf
Phi = calcPhi(x,m,var);
% 4)
w = Phi\f;
y = Phi*w;
rbfplot1(x,y,f,units);


%% On-line training using the delta rule

clear all, close all

fun = 'sin2x';
eta = 0.6;
units = 50;


% 1)
x = 0:0.1:2*pi;
x = x';

% 2)
f = sin(2*x);

makerbf
itermax = itermax*40; 
itersub = 1000;

diter


%% On-line training with other function

clear all, close all

fun = 'gaussian';
eta = 0.6;
units = 9;


% 1)
x = 0:0.1:5;
x = x';

makerbf
itermax = itermax*2; 
itersub = 100;

diter


%% RBF Placement by Self Organization
clear all, close all

plotinit 

data = read('cluster');

units = 5;

vqinit

singlewinner = 1;

%vqiter

%%
clear all, close all

plotinit 

data = read('cluster');

units = 5;

vqinit

singlewinner = 0;

vqiter


%%

clear all, close all

plotinit 

data = read('cluster');

units = 4;

vqinit

singlewinner = 0;

emiterb;


%%
clear all, close all


plotinit 

data = read('cluster');

units = 4
vqinit

singlewinner = 1;

emiterb


%% Function Approximation

clear all, close all





for j = 1:2

    lowpass = j-1;

    [xtrain ytrain]=readxy('ballist',2,2);
    [xtest ytest]=readxy('balltest',2,2);

    if lowpass
        xtrain(:,1) = smooth(xtrain(:,1));
        xtrain(:,2) = smooth(xtrain(:,2));
        ytrain(:,1) = smooth(ytrain(:,1));
        ytrain(:,2) = smooth(ytrain(:,2));
    end

    units=20;
    data=xtrain;
    rng(1)
    vqinit;
    singlewinner=1;
    emiterb

    Phi=calcPhi(xtrain,m,var);

    d1=ytrain(:,1);
    d2=ytrain(:,2);
    dtest1=ytest(:,1);
    dtest2=ytest(:,2);

    w1 = Phi\d1;
    w2 = Phi\d2;

    y1 = Phi*w1;
    y2 = Phi*w2;


    Phitest=calcPhi(xtest,m,var); 
    ytest1=Phitest*w1;
    ytest2=Phitest*w2;
    
    figure('Name',sprintf('Smooting = %0.1f', j-1),'NumberTitle','off')
    %figure(j);
    subplot(2,2,1)
    xyplot(d1,y1,'train1');
    subplot(2,2,2)
    xyplot(d2,y2,'train2');
    subplot(2,2,3)
    xyplot(dtest1,ytest1,'test1');
    subplot(2,2,4)
    xyplot(dtest2,ytest2,'test2');
end

















