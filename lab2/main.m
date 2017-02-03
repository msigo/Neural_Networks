
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


















