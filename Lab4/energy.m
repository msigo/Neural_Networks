clear all; close all;clc;
%% Energy Calculations
x1 = vm([0 0 1 0 1 0 0 1]);
x2 = vm([0 0 0 0 0 1 0 0]);
x3 = vm([0 1 1 0 0 1 0 1]);

x_start = [x1; x2; x3];

W = x_start'*x_start;

x1d = vm([1 0 1 0 1 0 0 1]);
x2d = vm([1 1 0 0 0 1 0 0]);
x3d = vm([1 1 1 0 1 1 0 1]);

%change more from attractors -> larger value of the energy. Energy get
%smaller for each iteration
x_input = [x1d; x2d; x3d];

x_old = x_input;
while 1
    x_slask = sgn(x_old*W);
    E = -sum(sum(W.*(x_old'*x_old)))

    if x_slask == x_old
        disp('converged');
        break
    end
    
    x_old = x_slask;
    
end

x_update = x_old;


%% Random weight matrix
W = 1.*randn(8,8);

% Make the weight matrix symmetric. (Which guaranties convergence)
W= 0.5*(W'*W);

x_input = [x1d; x2d; x3d];

x_old = x_input;
E = [];
while 1
    x_slask = sgn(x_old*W);
    E = [E -sum(sum(W.*(x_old'*x_old)))];

    if x_slask == x_old
        disp('converged');
        break
    end
    
    x_old = x_slask;
    
end
plot(E)


