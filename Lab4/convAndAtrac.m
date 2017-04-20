clear all, close all, clc

%% Translate the calculation of the weight matrix and the update
% rule into Matlab expressions
x1 = vm([0 0 1 0 1 0 0 1]);
x2 = vm([0 0 0 0 0 1 0 0]);
x3 = vm([0 1 1 0 0 1 0 1]);

x = [x1; x2; x3];

W = x'*x;

x_update = sgn(x*W);

%%

x1d = vm([1 0 1 0 1 0 0 1]);
x2d = vm([1 1 0 0 0 1 0 0]);
x3d = vm([1 1 1 0 1 1 0 1]);

x_input = [x1d; x2d; x3d];
x_old = x_input;

while 1
    x_slask = sgn(x_old*W);
   
    if x_slask == x_old
        disp('converged');
        break
    end
    
    x_old = x_slask;
    
end

x_update = x_old;


%% Search for additonal fix-points 
n = 8;
P = dec2bin(0:(2^n)-1)-'0'
x_old=P;

 while 1
        x_slask = sgn(x_old*W);

        if x_slask == x_old
            disp('converged');
            break
        end

        x_old = x_slask;
 end
 
 x_converged = x_old;
 
 fix_points = sum((sum(P-x_converged,2) == 0));

%% More distortion

x1d = vm([1 1 1 1 1 1 1 1]);
x2d = vm([1 1 1 1 0 0 0 0]);
x3d = vm([1 0 0 1 0 1 0 1]);

x_input = [x1d; x2d; x3d];
x_old = x_input;

while 1
    x_slask = sgn(x_old*W);
   
    if x_slask == x_old
        disp('converged');
        break
    end
    
    x_old = x_slask;
    
end

x_update = x_old;














