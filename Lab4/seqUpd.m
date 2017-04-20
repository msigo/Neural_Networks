clear all, close all, clc

pict;
w = p1'*p1 + p2'*p2 + p3'*p3;

x_rand = sgn(randn(2,1024));  

x_input = x_rand;
x_old = x_input;


vis(x_old);
counter = 0;
while 1
    
    x_slask = sgn(x_old*w);
   
    if x_slask == x_old
        disp('converged');
        break
    end
    
    figure(2)
    vis(x_slask);
    

    
    x_old = x_slask;
    counter = counter +1;
end

x_update = x_old;

figure(3)
vis(x_update);
counter





