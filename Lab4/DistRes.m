pict;
n=5;
p1dist = flipCustom(p1,n);
subplot(2,2,1)
vis(p1);
title('Original')
subplot(2,2,2)
vis(p1dist);
title(sprintf('%d pixels flipped',n));
subplot(2,2,3)
vis(p1dist-p1);
title('Indication of which pixels that were flipped')

%% How many percent of all pixels can be flipped? 
w = p1'*p1 + p2'*p2 + p3'*p3;

curr_pic  = p2;
for i = 1:9:100
    figure;
    curr_pic_dist = flipCustom(curr_pic,i*(1024/100));

    subplot(2,2,1)
    vis(curr_pic_dist);
    title('Flipped picture')
    
    subplot(2,2,3)
    vis(curr_pic_dist-curr_pic);
    
    x_old = curr_pic_dist;
   	for t = 1:30
    
        x_slask = sgn(x_old*w);

        if x_slask == x_old
            disp('converged');
            break
        end
        
        x_old = x_slask;    
    end
    subplot(2,2,2)
    vis(x_slask);
    title(sprintf('Restored for %d \%',i));
end


%% Capacity
close all;
n=100;
R = sign(randn(n,1024));
w = 0;

for i = 1:n
    w = w + R(i,:)'*R(i,:);
end

curr_pic  = R(1,:);
for i = 1:9:101
    figure;
    curr_pic_dist = flipCustom(curr_pic,(i-1)*(1024/100));

    subplot(2,2,1)
    vis(curr_pic_dist);
    title('Flipped Picture')
    
    subplot(2,2,3)
    title('Original Picture')

    vis(curr_pic);
    
    x_old = curr_pic_dist;
   	for t = 1:30
    
        x_slask = sgn(x_old*w);

        if x_slask == x_old
            disp('converged');
            break
        end
        
        x_old = x_slask;    
    end
    subplot(2,2,2)
    vis(x_slask);
    title(sprintf('Restored for %d \%',i));
    subplot(2,2,4)
    vis(sign(x_slask-curr_pic));

end


%% Capacity part 2

R = sign(0.5+randn(300,100));
w=0;
res = zeros(1,300);
for i = 1:300
    w = w + R(i,:)'*R(i,:);
    %w = w-diag(diag(w));

    iter  = 0;
    for t = 1:i-1
        x_slask = sgn(R(t,:)*w);

        if x_slask == R(t,:)
            iter = iter +1;
        end
        
    end 
    res(i) = iter;
end
plot(res)
hold on
plot(1:300)

%% Capacity Part 3


R = sign(0.5+randn(300,100));
w=0;
flipped = 5;
res = zeros(1,300);

for i = 1:300
    w = w + R(i,:)'*R(i,:);
    %w = w-diag(diag(w));
    iter  = 0;
    for t = 1:i-1
        R_flip = flipCustom(R(t,:),flipped);
        x_slask = sgn(R_flip*w);

        if x_slask == R(t,:)
            iter = iter +1;
        end
        
    end 
    res(i) = iter;
end
plot(res)
hold on
plot(1:300)
    
%% Sparse Patterns
N = 100;
P = 300;
rho = 0.1;
R = zeros(P*(1-rho),N);
R = [R;ones(P*rho,N)];
[p,q] = size(R);
ix = randperm(N*P);
R  = R(ix); 
R = reshape(R,p,q);

w=0;
flipped = 5;
res = zeros(1,P);
rho = 1/(N*P) * sum(sum(R))
R_rho = R - rho;
real_res = [];
theta_vec = 0:10:100;
for theta = theta_vec
    w = 0;
    res = 0;
    for i = 1:P
        w = w + R_rho(i,:)'*R_rho(i,:);
        iter  = 0;
        for t = 1:i-1
            x_slask = 0.5+0.5*sgn(R(t,:)*w- theta);

            if x_slask == R(t,:)
                iter = iter +1;
                
            end

        end 
        res(i) = iter;
    end
    real_res = [real_res max(res)];
end
plot(theta_vec,real_res)







