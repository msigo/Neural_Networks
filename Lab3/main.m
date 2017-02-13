%% Topological Ordering of Animal Species

animals

W = rand(100,84);
epochs = 20;
[nAnimals, ~] = size(props);
eta = 0.2;

for i = 1:epochs
    for j = 1:nAnimals
       
        p = props(j,:);
        
        dist = pdist2(p, W, 'euclidean'); 
        [mindist,min_idx] = min(dist);                                            %// minimum distance
        
        % Define neighbours to the index,
        neighbours = 40-i*2;
        
        lowerBound = min_idx - neighbours;
        if lowerBound < 1
            lowerBound = 1;
        end
        
        higherBound = min_idx + neighbours;
        if higherBound > size(W,1)
            higherBound = size(W,1);
        end
      
        
        % Update neighbours and the index
        W(lowerBound:higherBound,:) = W(lowerBound:higherBound,:) + eta*(repmat(p,higherBound-lowerBound+1,1) - W(lowerBound:higherBound,:)); %eta*(props(min_idx-neigbours:min_idx+neigbours,:)-W(min_idx-neigbours:min_idx+neigbours,:));%(pdist2(p, W(min_idx-neigbours:min_idx+neigbours), 'euclidean'))
        
    end
end



for j = 1:nAnimals

    p = props(j,:);

    dist = pdist2(p, W, 'euclidean'); 
    [mindist,min_idx] = min(dist);                                            %// minimum distance
    
    pos(j) = min_idx;
end

[dummy, order] = sort(pos);
snames(order)'

%% SALdlasd

cities

W = rand(10,2);
epochs = 20;
[nCities, ~] = size(props);
eta = 0.2;

for i = 1:epochs
    for j = 1:nCities
       
        p = props(j,:);
        
        dist = pdist2(p, W, 'euclidean'); 
        [mindist,min_idx] = min(dist);                                            %// minimum distance
        
        % Define neighbours to the index,
        neighbours = 2;
    
        bounds = [];
        lowerBound = min_idx - neighbours;
        higherBound = min_idx + neighbours;
        
        if lowerBound < 1
            bounds = [1:higherBound, (size(W,1) + lowerbound):size(W,1)]; 
        if higherBound > size(W,1)
            bounds = [1:higherBound-size(W,1), lowerBound:size(W,1)];
        end
      
        
        % Update neighbours and the index
        W(bounds,:) = W(bounds,:) + eta*(repmat(p,size(bounds,2)+1,1) - W(bounds,:)); %eta*(props(min_idx-neigbours:min_idx+neigbours,:)-W(min_idx-neigbours:min_idx+neigbours,:));%(pdist2(p, W(min_idx-neigbours:min_idx+neigbours), 'euclidean'))
        
    end
end


tour = [W;W(1,:)];
plot(tour(:,1),tour(:,2),'b-*',city(:,1), city(:,2),'+')



