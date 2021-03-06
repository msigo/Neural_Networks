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
[nCities, ~] = size(city);
eta = 1;

neighbours = 2;
for i = 1:epochs
    if i == epochs - 15
        neighbours = 1;
    end
    if i == epochs - 5
        neighbours = 0;
    end
    for j = 1:nCities
       
        p = city(j,:);
        
        dist = pdist2(p, W, 'euclidean'); 
        [mindist,min_idx] = min(dist);                                            %// minimum distance
        
        % Define neighbours to the index,
       
        
        lowerBound = min_idx - neighbours;
        higherBound = min_idx + neighbours;
        bounds = [lowerBound:higherBound];
        
        if lowerBound < 1
            bounds = [1:higherBound, (size(W,1) + lowerBound):size(W,1)];
        end
        if higherBound > size(W,1)
            bounds = [1:higherBound-size(W,1), lowerBound:size(W,1)]
        end
    
        % Update neighbours and the index
        W(bounds,:) = W(bounds,:) + eta*(repmat(p,size(bounds,2),1) - W(bounds,:)); %eta*(props(min_idx-neigbours:min_idx+neigbours,:)-W(min_idx-neigbours:min_idx+neigbours,:));%(pdist2(p, W(min_idx-neigbours:min_idx+neigbours), 'euclidean'))
        
    end
end


tour = [W;W(1,:)];
plot(tour(:,1),tour(:,2),'b-*',city(:,1), city(:,2),'r+')

%% Data Clustering: Votes of MPs

votes
x = meshgrid([1:10]);
xpos = reshape(x,1,100);

W = rand(100,31);
epochs = 20;
[nvoters, ~] = size(vote);
eta = 0.2;

for i = 1:epochs
    for j = 1:nvoters
       
        p = vote(j,:);
        
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
        
        xneighbours = [lowerBound:higherBound];
        
        yneighbours = [];
        k = 1
        while min_idx - k*9 > 1 & k < neighbours
            yneighbours = [yneighbours, min_idx - k*9];
            k = k+1;
        end
        
        k = 1
        while min_idx + k*9 < size(W,1) & k < neighbours
            yneighbours = [yneighbours, min_idx + k*9];
            k = k+1;
        end
        
        
        
        bounds = [xneighbours, yneighbours];
      
        
        % Update neighbours and the index
        W(bounds,:) = W(bounds,:) + eta*(repmat(p,size(bounds,2),1) - W(bounds,:));
        
    end
end




for j = 1:nvoters

    p = vote(j,:);

    dist = pdist2(p, W, 'euclidean'); 
    [mindist,min_idx] = min(dist);                                            %// minimum distance
    
    pos(j) = min_idx;
end


%% 
mpparty
a = ones(1,100)*350;
a(pos) = 1:349;

p = [mppartys;0];
image(p(reshape(a,10,10))+1);
descr = {'no party: Black';
    'M: Blue';
    'FP: LightBlue ';
    'S: Pink';
    'V: Red';
    'MP: Green';
    'KD: White';
    'C: Yellow'};
text(-0.5,3,descr)


%%
mpsex
a = ones(1,100)*350;
a(pos) = 1:349;

p = [mpsexs;0];
image(p(reshape(a,10,10))+1);

%%
mpdistrict
a = ones(1,100)*350;
a(pos) = 1:349;

p = [mpdistricts;0];
image(p(reshape(a,10,10))+1);







