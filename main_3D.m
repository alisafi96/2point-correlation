%% Calculation of MDF for given dx
ebsd_grid = ebsd.gridify;

%% Sampling algorithm

% define set of correlation vectors in pixels

dX = -1:1;
dY = -1:1;
dZ = -1:1;

% Create object to store the spatially correlated mdf
resolved_mdf = cell(length(dX)*length(dY)*length(dZ),4);

nsamples = 0.5*length(ebsd_grid);

[max_x,max_y] = size(ebsd_grid);

element = 0;

for dz = 1:length(dZ)
    
    for dy = 1:length(dY)
        
        for dx = 1:length(dX)
            element = element + 1;
% Create random set of samples
sampleset = sort(randperm(length(ebsd_grid),nsamples));
X = floor(sampleset ./ max_y)+1;
Y = rem(sampleset,max_y)+1;

%i = (x-1)*max_y + y;
% initialization to get a dummy ori object
ori_tail = ebsd_grid(1:length(sampleset)).orientations;
ori_head = ebsd_grid(1:length(sampleset)).orientations;

for i = 1:length(X)
       
        % current coordinates
        x = X(i);
        y = Y(i);
    
        try
        
        ori_tail(i) = ebsd_grid(x,y).orientations;
        

        %Case 1
        if dX(dx) < 0 && dY(dy) < 0
            if x + dX(dx) < 1 && y + dY(dy) < 1
                ori_head(i) = ebsd_grid(x + dX(dx) + max_x,y + dY(dy) + max_y).orientations;
            elseif x + dX(dx) < 1
                ori_head(i) = ebsd_grid(x + dX(dx) + max_x,y + dY(dy)).orientations;
            elseif y + dy < 1
                ori_head(i) = ebsd_grid(x + dX(dx),y + dY(dy) + max_y).orientations;
            else
                ori_head(i) = ebsd_grid(x + dX(dx),y + dY(dy)).orientations;
            end
        end
        
        
        % Case 2
        if dX(dx) < 0 && dY(dy) >= 0
            if x + dX(dx) < 1 && y + dY(dy) >= max_Y
                ori_head(i) = ebsd_grid(x + dX(dx) + max_x,dY(dy) - (max_y - y)).orientations;
            elseif x + dX(dx) < 1
                ori_head(i) = ebsd_grid(x + dX(dx) + max_x,y + dY(dy)).orientations;
            elseif y + dy >= max_Y
                ori_head(i) = ebsd_grid(x + dX(dx),dY(dy) - (max_y - y)).orientations;
            else
                ori_head(i) = ebsd_grid(x + dX(dx),y + dY(dy)).orientations;
            end
            
        end
        
        % Case 3
        if dX(dx) >= 0  && dY(dy) < 0
            if x + dX(dx) > max_x && y + dY(dy) < 1
                ori_head(i) = ebsd_grid(dX(dx) - (max_x - x),y + dY(dy) + max_y).orientations;
            elseif x + dX(dx) > max_x
                ori_head(i) = ebsd_grid(dX(dx) - (max_x - x),y + dY(dy)).orientations;
            elseif y + dy < 1
                ori_head(i) = ebsd_grid(x + dX(dx),dY(dy) - (max_y - y)).orientations;
            else
                ori_head(i) = ebsd_grid(x + dX(dx),y + dY(dy)).orientations;
            end 
        end
        
        % Case 4
        if dX(dx) >= 0 && dY(dy) >= 0
            if x + dx > max_x && y + dy > max_y
                ori_head(i) = ebsd_grid(dx - (max_x - x),dy - (max_y - y)).orientations;
            elseif x + dx > max_x 
                ori_head(i) = ebsd_grid(dx - (max_x - x),y + dy).orientations;
            elseif y + dy > max_y
                ori_head(i) = ebsd_grid(x + dx,dy - (max_y - y)).orientations;
            else
                ori_head(i) = ebsd_grid(x + dx,y + dy).orientations;
            end   
        end



        % if one of the elements are not indexed, remove them
        catch
            ori_tail(i) = 0;
            ori_head(i) = 0;
        end
        
       
end

% remove 
condition = ori_tail.phi1 == 0 & ori_tail.Phi == 0 & ori_tail.phi2 == 0;
ori_head(condition) = [];
ori_tail(condition) = [];

mori = inv(ori_tail) .* ori_head;

% mdf = calcDensity(mori,'exact','Fourier');
mdf = calcDensity(mori);

resolved_mdf{element, 1} = dX(dx);
resolved_mdf{element, 2} = dY(dy);
resolved_mdf{element, 3} = dZ(dz);
resolved_mdf{element, 4} = mdf;
        end
    end
end
