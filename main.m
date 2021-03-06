%% Calculation of MDF for given dx
ebsd_grid = ebsd.gridify;

%% Sampling algorithm

% define set of correlation vectors in pixels

dX = -50:50;
dY = -50:50;

% Create object to store the spatially correlated mdf
resolved_mdf = cell(length(dX)*length(dY),3);

nsamples = round(0.25*length(ebsd_grid));

[max_x,max_y] = size(ebsd_grid);

element = 0;
error_counter = 0;
    
   
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
                xdx = x + dX(dx) + max_x;
                ydy = y + dY(dy) + max_y;
            elseif x + dX(dx) < 1
                xdx = x + dX(dx) + max_x;
                ydy = y + dY(dy);
            elseif y + dy < 1
                xdx = x + dX(dx);
                ydy = y + dY(dy) + max_y;
            else
                xdx = x + dX(dx);
                ydy = y + dY(dy);
            end
        end
        
        % Case 2
        if dX(dx) < 0 && dY(dy) >= 0
            if x + dX(dx) < 1 && y + dY(dy) > max_y
                xdx = x + dX(dx) + max_x;
                ydy = dY(dy) - (max_y - y);
            elseif x + dX(dx) < 1
                xdx = x + dX(dx) + max_x;
                ydy = y + dY(dy);
            elseif y + dy > max_y
                xdx = x + dX(dx);
                ydy = dY(dy) - (max_y - y);
            else
                xdx = x + dX(dx);
                ydy = y + dY(dy);
            end   
        end
        
        % Case 3
        if dX(dx) >= 0  && dY(dy) < 0
            if x + dX(dx) > max_x && y + dY(dy) < 1
                xdx = dX(dx) - (max_x - x);
                ydy = y + dY(dy) + max_y;
            elseif x + dX(dx) > max_x
                xdx = dX(dx) - (max_x - x);
                ydy = y + dY(dy);
            elseif y + dy < 1
                xdx = x + dX(dx);
                ydy = y + dY(dy) + max_y;
            else
                xdx = x + dX(dx);
                ydy = y + dY(dy);
            end 
        end
        
        % Case 4
        if dX(dx) >= 0 && dY(dy) >= 0
            if x + dx > max_x && y + dy > max_y
                xdx = dX(dx) - (max_x - x);
                ydy = dY(dy) - (max_y - y);
            elseif x + dx > max_x
                xdx = dX(dx) - (max_x - x);
                ydy = y + dY(dy);
            elseif y + dy > max_y
                xdx = x + dX(dx);
                ydy = dY(dy) - (max_y - y);
            else
                xdx = x + dX(dx);
                ydy = y + dY(dy);
            end   
        end
        
        
        ori_head(i) = ebsd_grid(xdx, ydy).orientations;
        
        % if one of the elements are not indexed, remove them
        catch
            
            error_counter = error_counter + 1;
            ori_tail(i) = 0;
            ori_head(i) = 0;
      
        end
        
        if length(ori_tail) < 50
      
            disp(i)
  
        end
        
        
end
% remove 
condition = ori_tail.phi1 == 0 & ori_tail.Phi == 0 & ori_tail.phi2 == 0;
ori_head(condition) = [];
ori_tail(condition) = [];

mori = inv(ori_tail) .* ori_head;

% mdf = calcDensity(mori,'exact','Fourier');
% mdf = calcDensity(mori);
mdf = calcKernelODF(mori);

resolved_mdf{element, 1} = dX(dx);
resolved_mdf{element, 2} = dY(dy);
resolved_mdf{element, 3} = mdf;
        


        end
end
