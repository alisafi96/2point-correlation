%%=========================================================================
% MTEX (5.3.1) script for grain size, morphology and texture anaylsis of Addalloy
%
%
% Author: Ali Reza Safi
%         Steel Institute (IEHK), RWTH Aachen University
%
% Short overview of functionalities:
%    Section 1  Import variables
%    Section 2  Statistical analysis (grain size, aspect ratios and shape axis orientation)
%%=========================================================================
%% Section 1 Load workspace variables

load('ebsd-and-grains.mat')

%% Section 2.1.1 Grain size analysis (area and number fractions)

% Select the grain object you want to analyze
grains_selected = grains;

% Minimum grain size (euivalent diameter) to be considered (in micron)
min_grain_size = 0.1;

% Number of bins in grain analysis
n_bins = 70;

condition = grains_selected.equivalentRadius > (min_grain_size / 2.0);
grain_area = grains_selected(condition).area;
total_area = sum(grain_area);
grain_size = 2.0 * grains_selected(condition).equivalentRadius;
n_grains = length(grain_size);

[bin_counts,bin_centers] = hist(grain_size,n_bins);
bin_counts = bin_counts';
bin_centers = bin_centers';
bin_width = bin_centers(3) - bin_centers(2);

bin_low = zeros(n_bins,1);
bin_high = zeros(n_bins,1);
for i = 1 : n_bins
    bin_low(i) = bin_centers(i) - bin_width / 2.0;
    bin_high(i) = bin_centers(i) + bin_width / 2.0;
end
max_grain_size = max(grain_size);
min_grain_size = min(grain_size);
bin_low(1) = min_grain_size;
bin_high(n_bins) = max_grain_size;

area_per_bin = zeros(n_bins,1);
area_counts_per_bin = zeros(n_bins,1);
for i = 1 : n_grains
      for j = 1 : n_bins
           if (grain_size(i) >= bin_low(j)) && (grain_size(i) <= bin_high(j))
              area_per_bin(j) = area_per_bin(j) + grain_area(i);
              area_counts_per_bin(j) = area_counts_per_bin(j) + 1;
           end
      end
end

bin_area_fraction = zeros(n_bins,1);
bin_area_fraction = 100.0 * area_per_bin / total_area;
bin_number_fraction = zeros(n_bins,1);
bin_number_fraction = 100.0 * bin_counts / n_grains;

mean_grain_size_area   = sum(bin_centers .* bin_area_fraction / 100.0)
mean_grain_size_number = sum(bin_centers .* bin_number_fraction / 100.0)
effective_grain_size = (mean_grain_size_area + mean_grain_size_number) / 2.0

figure
h = bar(bin_centers,[bin_area_fraction,bin_number_fraction],'histc','FaceColor','flat');
h(1,1).FaceColor = 'red';
h(1,2).FaceColor = 'blue';

l = cell(1,2);
l{1} = 'Area fraction';
l{2} = 'Number fraction';

legend(h,l);
set(gca,'FontSize',15)
% xticks(0:20:(bin_centers(n_bins)+ xticks_resolution));
% ylim([0 30])
% yticks(0:5:30);

xlabel('Grain size (\it{d}\rm) [\mum]')
ylabel('Fraction [%]')

% Fitting a log-normal function to the histogram
logn_fit = lognfit(grain_size);
mu = logn_fit(1)
sigma = logn_fit(2)
x = 0:0.1:round(max_grain_size);
hold on
p = plot(x,100.0 * bin_width * lognpdf(x,mu,sigma))
hold off
set(p,'DisplayName','Lognormal fit','LineStyle','--','LineWidth',2.5)

clear min_grain_size n_bins condition grain_area total_area grain_size n_grains
clear bin_counts bin_centers bin_width bin_low bin_high area_per_bin area_counts_per_bin 
clear bin_area_fraction bin_number_fraction h l logn_fit x p i j grains_selected

%% Section 2.1 Grain size analysis (only area fractions)

% Select the grain object you want to analyze
grains_selected = grains;

% Minimum grain size (euivalent diameter) to be considered (in micron)
min_grain_size = 0.1;

% Number of bins in grain analysis
n_bins = 200;

condition = grains_selected.equivalentRadius > (min_grain_size / 2.0);
grain_area = grains_selected(condition).area;
total_area = sum(grain_area);
grain_size = 2.0 * grains_selected(condition).equivalentRadius;
n_grains = length(grain_size);

[bin_counts,bin_centers] = hist(grain_size,n_bins);
bin_counts = bin_counts';
bin_centers = bin_centers';
bin_width = bin_centers(3) - bin_centers(2);

bin_low = zeros(n_bins,1);
bin_high = zeros(n_bins,1);
for i = 1 : n_bins
    bin_low(i) = bin_centers(i) - bin_width / 2.0;
    bin_high(i) = bin_centers(i) + bin_width / 2.0;
end
max_grain_size = max(grain_size);
min_grain_size = min(grain_size);
bin_low(1) = min_grain_size;
bin_high(n_bins) = max_grain_size;

area_per_bin = zeros(n_bins,1);
area_counts_per_bin = zeros(n_bins,1);
for i = 1 : n_grains
      for j = 1 : n_bins
           if (grain_size(i) >= bin_low(j)) && (grain_size(i) <= bin_high(j))
              area_per_bin(j) = area_per_bin(j) + grain_area(i);
              area_counts_per_bin(j) = area_counts_per_bin(j) + 1;
           end
      end
end

bin_area_fraction = zeros(n_bins,1);
bin_area_fraction = 100.0 * area_per_bin / total_area;

mean_grain_size_area   = sum(bin_centers .* bin_area_fraction / 100.0)

figure
h = bar(bin_centers,[bin_area_fraction],'histc','FaceColor','flat');
h(1).FaceColor = 'red';

bin = zeros(n_bins, 2);
bin(:,1) = bin_centers;
bin(:,2) = bin_area_fraction;

l = 'Area fraction';

legend(h,l);
set(gca,'FontSize',15)

xlabel('Grain size (\it{d}\rm) [\mum]')
ylabel('Fraction [%]')

% Fitting a log-normal function to the histogram



% logn_fit = fitdist(bin,'Lognormal');
% mu = logn_fit(1)
% sigma = logn_fit(2)
% x = 0:0.1:round(max_grain_size);
% af_epdf = pdf(bin_area_fraction,x);
% 
% hold on
% p = plot(x,100.0 * bin_width * af_epdf)
% hold off
% set(p,'DisplayName','Lognormal fit','LineStyle','--','LineWidth',2.5)

clear min_grain_size n_bins condition grain_area total_area grain_size n_grains
clear bin_counts bin_centers bin_width bin_low bin_high area_per_bin area_counts_per_bin 
clear bin_area_fraction bin_number_fraction h l logn_fit x p i j grains_selected

%% Section 2.1.1 Grain size analysis (only number fractions)

% Select the grain object you want to analyze
grains_selected = grains;

% Minimum grain size (euivalent diameter) to be considered (in micron)
min_grain_size = 0.1;

% Number of bins in grain analysis
n_bins = 200;

condition = grains_selected.equivalentRadius > (min_grain_size / 2.0);
grain_area = grains_selected(condition).area;
total_area = sum(grain_area);
grain_size = 2.0 * grains_selected(condition).equivalentRadius;
n_grains = length(grain_size);

[bin_counts,bin_centers] = hist(grain_size,n_bins);
bin_counts = bin_counts';
bin_centers = bin_centers';
bin_width = bin_centers(3) - bin_centers(2);

bin_low = zeros(n_bins,1);
bin_high = zeros(n_bins,1);
for i = 1 : n_bins
    bin_low(i) = bin_centers(i) - bin_width / 2.0;
    bin_high(i) = bin_centers(i) + bin_width / 2.0;
end
max_grain_size = max(grain_size);
min_grain_size = min(grain_size);
bin_low(1) = min_grain_size;
bin_high(n_bins) = max_grain_size;

area_per_bin = zeros(n_bins,1);
area_counts_per_bin = zeros(n_bins,1);
for i = 1 : n_grains
      for j = 1 : n_bins
           if (grain_size(i) >= bin_low(j)) && (grain_size(i) <= bin_high(j))
              area_per_bin(j) = area_per_bin(j) + grain_area(i);
              area_counts_per_bin(j) = area_counts_per_bin(j) + 1;
           end
      end
end

bin_number_fraction = zeros(n_bins,1);
bin_number_fraction = 100.0 * bin_counts / n_grains;

mean_grain_size_number = sum(bin_centers .* bin_number_fraction / 100.0)

figure
h = bar(bin_centers,[bin_number_fraction],'histc','FaceColor','flat');
h(1).FaceColor = 'red';

l = 'Number fraction';

legend(h,l);
set(gca,'FontSize',15)

xlabel('Grain size (\it{d}\rm) [\mum]')
ylabel('Fraction [%]')

% Fitting a log-normal function to the histogram
% logn_fit = lognfit(grain_size);
% mu = logn_fit(1)
% sigma = logn_fit(2)
% x = 0:0.1:round(max_grain_size);
% hold on
% p = plot(x,100.0 * bin_width * lognpdf(x,mu,sigma))
% hold off
% set(p,'DisplayName','Lognormal fit','LineStyle','--','LineWidth',2.5)

clear min_grain_size n_bins condition grain_area total_area grain_size n_grains
clear bin_counts bin_centers bin_width bin_low bin_high area_per_bin area_counts_per_bin 
clear bin_area_fraction bin_number_fraction h l logn_fit x p i j grains_selected

%% Section 2.1.1 Grain size analysis (area) both modes superposed

% Number of bins in grain analysis
n_bins = 70;

grain_area_small = grains_small.area;
grain_area_large = grains_large.area;
grain_area = grains.area;
total_area = sum(grain_area);
grain_size = 2.0 * grains.equivalentRadius;
grain_size_small = 2.0 * grains_small.equivalentRadius;
grain_size_large = 2.0 * grains_large.equivalentRadius;
n_grains = length(grain_size);
n_grains_small = length(grain_size_small);
n_grains_large = length(grain_size_large);

[bin_counts,bin_centers] = hist(grain_size,n_bins);
bin_counts = bin_counts';
bin_centers = bin_centers';
bin_width = bin_centers(3) - bin_centers(2);

bin_low = zeros(n_bins,1);
bin_high = zeros(n_bins,1);
for i = 1 : n_bins
    bin_low(i) = bin_centers(i) - bin_width / 2.0;
    bin_high(i) = bin_centers(i) + bin_width / 2.0;
end

max_grain_size = max(grain_size);
min_grain_size = min(grain_size);
bin_low(1) = min_grain_size;
bin_high(n_bins) = max_grain_size;

area_per_bin_small = zeros(n_bins,1);
area_counts_per_bin_small = zeros(n_bins,1);
for i = 1 : n_grains_small
      for j = 1 : n_bins
           if (grain_size_small(i) >= bin_low(j)) && (grain_size_small(i) <= bin_high(j))
              area_per_bin_small(j) = area_per_bin_small(j) + grain_area_small(i);
              area_counts_per_bin_small(j) = area_counts_per_bin_small(j) + 1;
           end
      end
end

area_per_bin_large = zeros(n_bins,1);
area_counts_per_bin_large = zeros(n_bins,1);
for i = 1 : n_grains_large
      for j = 1 : n_bins
           if (grain_size_large(i) >= bin_low(j)) && (grain_size_large(i) <= bin_high(j))
              area_per_bin_large(j) = area_per_bin_large(j) + grain_area_large(i);
              area_counts_per_bin_large(j) = area_counts_per_bin_large(j) + 1;
           end
      end
end



bin_area_fraction_small = zeros(n_bins,1);
bin_area_fraction_small = 100.0 * area_per_bin_small / total_area;

bin_area_fraction_large = zeros(n_bins,1);
bin_area_fraction_large = 100.0 * area_per_bin_large / total_area;

mean_grain_size_area_small = sum(grain_area_small .* grain_size_small / sum(grain_area_small))
mean_grain_size_area_large = sum(grain_area_large .* grain_size_large / sum(grain_area_large))

figure
h = bar(bin_centers,[bin_area_fraction_small,bin_area_fraction_large],'histc','FaceColor','flat');
h(1,1).FaceColor = 'red';
h(1,2).FaceColor = 'blue';

l = cell(1,2);
l{1} = 'Area fraction small';
l{2} = 'Area fraction large';

legend(h,l);
set(gca,'FontSize',15)
% xticks(0:20:(bin_centers(n_bins)+ xticks_resolution));
% ylim([0 30])
% yticks(0:5:30);

xlabel('Grain size (\it{d}\rm) [\mum]')
ylabel('Fraction [%]')


% area
% mean_small = mean(droplet_diameter);
% variance= var(droplet_diameter);
% sigma= (variance)^0.5
%  
%  
% y_lognormal=lognpdf(x,log(mean),sigma)
% figure
% plot(x,y_lognormal)


% Fitting a log-normal function to the histogram
% logn_fit = lognfit(grain_area_large);
% mu = logn_fit(1)
% sigma = logn_fit(2)
% x = 0:0.1:round(max_grain_size);
% hold on
% p = plot(x,100.0 * bin_width * lognpdf(x,mu,sigma))
% hold off
% set(p,'DisplayName','Lognormal fit','LineStyle','--','LineWidth',2.5)

clear min_grain_size n_bins condition grain_area total_area grain_size n_grains
clear bin_counts bin_centers bin_width bin_low bin_high area_per_bin area_counts_per_bin 
clear bin_area_fraction bin_number_fraction h l logn_fit x p i j grains_selected

%% Section 2.2 Grain aspect ratio analysis 

% Select the grain object you want to analyze
grains_selected = grains_small;

% Minimum grain size (euivalent diameter) to be considered (in micron)
min_grain_size = 0.1;

% Number of bins in grain analysis
n_bins = 30;

condition = grains_selected.equivalentRadius > (min_grain_size / 2.0);
grain_area = grains_selected(condition).area;
total_area = sum(grain_area);
grain_size = 2.0 * grains_selected(condition).equivalentRadius;
n_grains = length(grain_size);
aspect_ratio = 1.0 ./ grains_selected(condition).aspectRatio;

[bin_counts,bin_centers] = hist(aspect_ratio,n_bins);
bin_counts = bin_counts';
bin_centers = bin_centers';
bin_width = bin_centers(3) - bin_centers(2);

bin_low  = zeros(n_bins,1);
bin_high = zeros(n_bins,1);
for i = 1 : n_bins
    bin_low(i) = bin_centers(i) - bin_width / 2.0;
    bin_high(i) = bin_centers(i) + bin_width / 2.0;
end
bin_low(1) = 0.0;
bin_high(n_bins) = 1.0;

area_per_bin = zeros(n_bins,1);
area_counts_per_bin = zeros(n_bins,1);
for i = 1 : n_grains
      for j = 1 : n_bins
           if (aspect_ratio(i) >= bin_low(j)) && (aspect_ratio(i) <= bin_high(j))
              area_per_bin(j) = area_per_bin(j) + grain_area(i);
              area_counts_per_bin(j) = area_counts_per_bin(j) + 1;
           end
      end
end

bin_area_fraction = zeros(n_bins,1);
bin_area_fraction = 100.0 * area_per_bin / total_area;
bin_number_fraction = zeros(n_bins,1);
bin_number_fraction = 100.0 * bin_counts / n_grains;

mean_aspect_ratio_area   = sum(bin_centers .* bin_area_fraction / 100.0)
mean_aspect_ratio_number = sum(bin_centers .* bin_number_fraction / 100.0)
effective_aspect_ratio = (mean_aspect_ratio_area + mean_aspect_ratio_number) / 2.0

figure
h = bar(bin_centers,[bin_area_fraction,bin_number_fraction],'histc','FaceColor','flat');
h(1,1).FaceColor = 'red';
h(1,2).FaceColor = 'blue';

l = cell(1,2);
l{1} = 'Area fraction';
l{2} = 'Number fraction';

legend(h,l);
set(gca,'FontSize',15)
xticks(0:0.1:1);
% ylim([0 35])
% yticks(0:5:35);

xlabel('Grain aspect ratio (\itm\rm) [-]')
ylabel('Fraction [%]')

% Fitting a normal function to the histogram
[mu,sigma] = normfit(aspect_ratio)
x = 0:0.01:1;
hold on
p = plot(x,100.0 * bin_width * normpdf(x,mu,sigma))
hold off
set(p,'DisplayName','Normal fit','LineStyle','--','LineWidth',2.5)



clear min_grain_size n_bins condition grain_area total_area grain_size n_grains aspect_ratio
clear bin_counts bin_centers bin_width bin_low bin_high area_per_bin area_counts_per_bin 
clear bin_area_fraction bin_number_fraction h l logn_fit x p i j grains_selected

%% Section 2.3 Grain shape (ellipse axis) orientation analysis

% Select the grain object you want to analyze
grains_selected = grains_small;

% Minimum grain size (euivalent diameter) to be considered (in micron)
min_grain_size = 0.1;

% Number of bins in grain analysis
n_bins = 30;

condition = grains_selected.equivalentRadius > (min_grain_size / 2.0);
grain_area = grains_selected(condition).area;
total_area = sum(grain_area);
grain_size = 2.0 * grains_selected(condition).equivalentRadius;
n_grains = length(grain_size);

[omega,~,~] = fitEllipse(grains_selected(condition));
omega = omega / degree;
theta = zeros(n_grains,1);
for i = 1:n_grains
    if (omega(i) >= 0.0) && (omega(i) <=90.0)
        theta(i) = 90.0 - omega(i);
    else
        theta(i) = 270.0 - omega(i);
    end
end


[bin_counts,bin_centers] = hist(theta,n_bins);
bin_counts = bin_counts';
bin_centers = bin_centers';
bin_width = bin_centers(3) - bin_centers(2);

bin_low  = zeros(n_bins,1);
bin_high = zeros(n_bins,1);
for i = 1 : n_bins
    bin_low(i) = bin_centers(i) - bin_width / 2.0;
    bin_high(i) = bin_centers(i) + bin_width / 2.0;
end
bin_low(1) = 0.0;
bin_high(n_bins) = 180.0;

area_per_bin = zeros(n_bins,1);
area_counts_per_bin = zeros(n_bins,1);
for i = 1 : n_grains
      for j = 1 : n_bins
           if (theta(i) >= bin_low(j)) && (theta(i) <= bin_high(j))
              area_per_bin(j) = area_per_bin(j) + grain_area(i);
              area_counts_per_bin(j) = area_counts_per_bin(j) + 1;
           end
      end
end

bin_area_fraction = zeros(n_bins,1);
bin_area_fraction = 100.0 * area_per_bin / total_area;
bin_number_fraction = zeros(n_bins,1);
bin_number_fraction = 100.0 * bin_counts / n_grains;

mean_angle_area   = sum(bin_centers .* bin_area_fraction / 100.0)
mean_angle_number = sum(bin_centers .* bin_number_fraction / 100.0)
effective_angle = (mean_angle_area + mean_angle_number) / 2.0

figure
h = bar(bin_centers,[bin_area_fraction,bin_number_fraction],'histc','FaceColor','flat');
h(1,1).FaceColor = 'red';
h(1,2).FaceColor = 'blue';

l = cell(1,2);
l{1} = 'Area fraction';
l{2} = 'Number fraction';

legend(h,l);
set(gca,'FontSize',15)
xticks(0:20:180);
% ylim([0 35])
% yticks(0:5:35);

xlabel('Grain major axis angle (\it\omega\rm) [ {\circ} ]')
ylabel('Fraction [%]')

% Fitting a normal function to the histogram
[mu,sigma] = normfit(theta)
x = 0:0.1:180;
hold on
p = plot(x,100.0 * bin_width * normpdf(x,mu,sigma))
hold off
set(p,'DisplayName','Normal fit','LineStyle','--','LineWidth',2.5)

clear min_grain_size n_bins condition grain_area total_area grain_size n_grains aspect_ratio
clear bin_counts bin_centers bin_width bin_low bin_high area_per_bin area_counts_per_bin 
clear bin_area_fraction bin_number_fraction h l logn_fit x p i j omega theta grains_selected

%% 3.1 Separate ebsd data to large and small regions

ebsd_large = ebsd(grains_large);
ebsd_small = ebsd(grains_small);

%% Define the pseudo grain boundaries for visualization
% Boundary is a pseudo boundary if it satisfies one of the two conditions
condition1 = ismember(grains.boundary.grainId(:,1),large_grains(:,2)) + ismember(grains.boundary.grainId(:,2),small_grains(:,2));
condition2 = ismember(grains.boundary.grainId(:,1),small_grains(:,2)) + ismember(grains.boundary.grainId(:,2),large_grains(:,2));
x1 = find(condition1 == 2);
x2 = find(condition2 == 2);
x = cat(1,x1,x2);
pseudo_grainboundary = grains.boundary(x);

clear condition1 condition2 x1 x2 x

%% 3.2 Texture analysis large grains

ebsd_selected = ebsd;

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

ipfKey = ipfColorKey(ebsd_selected('Aluminium'));
ipfKey.inversePoleFigureDirection = vector3d.X;

colors = ipfKey.orientation2color(ebsd_selected('Aluminium').orientations);
figure
plot(ebsd_selected('Aluminium'),colors)

hold on
plot(grains_large.boundary,'lineWidth',1)
plot(grains_small.boundary,'lineWidth',1)
plot(pseudo_grainboundary,'linewidth',3)
hold off


%% Plot ipf
odf = calcODF(ebsd_small('Aluminium').orientations);

% Define PF directions and plot PF
% h = [Miller(1,0,0,cs) Miller(0,1,0,cs) Miller(0,0,1,cs)];

figure
plotIPDF(odf,[xvector,yvector,zvector],'antipodal','contours',40)
mtexColorbar
CLim(gcm,[0.85,1.25])

clear ebsd_selected

%% Plot aspect ratio colouring

% normalize aspect ratios
large_aspectRatio = 1 ./ grains_large.aspectRatio;
small_aspectRatio = 1 ./ grains_small.aspectRatio;

% plot grains wrt aspect ratios
figure
plot(grains_large,large_aspectRatio);
hold on
plot(grains_small,small_aspectRatio);
CLim(gcm,[0 1]);
mtexColorMap hot
mtexColorbar
hold off
hold on
% superpose pseudo grain boundaries
plot(pseudo_grainboundary,'linewidth',3)
hold off

