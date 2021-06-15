%% Import data
mtexdata titanium;

%% Building the source Wigner D matrix to construct Fourier coefficients

% Define set of orientations & bandwidth
ori = ebsd.orientations;
bandwidth = 15;

C = sum(conj(WignerD(ori, 'bandwidth', bandwidth)),2);
C = C/length(ori);

%% Calculate and compare odfs

%odf_kernel = calcKernelODF(ebsd.orientations,'halfwidth',10*degree);

% Harmoic expansion of the unimodal odf
odf = calcDensity(ori,'Fourier','bandwidth', bandwidth, 'exact');

% Direct calculation
odf_C = FourierODF(C, ori.CS);


%% Evaluate Fourier coefficients
figure
plotFourier(odf, 'bandwidth', bandwidth)

figure
plotFourier(odf_C, 'bandwidth', bandwidth)



%% Evaluating odfs using pole figure plots

h = { ...
  Miller(0,0,0,2,ori.CS),...
  Miller(1,0,-1,0,ori.CS),...
  Miller(1,0,-1,1,ori.CS),...
  Miller(1,0,-1,2,ori.CS),...
  };

figure
plotPDF(odf,h,'antipodal');

mtexColorbar
%CLim(gcm,[0.0,2])

figure
plotPDF(odf_C,h,'antipodal');

mtexColorbar
%CLim(gcm,[0.0,2])


%% Calculate Misorientations

% Uncorrelated mdfs
mdf = calcMDF(odf, 'bandwidth', bandwidth);
mdf_C = calcMDF(odf_C);


figure
plotAxisDistribution(mdf)
figure
plotAxisDistribution(mdf_C)
%% Calculate correlated MDF

%grains = calcGrains(ebsd,'angle',10*degree);

% Find grain boundary misorientation
mori = grains.boundary('Titanium (Alpha)','Titanium (Alpha)').misorientation;

mdf_correlated  = calcMDF(mori,'bandwidth',bandwidth);

C_mori = sum(conj(WignerD(mori, 'bandwidth', bandwidth)),2);
C_mori = C_mori/length(mori);

mdf_C_correlated = FourierODF(C_mori,mori.CS,mori.SS,'antipodal');

%% Plot Axis Distribution

figure
plotAxisDistribution(mdf_correlated);

figure
plotAxisDistribution(mdf_C_correlated);


%% Calculate Kernel

% psi = calcKernel(ori);
odf_kernel = calcKernelODF(ebsd.orientations,'halfwidth',5*degree);
%% Evaluation of ODF for a single orientation (can be replaced by MTEX function)



% Wigner D of orientation of interest
% Evaluate for random orientation
D_g = WignerD(ebsd.orientations(20), 'bandwidth', bandwidth);

l = 0:bandwidth;
d = (2*l+1).^2;
cs = [0 cumsum(d)];

ODF_sum = 0;

for i = 1:length(C)
    
    % Check which l is currently evaluated
    if ismember(i, cs)
        [~,l] = ismember(i, cs);
        l = l-2;
    end
    
    element = (2*l+1) * C(i) * D_g(i);
    
    ODF_sum = ODF_sum + element;
end
