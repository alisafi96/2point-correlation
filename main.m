%% Import data
mtexdata titanium;

%% Calculate WignerD/GSH functions

%odf_kernel = calcKernelODF(ebsd.orientations,'halfwidth',10*degree);

%odf_psi = calcDensity(ebsd.orientations,'kernel',psi);

%% Building the source Wigner D matrix to construct Fourier coefficients

ori = ebsd(1:100).orientations;
bandwidth = 20;

C = sum(conj(WignerD(ori, 'bandwidth', bandwidth)),2);
C = C/length(ori);

%% Evaluation of ODF (equation 4.21) 

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
