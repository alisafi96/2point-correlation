%% Direct evaluation of ODF
bandwidth = 10;

C = sum(conj(WignerD(ori, 'bandwidth', bandwidth)),2);
C = C/length(ori);

%%
l = 0:bandwidth;
d = (2*l+1).^2;
cs = [0 cumsum(d)];

for i = 1:length(C)
    
    % Check which harmonic degree l is currently evaluated
    if ismember(i, cs)
        [~,l] = ismember(i, cs);
        l = l-2;
    end
        
    C(i) = C(i) * (2*l+1);
end


%% Calculate and compare odfs

odf_kernel = calcKernelODF(ori,'halfwidth',10*degree, 'bandwidth', bandwidth);
component = odf_kernel.components{1, 1};

% Harmoic expansion of the unimodal odf
odf = calcDensity(ori,'Fourier','bandwidth', bandwidth);

% Direct calculation
odf_C = FourierODF(f_hat, ori.CS);

%% Calculate fhats from unimodal ODF using NFSOFT

g = quaternion(component.center);
c = component.weights;
f_hat = gcA2fourier(g,c,A);

%% Test evaluations

f = eval(odf,ori(2500))
f = eval(odf_kernel,ori(2500))
f = eval(odf_C,ori(2500))

%% Calculate ODF from Unimodal orientations

C_uni = sum(conj(WignerD(odf_kernel.components{1, 1}.center, 'bandwidth', bandwidth)),2);
C_uni = C_uni/length(odf_kernel.components{1, 1}.center);

odf_uni = FourierODF(C_uni, ori.CS);
%% ODF calculation from harmonic expansion of kernel function

% Specify set of orientations from unimodal form
ori_uni = orientation.byEuler(component.center.phi1,component.center.Phi,component.center.phi2,component.CS);

% Compute complex coefficients in complex form
C = conj(WignerD(ori_uni, 'bandwidth', bandwidth));

% Multiplication of coefficients with weights
for i = 1:length(component.weights)
C(:,i) = component.weights(i,1) * C(:,i);
end

% Apply summation of 
C = sum(C,2);

% Multiply with Chebyshev coefficients

l = 0:bandwidth;
d = (2*l+1).^2;
cs = [0 cumsum(d)];

for i = 1:length(C)
    
    % Check which harmonic degree l is currently evaluated
    if ismember(i, cs)
        [~,l] = ismember(i, cs);
        l = l-2;
    end
        
    C(i) = C(i) * component.psi.A(l+1);
end

%% Evaluate Fourier coefficients
figure
plotFourier(odf, 'bandwidth', bandwidth)

figure
plotFourier(odf_C, 'bandwidth', bandwidth)

figure
plotFourier(odf_kernel, 'bandwidth', bandwidth)
%% Evaluating odfs using pole figure plots

h = { ...
  Miller(0,0,0,2,mori.CS),...
  Miller(1,0,-1,0,mori.CS),...
  Miller(1,0,-1,1,mori.CS),...
  Miller(1,0,-1,2,mori.CS),...
  };

figure
plotPDF(odf,h,'antipodal');

mtexColorbar
%CLim(gcm,[0.0,2])

figure
plotPDF(odf_C,h,'antipodal');

mtexColorbar
%CLim(gcm,[0.0,2])

figure
plotPDF(odf_kernel,h,'antipodal');

mtexColorbar


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

%% Function

function c_hat = gcA2fourier(g,c,A)

% 2^4 -> nfsoft-represent
% 2^2 -> nfsoft-use-DPT
nfsoft_flags = bitor(2^4,4);
plan = nfsoftmex('init',length(A)-1,length(g),nfsoft_flags,0,4,1000,2*ceil(1.5*(length(A)+1)));
nfsoftmex('set_x',plan,Euler(g,'nfft').');
nfsoftmex('set_f',plan,c(:));
nfsoftmex('precompute',plan);
nfsoftmex('adjoint',plan);
c_hat = nfsoftmex('get_f_hat',plan);
nfsoftmex('finalize',plan);

for l = 1:length(A)-1
  ind = (deg2dim(l)+1):deg2dim(l+1);
  c_hat(ind) = A(l+1)* reshape(c_hat(ind),2*l+1,2*l+1);
end
  
end
