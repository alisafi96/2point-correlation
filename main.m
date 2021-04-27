%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = crystalSymmetry('m-3m', [3.7 3.7 3.7], 'mineral', 'Iron');

% plotting convention
% setMTEXpref('xAxisDirection','east');
% setMTEXpref('zAxisDirection','outOfPlane');

plotx2east;
plotzOutOfPlane;

%% Specify File Names

% path to files
pname = '\\juno\Homes\user\asafi\2point-correlation\data';

% which files to be imported
fname = [pname '\LMD_EBSD.txt'];

%% Import the Data

% create an EBSD variable containing the data
ebsd = EBSD.load(fname,CS,'interface','generic',...
  'ColumnNames', { 'phi1' 'Phi' 'phi2' 'x' 'y'}, 'Bunge');

%%
plot(ebsd,ebsd.orientations)

%% Calculate WignerD/GSH functions

% bandwidth or degree
l = 1;
% D = WignerD(ebsd.orientations(1),'degree',l);
% D = WignerD(ebsd.orientations(1),'degree',l);
D = WignerD(ebsd.orientations(1));