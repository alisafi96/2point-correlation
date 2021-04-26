%% Import Script for EBSD Data
%
% This script was automatically created by the import wizard. You should
% run the whoole script or parts of it in order to import your data. There
% is no problem in making any changes to this script.

%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('321', [4.9 4.9 5.4], 'X||a', 'Y||b*', 'Z||c*', 'mineral', 'Quartz', 'color', [0.53 0.81 0.98])};

% plotting convention
% setMTEXpref('xAxisDirection','west');
% setMTEXpref('zAxisDirection','outOfPlane');
plotx2east;
plotzOutOfPlane;

%% Specify File Names

% path to files
pname = '/Users/ali/Documents/GitHub/2point-correlation/data';

% which files to be imported
fname = [pname '/test.txt'];

%% Import the Data

% create an EBSD variable containing the data
% ebsd = EBSD.load(fname,CS,'interface','ang',...
%  'convertEuler2SpatialReferenceFrame');

ebsd = EBSD.load(fname,'cs',CS,'ColumnNames',{'x','y','Euler1','Euler2','Euler3','phase'});

%% Calculate WignerD/GSH functions

D = WignerD(ebsd('Quartz').orientations);