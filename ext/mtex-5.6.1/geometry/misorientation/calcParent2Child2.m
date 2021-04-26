function [p2c, omega] = calcParent2Child2(mori,p2c,alpha,varargin)
%
% Syntax
%
%   p2c = calcParent(childOri,p2c)
%
%   p2c = calcParent(c2c,p2c)
%   p2c = calcParent(c2c,p2c,alpha)
%
%   [p2c, fit] = calcParent(c2c,p2c)
%
% Input
%  childOri - child @orientation
%  c2c      - child to child mis@orientation
%  p2c      - initial gues or list of measured parent to child boundary mis@orientations
%  alpha    - damping factor ( default - 0)
%
% Output
%  p2c      - parent to child mis@orientation
%  fit      - disorientation angle between all c2cs misorientations and the computed one
%
% Options
%  maxIteration - maximum number of iterations (default - 100)
%  threshold    - consider only misorientation within the threshold to the initial p2c (default - 5*degree)
%
% References
%
% * Tuomo Nyyssönen, <https://www.researchgate.net/deref/http%3A%2F%2Fdx.doi.org%2F10.1007%2Fs11661-018-4904-9?_sg%5B0%5D=gRJGzFvY4PyFk-FFoOIj2jDqqumCsy3e8TU6qDnJoVtZaeUoXjzpsGmpe3TDKsNukQYQX9AtKGniFzbdpymYvzYwhg.5jfOl5Ohgg7pW_6yACRXN3QiR-oTn8UsxZjTbJoS_XqwSaaB7r8NgifJyjSES2iXP6iOVx57sy8HC4q2XyZZaA
% Crystallography, Morphology, and Martensite Transformation of Prior
% Austenite in Intercritically Annealed High-Aluminum Steel>
%

% compute misorientations if pairs of orientations are given
if isa(mori.SS, 'specimenSymmetry'), mori = inv(mori(:,1)) .* mori(:,2); end

% third input is damping factor
if nargin<3, alpha = length(p2c) > 1; end

threshold = get_option(varargin,'threshold',5*degree);

% if p2c is a list of parent2child orientation relationships - take the mean first
if length(p2c) > 1
  alpha = alpha * length(p2c)/length(mori);
  p2c = mean(p2c,'robust');
  p2c0 = p2c;
end

% prepare iterative loop
maxIt = get_option(varargin,'maxIterarion',100);
diso = nan(maxIt,1);
p2cOld = p2c;

csChild = p2c.SS.properGroup;

% iterate until convergance
for k = 1:maxIt
  
  %progress(k,maxIt,'searching orientation relationship: ');
  diso(k) = angle(p2c,p2cOld)/degree;
  p2cOld = p2c;
    
  %check for convergence
  if k>5 && median(diso(k-5:k)) < 0.5, break; end
  
  % child to child misorientation variants
  p2cV = p2c.variants; p2cV = p2cV(:);
  c2c = p2c * inv(p2cV); %#ok<MINV>
  
  % misorientation to c2c variants
  omega = angle_outer(mori, c2c);

  % compute best fitting variant
  [omega, variant] = min(omega,[],2);
    
  % take only those c2c misorientations that are suffiently close to the
  % current candidate
  ind = omega(:) < threshold;
  
  %
  fit(k) = mean(omega(ind));
  fit(k)
  
  % c2c for the best variant
  c2cV = p2c .* inv(p2cV(variant(ind)));
  
  % next we determine the best fitting child symmetry
  omega = zeros(nnz(ind),numSym(csChild));
  for ics = 1:numSym(csChild)
    
    omega(:,ics) = min(angle_outer(mori(ind) .* csChild.rot(ics) .* inv(c2cV), ...
      csChild.rot,'noSymmetry'),[],2);    
  end
  
  [~,ics] = min(omega,[],2);
  
  %weigthts = [nnz(ind);ones(nnz(ind),1)];
  %weigthts = [100000;ones(nnz(ind),1)];
  %p2c = mean([p2c;inv(mori(ind)) .* inv(csChild.rot(ics)) .* p2cV(variant(ind))],...
  %  'weights',weigthts);
  p2c = mean(inv(mori(ind)) .* inv(csChild.rot(ics)) .* p2cV(variant(ind)));
  
end