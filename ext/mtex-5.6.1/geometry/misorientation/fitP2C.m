function p2c = fitP2C(c2c,p2c0)

% maybe c2c was list of pairs of orientations
if ~isa(c2c.SS,'crystalSymmetry')
  c2c = inv(c2c(:,1)) .* c2c(:,2);
end


radius = 10*degree;

p2c = p2c0;

for k = 1:5

  p2cGrid = localOrientationGrid(p2c.CS,p2c.SS,radius,'center',p2c,'resolution',radius/3);

  fit = calcP2CFit(p2cGrid,c2c);
  
  [~,i] = min(fit); fit(i)./degree
  
  p2c = p2cGrid(i);
  
  radius = radius / 2;
  
end

end