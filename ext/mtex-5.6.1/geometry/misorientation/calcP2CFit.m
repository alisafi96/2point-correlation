function fit = calcP2CFit(p2c,mori)

fit = zeros(size(p2c));

n = round(0.7 * length(mori));

for k = 1:length(p2c)

  % child to child misorientation variants
  p2cV = p2c(k).variants; p2cV = p2cV(:);
  
  c2c = p2c(k) * inv(p2cV);
  
  % misorientation to c2c variants
  omega = angle_outer(mori, c2c);
  
  % compute best fitting variant
  omega = sort(min(omega,[],2));
  
  fit(k) = mean(omega(1:n));
  
end

end
