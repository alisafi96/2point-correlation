function job = revert(job,ind)
% revert parent orientations to measured child orientations
%
% Syntax
%
%   % undo all reconstruction
%   job.revert
%
%   % undo grains by condition
%   ind = job.grains.fit > 5*degree
%   job.revert(ind)
%
%   % undo specific grains
%   job.revert(job.parentGrains(5))
%
% Input
%  job    - @parentGrainReconstructor
%  ind    - true/false of which parent grains should be reverted to child grains
%  grains - list of grains to be reverted
%
% Output
%  job - @parentGrainReconstructor
%
% See also
% MaParentGrainReconstruction
%

% revert everything
if nargin == 1, ind = true(size(job.grains)); end

% input should be index to job.grains
if isa(ind,'grain2d'), ind = id2ind(job.grains,ind.id); end

% make ind a logical vector
if ~islogical(ind)
  ind = full(sparse(ind,1,true,length(job.grains),1));
end

% which are the original grains to revert
% doRevert = ismember(job.mergeId, ind);
doRevert = ind(job.mergeId) & job.isTransformed;

if ~any(doRevert), return; end

% we may have two situations
% 1. grains to revert have not yet been merged
% 2. grains has already been merged

% no merge needs to be undone
if max(accumarray(job.mergeId(doRevert),1)) == 1

  % reset phase and rotation
  job.grains.phaseId(doRevert) = job.grainsMeasured.phaseId(doRevert);
  job.grains.prop.meanRotation(doRevert) = job.grainsMeasured.prop.meanRotation(doRevert);

elseif 1 % we do not undo the merge but we redo the merge of the remaining grains
  
  newMergeId = job.mergeId;
  
  % grains to revert a put first
  newMergeId(doRevert) = 1:nnz(doRevert);
  
  % shift everything that we do not touch at the end
  shift = nnz(doRevert) - cumsum(ind);
  newMergeId(~doRevert) = newMergeId(~doRevert) + shift(newMergeId(~doRevert));
  
  % remerge grains
  [grains, mergeId] = merge(job.grainsMeasured, newMergeId);
  
  
  % reassign phase and rotation to the reverted grains
  grains.phaseId(1:nnz(doRevert)) = job.grainsMeasured.phaseId(doRevert);
  grains.prop.meanRotation(1:nnz(doRevert)) = job.grainsMeasured.prop.meanRotation(doRevert);
  
  % reassign phase and rotation to the not reverted grains
  grains.phaseId((nnz(doRevert)+1):end) = job.grains.phaseId(~ind);
  grains.prop.meanRotation((nnz(doRevert)+1):end) = job.grains.prop.meanRotation(~ind);
  
  job.grains = grains;
  job.mergeId = mergeId;

end

% update grain boundaries and triple points
job.grains = job.grains.update;

