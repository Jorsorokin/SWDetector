function [bounds,dsplit] = splitdata_overlap(d,win,overlap)
% [bounds,dsplit] = splitdata_overlap(d,win,overlap)
%
% Splits vector/matrix "d" into segments of length "win", overlapping
% with the next by "overlap" amount of points. Note that "win" and
% "overlap" should be given in samples
%
% By JMS, 10/31/16

% get number of segments 
if numel(d) == 1
    N = d;
else
    [N,M] = size(d);
end
nseg = floor( N / (win-overlap) ); % since win + (overlap*nseg) = N

% preallocate matrix for storing the split data and indexing bounds
bounds(:,1) = (win-overlap)*((1:nseg)-1) + 1;
bounds(:,2) = bounds(:,1) + win-1;
bounds(bounds(:,2)>N,:) = [];
nseg = size( bounds,1 );
if nargout > 1
    dsplit = zeros(win,nseg,M);
    
    % loop and take out overlapped data segments
    for seg = 1:nseg
        dsplit(:,seg,:) = d(bounds(seg,1):bounds(seg,2),:);
    end
end

end