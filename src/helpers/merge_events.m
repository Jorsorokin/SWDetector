function events = merge_events( x,L )
    % events = merge_events( x,L )
    %
    % merges events in x that are closer together than "gap" points
    %
    % inputs:
    %   x - an n x 2 matrix, where 1st column is start, second is end,
    %       or an n x 1 vector with starting times. In both situations,
    %       rows of x must be monotonically increasing
    %
    %   L - maximum # of time allowed to combine events
    %
    % outputs:
    %   events - a new vector of merged event times
    %
    % by jordan sorokin, 2/3/19
    
    if size( x,2 ) == 2
        removeBool = [false; (x(2:end,1) - x(1:end-1,2)) <= L];
    else
        removeBool = [false; diff( x ) <= L];
    end
    
    within = [0;diff( removeBool )];
    start = (within~=1 & ~removeBool);
    stop = [xor( start(1:end-1), diff(start) ); true];
    
    events = [x(start,1),x(stop,2)];
end
        
    