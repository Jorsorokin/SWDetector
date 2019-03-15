function events = remove_short_events( x,L )
    % events = remove_short_events( x,L )
    %
    % removes short events from a list of event times
    %
    % inputs:
    %   x - an n x 2 matrix or n x 1 vector
    %       - if matrix, then first column is start time, second is end
    %       - if vetor, each element represents an event duration
    %
    %   L - scalar minimum event length
    %
    % outputs:
    %   events - m x 2 or m x 1 vector, with m <= n
    %
    % by jordan sorokin, 2/3/19
    
    if size( x,2 ) == 1
        events = x(x >= L);
    else
        d = x(:,2) - x(:,1);
        events = x(d >= L,:);
    end
end