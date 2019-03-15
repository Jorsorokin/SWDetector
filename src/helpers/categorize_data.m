function [events,bounds] = categorize_data( data,fs,window,varargin )
    % [events,bounds] = categorize_data( data,fs,window,(start,stop) )
    %
    % manually mark all windows of data with a specified tag (number)
    %
    % inputs:
    %   data - vector of eeg data
    %   fs - sampling rate
    %   window - window of time to segment data into (in sec)
    %   (start) - the starting time to begin searching
    %   (stop) - the ending time to begin searching
    %
    % outputs:
    %   events - an n-length vector of event classes (categories) defined by the user
    %   bounds - an nx2 matrix specifying the starting/ending times of each event
    %               - note, the bounds are relative to the start of the
    %               data, so even if you supply a different "start" point,
    %               the bounds will add this start time to the first event
    
    
    % get length for the data segments
    totalTime = length(data)/fs;
    x = 1/fs:1/fs:totalTime;
    
    % check optionals
    if nargin > 3 && ~isempty(varargin{1})
        start = varargin{1};
        startpt = round(start*fs);
    else
        start = x(1);
        startpt = 1;
    end
    
    if nargin > 4 && ~isempty(varargin{2})
        stop = varargin{2};
        stoppt = round(stop*fs);
    else
        stop = x(end);
        stoppt = length(x);
    end
    
    nSegs = floor( (stop-start) / window );
    nPts = floor( window*fs );
    
    % begin plots for visualization
    mu = mean(data(startpt:stoppt)); sd = std(data(startpt:stoppt));
    yl = [mu-3*sd, mu+3*sd];
    
    figure;
    ax1 = subplot(2,1,1); 
    plot( ax1,x(startpt:stoppt),data(startpt:stoppt),'k' );
    set( ax1,'nextplot','add','ylim',yl,'tickdir','out','box','off' );
    
    ax2 = subplot(2,1,2);
    set( ax2,'nextplot','add','ylim',yl,'xlim',[0 nPts] );
    
    events = zeros(1,nSegs);
    bounds = zeros(nSegs,2);
    
    for i = 1:nSegs
        t = ((i-1)*nPts+1:i*nPts) + startpt;
        bounds(i,:) = [t(1),t(end)];
        minX = x(t(1));
        maxX = x(t(end));
        
        l = plot( ax1,x(t),data(t),'r' );
        set( ax1,'xlim',[max(minX-window*20,x(1)), min(maxX+window*20,x(end))] );
        
        plot( ax2,data(t),'r' );
        events(i) = mark_event();
        
        cla(ax2); delete(l);
    end
    
    
    function class = mark_event()
        while true
            waitforbuttonpress;
            class = str2double( get(gcf,'currentcharacter') );
            if ~isnan(class)
                break
            end
        end
    end
end
    
    