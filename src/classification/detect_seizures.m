function S = detect_seizures( data,fs,varargin )
    % S = detect_seizures( data,fs,(pThresh,minLength,maxGap,W,classifier) )
    %
    % detects seizures by first detecting individual spike-waves, then
    % merging nearby spike-waves into "events", and finally rejecting 
    % the resulting events that are either (a) too short, or (b) have an overall
    % probablity of being a seizure less than the specified probability threshold
    %
    %   * Note: call optional inputs using name-value format
    % 
    % inputs: 
    %   data - N x M data matrix
    %
    %   fs - the sampling rate
    %
    %   (pThresh) - probability threshold for detecting swd, between [0, 1] 
    %               (default = mean( model.p1 ) - std( model.p1 )
    %
    %   (minLength) - minimum length for an event to be considered a valid seizure 
    %                 (default = 1 second)
    %
    %   (maxGap) - maximum amount of time between neighboring swds to
    %              consdier them as part of a single event
    %              (default = 1 second)
    %
    %   (refine) - boolean; if true, refines seizures by removing those
    %              with median probability < pThresh (default = true)
    %
    %   (classifier) - the filename containing the classifier params to use
    %                  (default = '/shared/Jordan/Seizure_detection/swd_detection/models/norm-v-swd.mat')
    %
    %   (W) - pre-computed wavelet coefficients. Note: if using the default
    %         classifier, these MUST be computed as follows:
    %           1) zdata = zscore( data );
    %               1a) resample data to 500 Hz if not already
    %           2) W = modwt( zdata,'sym4',8 );
    %   
    % outputs:
    %   S - structure containing following fields:
    %           start - starting times of seizures
    %           end - ending times of seizures
    %           P - probability of being a seizure (average P(swd) over entire seizure)
    %           X - features used (relative wavelet variance) for each event
    %           V - matrix of data snips, aligned to seizure onset
    %
    % By Jordan Sorokin, 2/3/19
    
    [options,model] = check_inputs( varargin );
    window = floor( model.window*model.fs / 2 ); % could probably be optimized
    overlap = floor( model.window*model.fs / 4 ); % ditto
    
    % resample if necessary
    if fs ~= model.fs
        den = gcd( fs,model.fs );
        data = resample( data, model.fs/den, fs/den );
    end
    
    % split X and compute features over each chunk of data
    if isempty( options.W )
        options.W = modwtDecomp( zscore( data ),model.fs,'sym4',8 );
    else
        assert( size( options.W,2 ) == size( data,1 ),'# of columns of W must equal # of rows of X' );
        assert( size( options.W,3 ) == size( data,2 ),'# of channels in X and W do not match' );
    end
    
    bounds = splitdata_overlap( data,window,overlap );
    X = zeros( size( bounds,1 ),model.nCoeffs );
    
    for m = 1:size( options.W,3 )
        [~,w] = splitdata_overlap( options.W(:,:,m)',window,overlap );
        v = squeeze( var( w,1 ) );
        if model.useSQRT
            v = sqrt( v );
        end
        
        if model.useInteractions
            v = quadratic_expansion( v );
        end
        
        X = X + (v ./ sum( v,2 ));
        X = X + v;
    end

    % classify events based on model params
    X = X ./ m; % average wavelet relative variance over channels
    [P,swdBool] = classify_event( X,[model.B0;model.B],options.pThresh );
    
    if nnz( swdBool ) == 0
        disp( 'no detected events' );
        S = [];
        return
    end
    
    % detects regions of of high probability using cubic spline smoothing
    % to remove noise from sparse classification
    x = bounds/model.fs;
    Pspline = csaps( x(:,1),P,0.08 );
    fit = sum( Pspline.coefs,2 );
    [~,loc] = findpeaks( fit,'minpeakheight',median( fit ) + mad( fit ),...
                             'minpeakdistance',round( options.maxGap/2*model.fs ) );
    events = zeros( numel( loc ),2 );
    for i = 1:numel( loc )
        events(i,:) = search_probability_space( x,loc(i) ); % search from inside-out for each event
    end
        
    if size( events,1 ) > 1
        events = merge_events( events,options.maxGap );
    end
    
    events = remove_short_events( events,options.minLength ); 
    n = size( events,1 );
    if n == 0
        disp( 'no detected events with length > %0.2f\n',options.minLength );
        S = [];
        return
    end
    
    % only store events with mean probability > pThresh
    %   ... although some events may have had a few data chunks with P > pThresh, 
    %       (which were then merged and considered, together, an event),
    %       the overall probability of the entirety of the event may be smaller ...
    keptStart = find( ismember( x(:,1),events(:,1) ) );
    keptEnd = find( ismember( x(:,2),events(:,2) ) );
    bins = 0:0.2:1;
    expectedProb_g0 = histcounts( model.p0,bins,'normalization','probability' );
    expectedProb_g1 = histcounts( model.p1,bins,'normalization','probability' );

    p = cell( n,1 );
    x = cell( n,1 );
    dKL = zeros( n,2 );
    v = cell( n,1 );
    w = cell( n,1 );
    for i = 1:n
        observedProb = histcounts( P(keptStart(i):keptEnd(i)),bins,'normalization','probability' );
        dKL(i,1) = kldiv( observedProb,expectedProb_g0 );
        dKL(i,2) = kldiv( observedProb,expectedProb_g1 );     
        p{i} = P(keptStart(i):keptEnd(i)); 
        x{i} = X(keptStart(i):keptEnd(i),:);
        v{i} = data(bounds(keptStart(i),1):bounds(keptEnd(i),2),:);
        w{i} = options.W(:,bounds(keptStart(i),1):bounds(keptEnd(i),2),:);
    end
    
    if options.refine
        likelySzr = cell2mat( cellfun( @(x)(median(x) > options.pThresh),p,'un',0 ) );
    else
        likelySzr = true( 1,length( p ) );
    end
    
    S = struct( 'fs',fs,...
                'nEvents',nnz(likelySzr),...
                'pThresh',options.pThresh,...
                'start',events(likelySzr,1),...
                'end',events(likelySzr,2),...
                'dKL',dKL(likelySzr,:),...
                'params',rmfield(options,'W') );
    S.P = p(likelySzr);
    S.X = x(likelySzr);
    S.V = v(likelySzr);
    S.W = w(likelySzr);
   
    %% helper function
    function [options,model] = check_inputs( inputs )
        
        if ~ispc
            classifierFile = '/shared/Jordan/swd_detection/models/norm-v-swd.mat';
        else
            classifierFile = 'Z:\Jordan\Seizure_detection\swd_detection\models\norm-v-swd.mat';
        end
        
        names = {'pThresh','minLength','maxGap','classifier','W','refine'};
        defaults = {nan,1,1,classifierFile,[],true};
        
        options = inputParser;
        for j = 1:numel( names )
            options.addParameter( names{j},defaults{j} );
        end
        
        parse( options,inputs{:} );
        options = options.Results;
        
        assert( isfile( options.classifier ),'supplied classifier file not found!' );
        assert( options.minLength > 0,'minimum length must be positive' );
        
        model = load( options.classifier,'model' );
        model = model.model;
        if isnan( options.pThresh )
            options.pThresh = median( model.p1 ) - 2*mad( model.p1 );
        else
            assert( (options.pThresh >= 0 & options.pThresh <= 1),'probability threshold must be between 0 and 1!' );
        end
    end

    function bound = search_probability_space( x,start )
        % performs a forward and backward search of probabilities in "p" 
        % starting from "start", and connects regions with high probability
        
        % forward search
        stop = start;
        while stop < length( x )
            newStop = find( swdBool(stop+1:end),1 ) + stop;
            if x(newStop,1) - x(stop,2) <= options.maxGap
                stop = newStop;
            else
                break
            end
        end
        
        % backward search
        while start > 1
            newStart = find( swdBool(1:start-1),1,'last' );
            if x(start,1) - x(newStart,2) <= options.maxGap
                start = newStart;
            else
                break
            end
        end
        
        bound = [x(start,1),x(stop,2)];
    end
end
    