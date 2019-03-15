function [X,Y] = create_wavelet_features( R,varargin )
    % [X,Y] = create_wavelet_features( R,(sqrtTransform,useInteractions) )
    %
    % creates features for classifier from the R-struct
    %   - currently features are: 
    %           x = v_i / sum_j{ v_j }
    %     
    %     where v = var( W_i ), the variance of the ith wavelet coefficient 
    %     over the duration of the data chunk
    %
    %   - if "sqrtTransform" is true, features are: 
    %           x = sqrt( v_i ) / sum_j{ sqrt( v_j )}
    %
    %   - if "useInteractions" is true, the relative wavelet variances are
    %     expanded using a quadratic basis:
    %
    %           i.e. [x1,x2] --> [x1,x2,x1^2,x2^2,x1*x2]
    %    
    %     this results in possibly non-linear interactions between wavelet
    %     variances that may improve classification; however, the
    %     dimensionality of the dataset becomes much larger, and thus we
    %     need quadratically-more datapoints to learn the features...
    % 
    % Because different animals/files in the multi-dimensional R-struct 
    % (created via "mark_swd.m") may have different numbers of seizures, 
    % this function extracts the same number of samples from each file
    % to equally represent the features among different animals/files
    
    % check if we want to log-transform the wavelet variances
    if nargin > 1 && ~isempty( varargin{1} )
        useLog = varargin{1};
    else
        useLog = false;
    end
    
    if nargin > 2 && ~isempty( varargin{2} )
        useInteractions = varargin{2};
    else
        useInteractions = false;
    end
    
    nSWD = [R.nSWD];
    nArtifact = [R.nArtifact];
    nNormal = [R.nNormal]; 
    nWave = size( R(1).data.W,1 ); % number of wavelet bases 
    if useInteractions
        nWave = nWave + nWave^2;
    end
    
    % extract the minimum # of SWD from the animals, and the same number of
    % normal data points as that of SWD, and min( nArtifact,nSWD ) number
    % of artifactual data points
    nSWD_max = min( nSWD );
    nArtifact = arrayfun( @(x)(min(x,nSWD_max)),nArtifact ); % no more than n per animal, but variable between animals
    nAnimals = numel( R );
    N = (nSWD_max*2)*nAnimals + sum(nArtifact);
    
    X = nan( N,nWave );
    Y = nan( N,1 );
    count = 0;
    
    for i = 1:nAnimals
        if useLog
            w = squeeze( sqrt( var( R(i).data.W,[],2 ) ) );
        else
            w = squeeze( var( R(i).data.W,[],2 ) );
        end
        
        if useInteractions
            w = quadratic_expansion( w' )';
        end
        
        x = w ./ sum( w,1 );
        x = w;
        y = R(i).data.label;
        
        % extract random samples for the different labels
        for class = unique( y )
            label = class{:};
            switch label
                case 'normal'
                    inds = randperm( nNormal(i),nSWD_max );
                    classVal = 0;
                    n = nSWD_max;
                case 'swd'
                    inds = randperm( nSWD(i),nSWD_max );
                    classVal = 1;
                    n = nSWD_max;
                case 'artifact'
                    inds = 1:nArtifact(i);
                    classVal = -1;
                    n = nArtifact(i);
            end
            
            idx = contains( y,label );
            xsub = x(:,idx);
            
            X(count+1:count+n,:) = xsub(:,inds)';
            Y(count+1:count+n) = classVal;
            count = count + n;
        end
    end
end 