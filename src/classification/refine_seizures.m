function S = refine_seizures( S,tolerance )
    % S = refine_seizures( S,tolerance )
    %
    % refines seizures by iteratively removing small chunks of 
    % data in an outward -> in fashion, each time checking if the
    % there are at least 3 sequential points with p(swd) > threshold
    % from the current starting point. If not, the current starting point
    % is discarded, and the next time point is considered. 
    %
    % Inputs:
    %   S - the structure obtained via "detect_seizures.m"
    %
    %   tolerance - a positive real number. This defines the probability
    %               threshold as the value where at least 4 consecutive
    %               data segments have a probability greater than this value
    %               to be considered the start of the seizure
    % Output:
    %   S - a refined seizure structure
    %
    % By Jordan Sorokin, 2/8/19
    
    removeSzr = false( S.nEvents,1 );
    nAdvance = 4;
    for szr = 1:S.nEvents
        P = S.P{szr};
        nSegs = length( P );
        
       % smallProb = P < (median( P ) - tolerance*mad( P ));
        smallProb = P < tolerance;
        if all( smallProb )
            removeSzr(szr) = true;
        else
            badSegs = diff( smallProb );
            
            if smallProb(1) == 0
                start = 1;
            else
                start = find( badSegs==-1,1 );
            end

            % iterate over starting locations
            while start <= length( smallProb ) - 4
                if nnz( smallProb(start:start+nAdvance) ) == 0
                    break
                else
                    start = start + 1;
                end
            end
            
            if start > length( smallProb )-4
                removeSzr(szr) = true;
                continue
            end
            
            % find new ending location
            if smallProb(end) == 0
                stop = nSegs;
            else
                stop = find( badSegs,1,'last' );
            end
        end
        
        % check if we should remove this seizure entirely
        if ~removeSzr(szr)
            chunkLength = (S.end(szr) - S.start(szr)) / length(P);
            S.start(szr) = S.start(szr) + chunkLength*(start-1);
            S.end(szr) = S.end(szr) - chunkLength*(nSegs-stop);
            S.P{szr} = P(start:stop);
            S.X{szr} = S.X{szr}(start:stop,:);

            idx = ceil(chunkLength*S.fs*(start-1)+1):floor(chunkLength*S.fs*stop);
            idx(end) = min( idx(end),length( S.V{szr} ) );
            S.V{szr} = S.V{szr}(idx);
            S.W{szr} = S.W{szr}(:,idx);
        end
    end
        
    % now remove any seizures that are too small
    removeSzr =  removeSzr | (S.end - S.start) < S.params.minLength;
    S.start(removeSzr) = [];
    S.end(removeSzr) = [];
    S.P(removeSzr) = [];
    S.X(removeSzr) = [];
    S.V(removeSzr) = [];
    S.W(removeSzr) = [];
    S.dKL(removeSzr,:) = [];
    S.nEvents = nnz( ~removeSzr );
    S.pThresh_refined = tolerance;
end
    