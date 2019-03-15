function partitions = partition_dataset( N,P,X,Y )
    % partitions = partition_dataset( N,P,(X,Y))
    %
    % partitions a dataset into training, validation, and test
    % sets based on user-specified percentages of total N
    %
    % inputs:
    %   N - the total # of samples in the dataset
    %
    %   P - a 3-element vector with each element 0 < x < 1, 
    %       specifying the percent of data to partition into 
    %       the training (1st pos), validation (2nd), and test (3rd)
    %
    %   (X) - optional N x D matrix
    %
    %   (Y) - optional N x 1 class labels
    %
    % outputs:
    %   partitions - structure containing the INDICIES of the training,
    %                validation, and test sets. If the matrix X is
    %                provided, then included as well is the actual
    %                partition of the dataset (for convenience)
    %
    % By Jordan Sorokin, 2/1/19
    
    assert( sum( P ) == 1,'total percent of partitions must equal 1' );
    nTrain = floor( N*P(1) );
    nVal = floor( N*P(2) );    
    allInds = 1:N;

    partitions.train.inds = randperm( N,nTrain );
    allInds(partitions.train.inds) = [];
    
    partitions.val.inds = allInds(randperm(N-nTrain,nVal));
    allInds(ismember( allInds,partitions.val.inds )) = [];
    
    partitions.test.inds = allInds;
    
    if nargin > 2 
        assert( size(X,1)==N,'X must contain N samples as well' );
        
        partitions.train.X = X(partitions.train.inds,:);
        partitions.val.X = X(partitions.val.inds,:);
        partitions.test.X = X(partitions.test.inds,:);
    end
    
    if nargin > 3
        assert( size(Y,1)==N,'Y must contain N samples as well' );
        
        partitions.train.Y = Y(partitions.train.inds,:);
        partitions.val.Y = Y(partitions.val.inds,:);
        partitions.test.Y = Y(partitions.test.inds,:);
    end
end