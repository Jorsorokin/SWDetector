%% HYPERPARAMS
ALPHA = 0.05:0.05:1; % balance between ridge (0) and lasso (1) for elastic net
nALPHA = length( ALPHA );
modelType = 'logistic'; 
regularization = 'elastic net';
eventType = "stereotyped absence seizures";                             % can change this as desired

% if you used "get_swd_events_from_files", leave this as is
comparisons = {[0,1]};%, [0,-1], [-1,1]};
names = {'norm-v-swd_noRelativeVar'};%, 'norm-v-artifact', 'artifact-v-swd'};

% change your input/output folders here:
inputFolder = 'Z:\Jordan\Seizure_detection\swd_detection\data\';        % folder containing results from "get_swd_events_from_files"
inputFileName = 'data_segments_swd_noRelativeVar.mat';                             % contains your training data and parameters
modelOutFolder = 'Z:\Jordan\Seizure_detection\swd_detection\models';    % saves the model hyper parameters and coefficients
figOutFolder = 'Z:\Jordan\Seizure_detection\swd_detection\figs';        % saves some figures with parameter search info

%% partition the data once, so we can compare the same training sets
load( [inputFolder,filesep,inputFileName],'X','Y','params' ); % contains features, labels, and parameters
partitions = partition_dataset( size( X,1 ),[0.75 0.2 0.05],X,Y );
M = size( X,2 );

%% loop over comparisons
for comp = 1:length( comparisons )
    
    % create training and val sets for this set of groups
    IDX = ismember( partitions.train.Y,comparisons{comp} );
    trainX = partitions.train.X(IDX,:);
    trainY = partitions.train.Y(IDX,:);
    
    IDX = ismember( partitions.val.Y,comparisons{comp} );
    valX = partitions.val.X(IDX,:);
    valY = partitions.val.Y(IDX,:);
       
    if comp == 2
        trainY(trainY == -1) = 1;
        valY(valY == -1) = 1;
    elseif comp == 3
        trainY(trainY == -1) = 0;
        valY(valY == -1) = 0;
    end
    
    % create vecs for storing results
    ACCURACY = zeros( 1,nALPHA );
    PPV = zeros( 1,nALPHA );
    NPV = zeros( 1,nALPHA );
    DF = zeros( 1,nALPHA );
    B = zeros( M,nALPHA );
    B0 = zeros( 1,nALPHA ); 
    LAMBDA = zeros( 1,nALPHA );
    DEV = zeros( 1,nALPHA );
    count = 1;
    
    % loop over alpha to balance lasso with ridge
    for alpha = ALPHA
        
        % train an elastic net
        [b,stats] = lassoglm( trainX,trainY,'binomial','Alpha',alpha,'CV',10 );
        
        % predict validation set response
        [~,yhat] = classify_event( valX,[stats.Intercept(stats.Index1SE); b(:,stats.Index1SE)],0.5 );
        
        % store the results
        C = confusion_matrix( valY,yhat );
        [ACCURACY(count),PPV(count),NPV(count),~] = score_classification( C );
        DF(count) = stats.DF(stats.Index1SE);
        B(:,count) = b(:,stats.Index1SE);
        B0(count) = stats.Intercept(:,stats.Index1SE);
        LAMBDA(count) = stats.Lambda1SE;
        DEV(count) = stats.Deviance(stats.Index1SE);
        count = count + 1;
    end
    
    % store the hyperparams for this comparison
    [~,bestInd] = min( DEV );
    hyperparams = struct( 'alpha',ALPHA,'accuracy',ACCURACY,'ppv',PPV,'npv',NPV,...
                          'DF',DF,'B',B,'B0',B0,'lambda',LAMBDA,...
                          'deviance',DEV,'optimalIDX',bestInd,...
                          'trainingFile',[inputFolder,filesep,inputFileName] );
                      
    % now plot the results
    figure( 'Name',names{comp} );
    subplot( 3,2,1 ); hold on;
    plot( ALPHA,ACCURACY );
    plot( ALPHA(bestInd),ACCURACY(bestInd),'ro' );
    ylabel( 'accuracy' );
    
    subplot( 3,2,2 ); hold on;
    plot( ALPHA,DEV );
    plot( ALPHA(bestInd),DEV(bestInd),'ro' );
    ylabel( 'deviance' );  
    
    subplot( 3,2,3 ); hold on;
    plot( ALPHA,PPV,ALPHA,NPV );
    plot( ALPHA(bestInd),PPV(bestInd),'ko',ALPHA(bestInd),NPV(bestInd),'ko' );
    legend( {'PPV','NPV'},'box','off' );
    
    subplot( 3,2,4 ); hold on;
    plot( ALPHA,LAMBDA );
    plot( ALPHA(bestInd),LAMBDA(bestInd),'ro' );
    ylabel( 'lambda' );
    
    subplot( 3,1,3 ); hold on
    errorbar( mean(B,2),std(B,[],2),'k.' );
    plot( B(:,bestInd),'ro' );
    ylabel( 'coeff value' );

    % save the results
    set( gcf,'Renderer','painters' );
    saveas( gcf, [figOutFolder,filesep,names{comp},'_hyperparams'],'fig' );    

    % finally, save the hyperparameters, and a copy of the parameters
    % plus optimal coefficients into a variable "model" 
    % so that we can quickly load in the model for classification
    model = params;
    model.type = modelType;
    model.regularization = regularization;
    model.eventType = eventType;
    model.B0 = hyperparams.B0(bestInd);
    model.B = hyperparams.B(:,bestInd);
    model.nCoeffs = size( model.B,1 );
    model.df = hyperparams.DF(bestInd);
    
    % compute probabilities and predicted response from best model and
    % store the probs into the "model" param
    [P,~] = classify_event( valX,[B0(bestInd); B(:,bestInd)],0.5 );
    model.p0 = P(valY == 0);
    model.p1 = P(valY == 1);
    
    save( [modelOutFolder,filesep,names{comp},'.mat'],'hyperparams','model','inputFileName' );
end
    
    