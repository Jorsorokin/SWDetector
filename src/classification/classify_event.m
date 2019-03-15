function [P,Yhat] = classify_event( X,coeffs,thresh )
    % [P,Yhat] = classify_event( X,coeffs,thresh )
    %
    % classify an event using the logistic function based on learned
    % coefficients. Thresh is between [0,1] representing the probability
    % threshold to count an event as group 0 (Yhat = 0) or group 1 (Yhat =
    % 1). So, thresh = 0.7 means we must be 70% sure to conclude Yhat = 1. 
    %       - i.e. P(Yhat = 1) = 0.3 on the average
    
    logodds = X*coeffs(2:end) + coeffs(1);
    P = 1 ./ (1 + exp( -logodds ));
    Yhat = P >= thresh; 
end