function [accuracy,ppv,npv,F] = score_classification( C )
    % [accuracy,ppv,npv,F] = score_classification( C )
    %
    % computes classification scores from a given confusion matrix
    %
    % inputs:
    %   C - confusion matrix (see "confusion_matrix.m")
    %   
    % outputs:
    %   accuracy - (TP + TN) / (TP + FP + TN + FN)
    %   ppv - TP / (TP + FP)
    %   npv - TN / (TN + FN)
    %   F - 2 * (TP * ppv)/(TP + ppv)
    % 
    % By Jordan Sorokin, 2/2/19
    
    accuracy = trace( C ) / sum( C(:) );
    ppv = C(1,1) / sum( C(:,1) );
    npv = C(2,2) / sum( C(:,2) );
    F = 2 * (C(1,1)*ppv) / (C(1,1) + ppv);
end
    