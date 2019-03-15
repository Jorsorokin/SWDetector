function ax = plot_detection_results( data,fs,S,W )
    % ax = plot_detection_results( data,fs,S,(W) )
    %
    % plots a channel of eeg + detected seizures on top
    % Optionally plots the wavelet coefficients if they are 
    % provided, or if W = true (compute wavelet and plot)
    
    if nargin < 4
        W = modwt( data,'sym4',8 );
    end
    
    t = (1:length(data))/fs;
    kernel = gausswin( 25,2 ); 
    kernel = kernel / sum( kernel );

    figure;    
    ax(1) = subplot( 2,1,1 );
    imagesc( t,1:size(W,1),conv2( abs(W)',kernel,'same' )' );
    
    ax(2) = subplot( 2,1,2 ); hold on;
    plot( t,data ); xlim( [1 t(end)] );
    for i = 1:length( S.start )
        plot( [S.start(i),S.end(i)], [600 600],'linewidth',3,'color',repmat( 1-mean(S.P{i}),1,3 ) );
    end
    
    linkaxes( ax,'x' );
end
    
    
    
            
           