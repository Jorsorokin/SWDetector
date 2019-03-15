function plot_each_detected_seizure( S )
    % plot_each_detected_seizure( S )
    % 
    % plots each seizure along with its wavelet representation and
    % probability vector and allows user to cycle through each by pressing
    % the space bar to advance to the next seizure
    
    ax(1) = subplot( 3,7,1:6 ); hold on
    ax(2) = subplot( 3,7,8:13 ); hold on
    ax(3) = subplot( 3,7,15:20 ); hold on
    linkaxes( ax,'x' );
    
    sideax(1) = subplot( 3,7,7 ); hold on;
    sideax(2) = subplot( 3,7,14 ); hold on;
    sideax(3) = subplot( 3,7,21 ); hold on;
    
    ax1BinCenters = 1:size( S.W{1},1 );
    ax2Bins = linspace( min( cellfun(@min,S.V) ),max( cellfun(@max,S.V) ),21 );
    ax2BinCenters = mean([ax2Bins(2:end);ax2Bins(1:end-1)]);
    ax3Bins = 0:0.1:1;
    ax3BinCenters = mean([ax3Bins(2:end);ax3Bins(1:end-1)]);    
    
    ynames = {'wavelet scale','mV','P(szr)'};
    ylims = {ax1BinCenters([1,end]),ax2BinCenters([1,end]),ax3BinCenters([1,end])};
    for i = 1:3
        set( ax(i),'tickdir','out','box','off','ylim',ylims{i} );
        ylabel( ax(i),ynames{i} );
        
        set( sideax(i),'tickdir','out','box','off',...
            'ylim',ylims{i},'xlim',[0 1],'ytick',[],'xtick',[0 0.5 1] );
    end

    for i = 1:S.nEvents
        arrayfun( @cla,ax ); arrayfun( @cla,sideax );
        
        tSzr = linspace( S.start(i),S.end(i),length(S.V{i}) );
        tProb = linspace( S.start(i),S.end(i),length(S.P{i}) );
        
        imagesc( ax(1),tSzr,ax1BinCenters,flipud(S.W{i}) ); 
        plot( ax(2),tSzr,S.V{i},'k' );
        plot( ax(3),tProb,S.P{i},'o-' );
        set( ax(3),'xlim',[S.start(i),S.end(i)],'ylim',[0,1] );
        
        % side plots
        p = var( S.W{i},[],2 ) ./ sum( var( S.W{i},[],2 ) );
        barh( sideax(1),ax1BinCenters,p,'FaceColor','k','EdgeColor','k' );
        
        p = histcounts( S.V{i},ax2Bins,'normalization','probability' );
        barh( sideax(2),ax2BinCenters,p,'FaceColor','k','EdgeColor','k' );
        
        p = histcounts( S.P{i},ax3Bins,'normalization','probability' );
        barh( sideax(3),ax3BinCenters,p,'FaceColor','k','EdgeColor','k' );
        
        pause;
    end 
end