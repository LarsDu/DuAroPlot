figAx18cc = figure('Position',[0,0,500,500]);

ax18cc = Data18.plotCountsVsCounts(Data18.CY5_IND,Data18.TMR_IND,'r');
set(ax18cc,'Xlim',[0 1200],'Ylim',[0 1200]);

ylabel(ax18cc,['gfp transcripts']);
xlabel(ax18cc,['elt-2 transcripts']);
titleh = title('Hello')

set(ax18cc,'Title',titleh)