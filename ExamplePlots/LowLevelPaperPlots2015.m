%%%%%%%%%%Low level transcription plot

axLow = Data23.plotNucVsCountScatter('cy5','','Nuclei count vs. transcripts',ColorElt2);
hold(axLow,'on');
Data26.overlayNucVsCountScatter(axLow,'cy5','', ColorElt2);
Data26.overlayNucVsCountScatter(axLow,'tmr','',Color26);
hold(axLow,'on');
Data23.overlayNucVsCountScatter(axLow,'tmr','',Color23);
hold(axLow,'on');
axis(axLow,[0 160 0 80])

lowLeg = legend(axLow,'elt-2','1879:1A gfp','613:1A gfp');
%set(lowLeg,'Position',[0.75,0.75,0.1,0.1],'FontSize',10);

figLow = get(axLow,'Parent');
set(figLow,'Position',[0 0 500,300]);

set(lowLeg,'FontSize',10,'Position',[0.75 0.75 0.1 0.1]);
%%Add E-cell annotations
DuPlotWormData.addECellAnnotationsNucs(axLow,1);


%%%%%%%%%%Alternate Low level transcription plot

axLow2 = Data05.plotNucVsCountScatter('tmr','gfp','Nuclei count vs. transcripts',Color05);
hold(axLow2,'on');
Data18.overlayNucVsCountScatter(axLow2,'tmr','',Color18);
Data26.overlayNucVsCountScatter(axLow2,'tmr','',Color26);
Data23.overlayNucVsCountScatter(axLow2,'tmr','',Color23);
axis(axLow2,[0 160 0 80])

lowLeg2 = legend(axLow2,'1879:WT gfp','613:WT gfp','1879:1A gfp','613:1A gfp');
%set(lowLeg,'Position',[0.75,0.75,0.1,0.1],'FontSize',10);

figLow2 = get(axLow2,'Parent');
DuPlotWormData.addECellAnnotationsNucs(axLow2,1);


set(figLow2,'Position',[0 0 500,300]);
set(lowLeg2,'FontSize',10,'Position',[0.75 0.75 0.1 0.1]);

%%%%%%%%%%%%%%%%%%%Make markers clear for low level transcription plots

set(findobj(axLow,'Type','hggroup'),'MarkerFaceColor','none');
set(findobj(axLow2,'Type','hggroup'),'MarkerFaceColor','none');
