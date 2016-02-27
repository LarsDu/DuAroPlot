%Color spec for different strains
ColorElt2 = [.9 0 0]; %'dark red'
Color05 = 'magenta';
Color18 = 'green';
Color23 = 'cyan';
Color26 = 'blue';
Color06 = 'black';
Color12 = [1 0.5 0]; %orange
Color07 = [0.5 0.17 0.85]; %purple
Color28 = [.6 .6 .6]; %sixty shades of gray
Color25 = [0, .39,0]; %dark green


%%%%%%%%%%Aggregate elt-2 data versus 1879:WT
fig2A = figure('Position',[0,0,500,400])


hax0 = DataElt2.plotNucVsCountScatterSpline('cy5','','','red');
lineh0= findobj(hax0,'Type','Line');
hax1 = Data05.plotNucVsCountScatterSpline('tmr','','','blue');
lineh1= findobj(hax1,'Type','Line');


%Transfer plot data to ax0

DuPlotWormData.addPlotToAxes(hax1,hax0);





%Set legend
splineLegend = legend([lineh0,lineh1],...
               'elt-2',...
               '1879:WT gfp');
title(hax0,['Mean Expression Trajectory'],'FontSize',14,...
    'FontWeight','bold')

set(splineLegend,'Position',[0.25,0.68,0.1,0.1])

set(hax0,'Parent',fig2A);
set(splineLegend,'Parent',fig2A);
%%Add E-cell annotations
DuPlotWormData.addECellAnnotationsNucs(hax0);

close (ancestor(hax1,'figure'));


%%%%%%%%%%Aggregate elt-2 data versus 1879:WT
fig2Ab = figure('Position',[0,0,500,400])


hhax0 = DataElt2.plotNucVsCountScatter('cy5','','','red');
%lineh0= findobj(hax0,'Type','Line');
hhax1 = Data05.plotNucVsCountScatter('tmr','','','blue');
%lineh1= findobj(hax1,'Type','Line');


%Transfer plot data to ax0

DuPlotWormData.addPlotToAxes(hhax1,hhax0);





%Set legend
splineLegend2 = legend('elt-2',...
               '1879:WT gfp');
title(hhax0,['Mean Expression Trajectory'],'FontSize',14,...
    'FontWeight','bold')

set(splineLegend2,'Position',[0.25,0.68,0.1,0.1])

set(hhax0,'Parent',fig2Ab);
set(splineLegend2,'Parent',fig2Ab);
%%Add E-cell annotations
DuPlotWormData.addECellAnnotationsNucs(hhax0);

close (ancestor(hhax1,'figure'));

