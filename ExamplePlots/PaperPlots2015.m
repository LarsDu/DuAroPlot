%%Note: If you update your classes, type 'clear classes' and reload
%%class objects into the workspace. They should update with no problem


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Nuc vs. count plot for most strains
axLenVar = Data18_failure_purge.plotNucVsCountScatter('cy5','',['Nuclei ' ...
                    'count vs.mRNA'],ColorElt2)
figLenVar = ancestor(axLenVar,'Figure');
Data05.overlayNucVsCountScatter(axLenVar,'tmr','',Color05);
Data18_failure_purge.overlayNucVsCountScatter(axLenVar,'tmr','',Color18);
Data23.overlayNucVsCountScatter(axLenVar,'tmr','',Color23);
Data26.overlayNucVsCountScatter(axLenVar,'tmr','',Color26);

%axis(axesGfpMut,[0 160 0 1000]);
set(figLenVar,'Position',[0,0,640,480])

%Set legend
s1leg = legend(axLenVar,...
               '613:WT elt-2',...
               '1879:WT gfp',...   
               '613:WT gfp',... 
               '613:A gfp',...
               '1879:A gfp');			 			   
set(s1leg,'Position',[0.25,0.68,0.01,0.01],'FontSize',10);

DuPlotWormData.addECellAnnotationsNucs(axLenVar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

%%Plot Mean Trajectory using Savitsky Golay at 50% data span
axesMeanTraj = axes
Data18_failure_purge.plotMeanTrajSavitskyGolay(axesMeanTraj,9,5,ColorElt2);
hold(axesMeanTraj,'on')
Data05.plotMeanTrajSavitskyGolay(axesMeanTraj,9,6,Color05)
hold(axesMeanTraj,'on');

Data18_failure_purge.plotMeanTrajSavitskyGolay(axesMeanTraj,9,6,Color18);
hold(axesMeanTraj,'on')

%Data06.plotMeanTrajSavitskyGolay(axesMeanTraj,9,6,'k')
%hold(axesMeanTraj,'on')
%Data07.plotMeanTrajSavitskyGolay(axesMeanTraj,9,6,'m')
%hold(axesMeanTraj,'on')
%Data12.plotMeanTrajSavitskyGolay(axesMeanTraj,9,6,'r')
%hold(axesMeanTraj,'on')
Data23.plotMeanTrajSavitskyGolay(axesMeanTraj,9,6,Color23);
hold(axesMeanTraj,'on')
Data26.plotMeanTrajSavitskyGolay(axesMeanTraj,9,6,Color26);
hold(axesMeanTraj,'on')
axis(axesMeanTraj,[0 175 0 800])

%Set linewidth on axesMeanTraj plot
set(axesMeanTraj, 'LineWidth', 2);
a=findobj(axesMeanTraj);
alllines=findall(a,'Type','line');
set(alllines, 'lineWidth', 2);
title(axesMeanTraj,['Mean Trajectory over time at 50% data span'])
xlabel(axesMeanTraj,['Time (min interpolated from # nuclei)'])
ylabel(axesMeanTraj,['# gfp transcripts'])

%Set legend
s2leg = legend(axesMeanTraj,...
               '613:WT elt-2',...
               '1879:WT gfp',...
               '613:WT gfp',... 
               '613:A gfp', ...
               '1879:A gfp');
               
%Set fontsize
set(axesMeanTraj,'fontsize', 12);
set(s2leg,'Position',[0.25,0.68,0.1,0.1]);

DuPlotWormData.addECellAnnotationsNucs(axesMeanTraj);

%%%%%%%%%%Plot mean trajectories using splines
%First plot separately
figSpline = figure('Position',[0,0,500,400])


ax0 = Data18_failure_purge.plotNucVsCountSpline('cy5','','',ColorElt2);
ax1 = Data05.plotNucVsCountSpline('tmr','','',Color05);
ax2 = Data18_failure_purge.plotNucVsCountSpline('tmr','','',Color18);
ax3 = Data23.plotNucVsCountSpline('tmr','','',Color23);
ax4 = Data26.plotNucVsCountSpline('tmr','','',Color26);

%Transfer plot data to ax0

DuPlotWormData.addPlotToAxes(ax1,ax0);
DuPlotWormData.addPlotToAxes(ax2,ax0);
DuPlotWormData.addPlotToAxes(ax3,ax0);
DuPlotWormData.addPlotToAxes(ax4,ax0);



%Set legend
splineLegend = legend(ax0,...
               '613:WT elt-2',...
               '1879:WT gfp',...
               '613:WT gfp',... 
               '613:A gfp', ...
               '1879:A gfp');
title(ax0,['Mean Expression Trajectory'],'FontSize',14,...
    'FontWeight','bold')
%set(splineLegend,'Location','Northwest')
set(splineLegend,'Position',[0.25,0.68,0.1,0.1])

set(ax0,'Parent',figSpline);
set(splineLegend,'Parent',figSpline);
%%Add E-cell annotations
DuPlotWormData.addECellAnnotationsNucs(ax0);


close (ancestor(ax1,'figure'));
close (ancestor(ax2,'figure'));
close (ancestor(ax3,'figure'));
close (ancestor(ax4,'figure'));

%%%%%%%%%%Low level transcription plot

axLow = Data18_failure_purge.plotNucVsCountScatter('cy5','gfp','Nuclei count vs. transcripts',ColorElt2);
hold(axLow,'on');
Data26.overlayNucVsCountScatter(axLow,'tmr','',Color26);
hold(axLow,'on');
Data23.overlayNucVsCountScatter(axLow,'tmr','',Color23);
hold(axLow,'on');
axis(axLow,[0 160 0 80])

lowLeg = legend(axLow,'613:WT elt-2','1879:1A gfp','613:1A gfp');
%set(lowLeg,'Position',[0.75,0.75,0.1,0.1],'FontSize',10);

figLow = get(axLow,'Parent');
set(figLow,'Position',[0 0 500,300]);

set(lowLeg,'FontSize',10,'Position',[0.75 0.75 0.1 0.1]);
%%Add E-cell annotations
DuPlotWormData.addECellAnnotationsNucs(axLow,1);


%%%%%%%%%%Alternate Low level transcription plot

axLow2 = Data05.plotNucVsCountScatter('tmr','gfp','Nuclei count vs. transcripts',Color05);
hold(axLow2,'on');
Data12.overlayNucVsCountScatter(axLow2,'tmr','',Color12);
Data26.overlayNucVsCountScatter(axLow2,'tmr','',Color26);
Data23.overlayNucVsCountScatter(axLow2,'tmr','',Color23);
axis(axLow2,[0 160 0 80])

lowLeg2 = legend(axLow2,'1879:WT gfp','613:WT gfp1212','1879:1A gfp','613:1A gfp');
%set(lowLeg,'Position',[0.75,0.75,0.1,0.1],'FontSize',10);

figLow2 = get(axLow2,'Parent');
set(figLow2,'Position',[0 0 500,300]);
set(lowLeg2,'FontSize',10,'Position',[0.75 0.75 0.1 0.1]);
%%Add E-cell annotations
DuPlotWormData.addECellAnnotationsNucs(axLow2,1);

%%%%%%%%%%%%%%%%%%%Make markers clear for low level transcription plots

set(findobj(axLow,'Type','hggroup'),'MarkerFaceColor','none');
set(findobj(axLow2,'Type','hggroup'),'MarkerFaceColor','none');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nuc vs. count plot for secondary HGATAR mutants

axSecondary = Data05.plotNucVsCountScatter('tmr', 'gfp', ['Nuclei ' ...
                    'count vs.mRNA'],Color05)
figSecondary = ancestor(axSecondary,'Figure');
Data06.overlayNucVsCountScatter(axSecondary,'tmr','',Color06);
Data12.overlayNucVsCountScatter(axSecondary,'tmr','',Color12);


set(figSecondary,'Position',[0,0,640,480])

%Set legend
secondaryLeg = legend(axSecondary,...
               '1879:WT gfp',...
               '1879:4G gfp',...   
               '1879:11G gfp');			 			   
set(secondaryLeg,'Position',[0.25,0.68,0.05,0.05],'FontSize',10);

DuPlotWormData.addECellAnnotationsNucs(axSecondary);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nuc vs. count plot for secondary HGATAR mutants part 2
figure;
DuPlotWormData.plotCloudsAndMeans(Data18_failure_purge,'tmr',Data25,'tmr',Color18,Color25, 0.0001,'gfp','gfp');


%%%%%%%%%%Plot mean trajectories using splines for secondary HGATAR mutants
%First plot separately
figSpline2 = figure('Position',[0,0,500,400])


sax0 = Data05.plotNucVsCountSpline('tmr','','',Color05);
sax1 = Data18_failure_purge.plotNucVsCountSpline('tmr','','',Color18);
sax2 = Data06.plotNucVsCountSpline('tmr','','',Color06);
sax3 = Data12.plotNucVsCountSpline('tmr','','',Color12);
sax4 = Data25.plotNucVsCountSpline('tmr','','',Color25);

%Transfer plot data to ax0

DuPlotWormData.addPlotToAxes(sax1,sax0);
DuPlotWormData.addPlotToAxes(sax2,sax0);
DuPlotWormData.addPlotToAxes(sax3,sax0);
DuPlotWormData.addPlotToAxes(sax4,sax0);



%Set legend
splineLegend2 = legend(sax0,...
               '1879:WT gfp',...
               '613:WT gfp',... 
               '1879:4G gfp', ...
               '1879:11G gfp',...
               '613:4G gfp');
title(sax0,['Mean Expression Trajectory'],'FontSize',14,...
    'FontWeight','bold')

set(splineLegend2,'Position',[0.25,0.68,0.1,0.1])

set(sax0,'Parent',figSpline2);
set(splineLegend2,'Parent',figSpline2);
%%Add E-cell annotations
DuPlotWormData.addECellAnnotationsNucs(sax0);


close (ancestor(sax1,'figure'));
close (ancestor(sax2,'figure'));
close (ancestor(sax3,'figure'));
close (ancestor(sax4,'figure'));

%%%%%%%%%%Secondary HGATAR mutants with splines and points
%First plot separately
figSpline3 = figure('Position',[0,0,500,400])


ssax0 = Data05.plotNucVsCountScatterSpline('tmr','','',Color05);
lineh0= findobj(ssax0,'Type','Line');
ssax1 = Data06.plotNucVsCountScatterSpline('tmr','','',Color06);
lineh1= findobj(ssax1,'Type','Line');
ssax2 = Data12.plotNucVsCountScatterSpline('tmr','','',Color12);
lineh2= findobj(ssax2,'Type','Line');

%Transfer plot data to ax0

DuPlotWormData.addPlotToAxes(ssax1,ssax0);
DuPlotWormData.addPlotToAxes(ssax2,ssax0);




%Set legend
splineLegend3 = legend([lineh0,lineh1,lineh2],...
               '1879:WT gfp',...
               '1879:4G gfp', ...
               '1879:11G gfp');
title(ssax0,['Mean Expression Trajectory'],'FontSize',14,...
    'FontWeight','bold')

set(splineLegend3,'Position',[0.25,0.68,0.1,0.1])

set(ssax0,'Parent',figSpline3);
set(splineLegend3,'Parent',figSpline3);
%%Add E-cell annotations
DuPlotWormData.addECellAnnotationsNucs(ssax0);

close (ancestor(ssax1,'figure'));
close (ancestor(ssax2,'figure'));



%%%Plot coeffVar
newAx = axes;
Data05.plotCoeffVar(newAx,Data05.NUCLEI_IND,Data05.TMR_IND,Color05);
Data18_failure_purge.plotCoeffVar(newAx,Data18_failure_purge.NUCLEI_IND,Data18_failure_purge.TMR_IND,Color18);
Data06.plotCoeffVar(newAx,Data06.NUCLEI_IND,Data06.TMR_IND,Color06);
Data12.plotCoeffVar(newAx,Data12.NUCLEI_IND,Data12.TMR_IND,Color12);
Data25.plotCoeffVar(newAx,Data25.NUCLEI_IND,Data25.TMR_IND,Color25);
%Data26.plotCoeffVar(newAx,Data26.NUCLEI_IND,Data26.TMR_IND,Color26);
%Data23.plotCoeffVar(newAx,Data23.NUCLEI_IND,Data23.TMR_IND,Color23);
%Data07.plotCoeffVar(newAx,Data07.NUCLEI_IND,Data07.TMR_IND,Color07);
%Data28.plotCoeffVar(newAx,Data28.NUCLEI_IND,Data28.TMR_IND,Color28);

%%%%%%%%%%Length Variants with splines and points
%First plot separately
figSpline4 = figure('Position',[0,0,500,400])


axh0 = Data05.plotNucVsCountScatterSpline('tmr','','',Color05);
lineaxh0= findobj(axh0,'Type','Line');
axh1 = Data18_failure_purge.plotNucVsCountScatterSpline('tmr','','',Color18);
lineaxh1= findobj(axh1,'Type','Line');
axh2 = Data26.plotNucVsCountScatterSpline('tmr','','',Color26);
lineaxh2= findobj(axh2,'Type','Line');
axh3 = Data23.plotNucVsCountScatterSpline('tmr','','',Color23);
lineaxh3= findobj(axh3,'Type','Line');

%Transfer plot data to ax0

DuPlotWormData.addPlotToAxes(axh1,axh0);
DuPlotWormData.addPlotToAxes(axh2,axh0);
DuPlotWormData.addPlotToAxes(axh3,axh0);


%Set legend
splineLegend4 = legend([lineaxh0,lineaxh1,lineaxh2,lineaxh3],...
               '1879:WT gfp',...
               '613:WT gfp',...
               '1879:A gfp',...
               '613:A gfp');
title(axh0,['Mean Expression Trajectory'],'FontSize',14,...
    'FontWeight','bold')

set(splineLegend4,'Position',[0.25,0.68,0.1,0.1])

set(axh0,'Parent',figSpline4);
set(splineLegend4,'Parent',figSpline4);
%%Add E-cell annotations
DuPlotWormData.addECellAnnotationsNucs(axh0);


close (ancestor(axh1,'figure'));
close (ancestor(axh2,'figure'));
close (ancestor(axh3,'figure'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figAx18cc = figure('Position',[0,0,500,500]);

ax18cc = Data18_failure_purge.plotCountsVsCounts(Data18_failure_purge.CY5_IND,Data18_failure_purge.TMR_IND,'r');
axis(ax18cc,[0 1200 0 1200]);
%title(ax18cc, ['613:WT elt-2 vs. gfp transcripts'],'FontSize',14,'FontWeight', 'bold');
ylabel(ax18cc,['gfp transcripts']);
xlabel(ax18cc,['elt-2 transcripts']);

