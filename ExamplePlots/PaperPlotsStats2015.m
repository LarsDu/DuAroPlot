%Color spec for different strains

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

%%%%%%Plot mean Trajectory via splines

%Data05 vs 4G,11G,613:WT

DuPlotWormData.plotCloudsAndMeans(Data05,'tmr',Data06,'tmr',...
    Color05, Color06 ,0.0001,'gfp','gfp')

DuPlotWormData.plotCloudsAndMeans(Data05,'tmr',Data12,'tmr',...
    Color05, Color12 ,0.0001,'gfp','gfp')

DuPlotWormData.plotCloudsAndMeans(Data05,'tmr',Data18_failure_purge,'tmr',...
    Color05, Color18 ,0.0001,'gfp','gfp')



%613:WT vs 4G,11G,613:WT

DuPlotWormData.plotCloudsAndMeans(Data18_failure_purge,'tmr',Data25_failure_purge,'tmr',...
    Color18, Color25 ,0.0001,'gfp','gfp')



%%Additional plots
%1879:WT vs 1879:1S

DuPlotWormData.plotCloudsAndMeans(Data05,'tmr',Data07,'tmr',...
    Color05, Color07 ,0.0001,'gfp','gfp')


%613:WT vs 613:1T
DuPlotWormData.plotCloudsAndMeans(Data18_failure_purge,'tmr',Data28,'tmr',...
    Color18, Color28 ,0.0001,'gfp','gfp')


%Pseudo replicates -- These strains suggest more tightly controlled
%replicates would reveal no significant differences between data sets.

DuPlotWormData.plotCloudsAndMeans(Data05,'tmr',Data28,'tmr',...
    Color05, Color28 ,0.0001,'gfp','gfp')

%DuPlotWormData.plotCloudsAndMeans(Data18_failure_purgept2,'tmr',Data05,'tmr',...
%    Color05, Color28 ,0.0001,'gfp','gfp')

DuPlotWormData.plotCloudsAndMeans(DataElt2,'cy5',Data05,'tmr',ColorElt2,Color05,0.0001,'elt-2','gfp')




