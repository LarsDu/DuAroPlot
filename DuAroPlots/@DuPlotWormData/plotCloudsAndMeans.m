function [hax1,cloudAx] = plotCloudsAndMeans(Data1,channel1,...
                                              Data2,channel2,...
                                              varargin)
%SUBPLOT CLOUDS UNDER MEAN     Subplots plotNucVsCountSpline plot
% on top of a compareTwoDataClouds plot, with full E cell annotations
% The optional inputs are as follows:
% color1, color2, smoothingParam(default 0.0001),transcriptName1,
% transcriptName2, 'titleCard'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set default parameters for optional inputs
%Good instructions at
%http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
numvarargs = length(varargin);
if numvarargs > 6
    error('plotCloudsAndMeans', ...
        'requires at most 6 optional inputs: ',...
        'color1','color2','smoothingParam','transcriptName1',...
        'transcriptName2');
end



optargs = {'cyan','magenta',0.0001,'',''};

% now put these defaults into the valuesToUse cell array,
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;
% or ...
% [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names
[color1,color2,smoothingParam,transcriptName1,...
    transcriptName2] = optargs{:};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

titleCard = [Data1.strainName ' ' transcriptName1 ' vs. ',...
    Data2.strainName ' ' transcriptName2 ' transcripts'];
    
titleCardPerm = [Data1.strainName ' ' transcriptName1 ' vs. ',...
    Data2.strainName ' ' transcriptName2 ' expression difference'];


hax1 = Data1.plotNucVsCountScatterSpline(channel1,'',...
                         '',color1,smoothingParam);
                     
lineHand1= findobj(hax1,'Type','Line')

data2ax = Data2.plotNucVsCountScatterSpline(channel2,'',...
                         '',color2,smoothingParam);
lineHand2 = findobj(data2ax,'Type','Line')
                     
DuPlotWormData.addPlotToAxes(data2ax,hax1);


%Note: lineHandles is necessary for ensuring legend gets drawn in 
%correct order
lineHandles = [lineHand1,lineHand2];


leg1 = legend(lineHandles,[Data1.strainName ' ' transcriptName1],...
    [Data2.strainName ' ' transcriptName2],'FontSize',10)

set(leg1, 'Position',[0.25,0.75,0.1,0.1]);


DuPlotWormData.addECellAnnotationsNucs(hax1);

close (ancestor(data2ax,'figure'));

cloudAx = DuPlotWormData.compareTwoDataClouds(Data1,channel1,Data2,channel2);
DuPlotWormData.addECellAnnotationsNucs(cloudAx);

title(hax1,titleCard,'FontSize',14,'FontWeight','bold');
title(cloudAx,titleCardPerm,'FontSize',12,'FontWeight','bold');



end

