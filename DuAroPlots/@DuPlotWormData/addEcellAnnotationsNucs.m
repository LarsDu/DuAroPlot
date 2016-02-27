function [  ] = addEcellAnnotationsNucs(axesHandle,varargin)
%ADD E CELL ANNOTATIONS
%This static function will draw Ecell annotations on top of an existing
%plot where the X-axis specifies nuclei counts.
%It will draw a lines and captions at X EQUALS 7, 14, 44, and 100

%Set default parameters for optional inputs
%Good instructions at
%http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
numvarargs = length(varargin);
if numvarargs > 2
    error('addEcellAnnotationsNucs', ...
        'requires at most 3 optional inputs: ',...
        'line width, figure handle');
end

%Get the default figure handle from the axeHandle input
defaultFigHandle = ancestor(axesHandle,'figure')

optargs = {1,defaultFigHandle};

% now put these defaults into the valuesToUse cell array,
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;
% or ...
% [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names
%the user can also specify a figure handle
[lineWidth,figHandle] = optargs{:};


%%%%%%%%%Change Settings Here:%%%%%%%%%%%%%%%%
% y-Position
yPosE = 0.98;

%Textbox dimensions
w=0.005;
h=0.005;

%Manual offset
manXOff = 0.005; 

%Text color
colorEcell = [0.7 0.7 0.7]; %Light gray
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axesHandle1Pos = get(axesHandle,'Position');
xLimits = get(axesHandle,'xLim');
yLimits = get(axesHandle,'yLim');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%Draw E stage annotationlines
hold(axesHandle,'on');
%E-cell marker lines


yMin = min(yLimits)
yMax = max(yLimits)

line([7,7],[yMin,yMax],'Parent',axesHandle,'LineStyle','--','Color',colorEcell, ...
 'LineWidth',lineWidth); %1E

line([14,14],[yMin,yMax],'Parent',axesHandle,'LineStyle','--','Color',colorEcell,...
 'LineWidth',lineWidth); %2E

line([44,44],[yMin,yMax],'Parent',axesHandle,'LineStyle','--','Color',colorEcell,...
 'LineWidth',lineWidth); %4E


line([100,100],[yMin,yMax],'Parent',axesHandle,'LineStyle','--','Color',colorEcell,...
 'LineWidth',lineWidth); %8E



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Handle annotations
%Collect axesHandle position [x,y,w,h] in normalized figure units
xAxOff = axesHandle1Pos(1);
yAxOff = axesHandle1Pos(2);
widthAx= axesHandle1Pos(3);
heightAx = axesHandle1Pos(4);


%Calculate axesNormalized positions for each annotation
xScale = max(xLimits);
xPos1E = 7/xScale;
xPos2E = 14/xScale;
xPos4E = 44/xScale;
xPos8E = 100/xScale;

%Convert from axes normalized units to figure normalized units
%Apply manual offset as well

xPos1E = xAxOff+xPos1E*widthAx -manXOff;
xPos2E = xAxOff+xPos2E*widthAx -manXOff;
xPos4E = xAxOff+xPos4E*widthAx -manXOff;
xPos8E = xAxOff+xPos8E*widthAx -manXOff;


yPosE = +yAxOff+ yPosE*heightAx;

hold(axesHandle,'on');

an1 = annotation(figHandle,'textbox', [xPos1E,yPosE,w,h],...
           'String', '1E','LineStyle','none',...
           'FontSize', 9,...
           'Color', colorEcell,...
           'FontWeight','bold');

an2 = annotation(figHandle,'textbox', [xPos2E,yPosE,w,h],...
           'String', '2E','LineStyle','none',...
           'FontSize', 9,...
           'Color', colorEcell,...
           'FontWeight','bold');

an4 = annotation(figHandle,'textbox', [xPos4E,yPosE,w,h],...
           'String', '4E','LineStyle','none',...
           'FontSize', 9,...
           'Color', colorEcell,...
           'FontWeight','bold');

an8 = annotation(figHandle,'textbox', [xPos8E,yPosE,w,h],...
           'String', '8E','LineStyle','none',...
           'FontSize', 9,...
           'Color', colorEcell,...
           'FontWeight','bold');

end


%%Activate this feature by adding the line
%DuPlotWormData.addECellAnnotationsNucs(figHandle,axesHandle);
