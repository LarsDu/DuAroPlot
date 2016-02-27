function [figHandle] = subplotTwoFigures(axes1,axes2)
%SUBPLOT TWO FIGURES
%
%Takes the axes from two open plots and subplots them onto a
%single figure with axes1 on top of axes2.
% Usage:
% myDataObj.subplotTwoFigures(ax1,ax2)

%New figure
figHandle = figure;
sax1 = subplot(2,1,1)
sax2 = subplot(2,1,2)
%Set axes
%axis(sax1,[myPlot.xmin,myPlot.xmax,myPlot.ymin, ...
%    myPlot.ymax]);
%axis(sax2,[myPlot.xmin,myPlot.xmax,myPlot.ymin,myPlot.ymax]);

%Transfer plot data
copyobj(allchild(axes1),sax1);
copyobj(allchild(axes2),sax2);


end