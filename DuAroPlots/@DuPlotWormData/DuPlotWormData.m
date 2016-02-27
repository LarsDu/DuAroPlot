classdef DuPlotWormData < DuProcessWormData
    %DUPLOTWORMDATA     A class for encapsulating smFISH table data for
    %plotting, sorting, rearrangement, and manual curation.
    %
    %Requires DuProcessWormData, DuSpotViz (for some methods)
    %
    % Example usage:
    % MyPlot = DuPlotWormData(wormData)
    % myPlotAxesHandle = MyPlot.plotNucVsCount('tmr','gfp','Nuclei count vs. Pelt-2(613 bp promoter)::gfp','g')
    % MyPlot.plotMeanTraj(myPlotAxesHandle, MyPlot.NUCLEI_IND,MyPlot.TMR_IND,'blue')
    %
    %
    %
    %Helpful notes on subplotting:
    % Sometimes you will want to move some of these plots into
    % subplots, stacking them right on top of each other.
    % This can be achieved easily by the following:
    %
    %     %First plot your data, saving your axesHandles:
    %     ax1 = myData.plotNucVsCountScatter('tmr','gfp','title1','green')
    %     ax2 = myData.plotNucVsCountScatter('cy5','elt-2','title2,'red')
    %     %Open new figure and make subplots
    %     figure
    %     sax1 = subplot(2,1,1)
    %     sax2 = subplot(2,1,2)
    %     %Copy objects from ax1 and ax2 into sax1 and sax2
    %     copyobj(allchild(ax1),sax1)
    %     copyobj(allchild(ax2),sax2)
    
    
    properties (Access = public)
        
        strainName='';
        
        %Default Axes settings
        xmin=0;
        xmax=160;
        ymin=0;
        ymax=1200
        xtmin=0
        xtmax=175
        fontsize = 18;
    end %public properties
    
    
    
    methods
        function myPlot = DuPlotWormData(wormDataStruct,varargin)
            %Constructor
            %Initialize superclass constructor
            myPlot@DuProcessWormData(wormDataStruct);
            
            %%Handle default arguments
            numvarargs = length(varargin)
            optargs = {''}
            optargs(1:numvarargs) = varargin;
            myPlot.strainName = optargs{:};
            
        end %Constructor
        
        
        function [axesHandle] = plotNucVsCountErrorBar(myPlot,channel,...
                varargin)
            %PLOTNUCVSCOUNT     Plot nuclei counts versus
            %transcript counts with error bars
            %
            % Usage:
            % myPlotObj.plotNucVsCount('tmr', 'gfp', 'title', 'g');
            
            channelIndex = myPlot.channel.(channel); %IE: type 'tmr' to
            %get the tmr channel
            %Set default parameters for optional inputs
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 3
                error('Function', ...
                    'requires at most 3 optional inputs',...
                    '(transcriptName (for title),figureTitle, and color)');
            end
            
            optargs = {'','','g'};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            
            
            % Place optional args in memorable variable names
            [transcriptName,titleCard, color] = optargs{:};
            
            
            %Set up figure
            figHandle = figure('Position',[0,0,500,400]);
            axesHandle = axes('Parent', figHandle)
            
          
            
            %Draw the plot
            
            
            errorBarSeries = ...
                errorbar(myPlot.myWormDataStruct.spotNum(1:end, myPlot.NUCLEI_IND), ...
                myPlot.myWormDataStruct.spotNum(1:end,channelIndex), ...
                myPlot.myWormDataStruct.L(1:end, channelIndex-3), ...
                myPlot.myWormDataStruct.U(1:end, channelIndex-3),['g','.'] );
            
            errorbar_tick(errorBarSeries,1000);
            
            
            %Make sure the errorBar Graphics Object is attached to
            %the current axesHandle
            set(errorBarSeries,'Parent',axesHandle,'Color',color);
            
            %Set axes for specified axesHandle
            axis(axesHandle,[myPlot.xmin myPlot.xmax myPlot.ymin myPlot.ymax]);
            %Set fontsize
            set(axesHandle,'fontsize', myPlot.fontsize);
            
            %Label title and axes
            title(axesHandle,[titleCard],'FontSize',14,'FontWeight','bold');
            xlabel(axesHandle,'Number of nuclei','FontSize',14,'FontWeight','bold');
            ylabel(axesHandle,['Number of ' transcriptName ' transcripts'],...
                'FontSize',14,'FontWeight','bold');
            hold(axesHandle, 'on');
            
            %%Add E-cell annotations
            %DuPlotWormData.addECellAnnotationsNucs(figHandle,axesHandle);
            
            
        end % plotNucVsCount
        
        
        
        
        function [axesHandle] = plotNucVsCountScatter(myPlot,channel, varargin)
            %PLOT NUC VS COUNT SCATTER
            %Plot nuclei counts versus transcript
            %counts as a scatter plot
            % Usage:
            % myPlotObj.plotNucVsCount('tmr', 'gfp', 'g');
            
            channelIndex = myPlot.channel.(channel); %IE: type 'tmr' to
            %get the tmr channel
            
            
            %Set default parameters for optional inputs
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 3
                error('Function', ...
                    'requires at most 3 optional inputs',...
                    '(transcriptName (for title),figureTitle, and color)');
            end
            
            optargs = {'','','g'};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            
            
            % Place optional args in memorable variable names
            [transcriptName,titleCard, color] = optargs{:};
            
            
            figHandle = figure('Position',[0,0,500,400]);
            axesHandle = axes('Parent', figHandle)
            

            
            %Draw the plot
            %Note this error bar command should be altered to be more
            %robust
            
            scatterSeriesObj = ...
                scatter(myPlot.myWormDataStruct.spotNum(1:end, myPlot.NUCLEI_IND), ...
                myPlot.myWormDataStruct.spotNum(1:end,channelIndex), ...
                15,...
                'MarkerEdgeColor',color,...
                'MarkerFaceColor',color);
     
            
            
            %Make sure the errorBar Graphics Object is attached to
            %the current axesHandle
            set(scatterSeriesObj,'Parent',axesHandle);
            
            %Set axes for specified axesHandle
            axis(axesHandle,[myPlot.xmin myPlot.xmax myPlot.ymin myPlot.ymax]);
            
            
            %Label title and axes
            title(axesHandle,[titleCard],'FontSize',14,'FontWeight','bold');
            xlabel(axesHandle,'Number of nuclei','FontSize',14,'FontWeight','bold');
            ylabel(axesHandle,['Number of ' transcriptName ' transcripts'],...
                'FontSize',14,'FontWeight','bold');
            hold(axesHandle, 'on');
          
            %%Add E-cell annotations
            %DuPlotWormData.addECellAnnotationsNucs(figHandle,axesHandle);
            
        end % plotNucVsCountScatter
        
        
        
        
        function [axesHandle] = overlayNucVsCount(myPlot,axesHandle, ...
                channel,varargin)
            %OVERLAY NUC VS COUNT
            %Overlay nuclei vs count data onto an existing axes
            %This plot will have error bars
            %
            %Usage:
            %myPlotObj.overlayNucVsCount(gca,'tmr','gfp','green')
            
            channelIndex = myPlot.channel.(channel);
            
            
            %Set default parameters for optional inputs
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 2
                error('Function', ...
                    'requires at most 2 optional inputs',...
                    '(transcriptName (for title) and color)');
            end
            
            optargs = {'','g'};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            
            
            % Place optional args in memorable variable names
            [transcriptName, color] = optargs{:};
            
            
            %Draw the plot
            %Note this error bar command should be altered to be more
            %robust
            
            errorBarSeries = ...
                errorbar(myPlot.myWormDataStruct.spotNum(1:end, myPlot.NUCLEI_IND), ...
                myPlot.myWormDataStruct.spotNum(1:end,channelIndex), ...
                myPlot.myWormDataStruct.L(1:end, channelIndex-3), ...
                myPlot.myWormDataStruct.U(1:end, channelIndex-3),[color,'.'] );
            
            errorbar_tick(errorBarSeries,1000);
            
            %Make sure the errorBar Graphics Object is attached to
            %the current axesHandle
            set(errorBarSeries,'Parent',axesHandle);
            
            %Set axes for specified axesHandle
            axis(axesHandle,[myPlot.xmin myPlot.xmax myPlot.ymin myPlot.ymax]);
            %Set fontsize
            set(axesHandle,'FontSize', myPlot.fontsize,'MarkerEdgeColor',color);
            
            %Label title and axes
            %title(axesHandle,[titleCard],'FontSize',14,'FontWeight','bold');
            xlabel(axesHandle,'Number of nuclei','FontSize',14,'FontWeight','bold');
            ylabel(axesHandle,['Number of ' transcriptName ' transcripts'],'FontSize',14,'FontWeight','bold');
            hold(axesHandle, 'on');
            
        end % overlayNucVsCount
        
        
        
        
        
        function [axesHandle] = overlayNucVsCountScatter(myPlot,axesHandle, ...
                channel, varargin)
            %OVERLAY NUC VS COUNT SCATTER
            %Overlay nuclei vs count data onto an existing axes
            %This plot will be in the form of a scatterplot
            %
            % Usage:
            % myPlotObj.overlayNucVsCount(gca,'tmr','gfp','green')
            
            
            channelIndex = myPlot.channel.(channel); %IE: type 'tmr' to
            %get the tmr channel
            
            
            %Set default parameters for optional inputs
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 2
                error('Function', ...
                    'requires at most 2 optional inputs ',...
                    '(transcriptName, color)');
            end
            
            optargs = {'','g'};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            
            
            % Place optional args in memorable variable names
            [transcriptName, color] = optargs{:};
            
            
            
            
            %Draw the plot
            %Note this error bar command should be altered to be more
            %robust
            
            scatterSeriesObj = ...
                scatter(myPlot.myWormDataStruct.spotNum(1:end, myPlot.NUCLEI_IND), ...
                myPlot.myWormDataStruct.spotNum(1:end,channelIndex), ...
                15,'MarkerEdgeColor',color,'MarkerFaceColor',color);
            
            %Make sure the scatterSeriesObj is attached to
            %the current axesHandle
            set(scatterSeriesObj,'Parent',axesHandle);
            
            %Set axes for specified axesHandle
            axis(axesHandle,[myPlot.xmin myPlot.xmax myPlot.ymin myPlot.ymax]);
            %Set fontsize
            set(axesHandle,'fontsize');
            
            
            
            %Label title and axes
            %title(axesHandle,[titleCard]);
            xlabel(axesHandle,'Number of nuclei','FontSize',14,'FontWeight','bold');
            ylabel(axesHandle,['Number of ' transcriptName ' transcripts'],'FontSize',14,'FontWeight','bold');
            hold(axesHandle, 'on');
            
        end % overlayNucVsCountScatter
        
        
        
        
        
        function [axesHandle] = plotTimeVsCount(myPlot,channel,varargin)
            %PLOT TIME VS COUNT
            %Plot developmental time in minutes versus transcript
            %count with error bars.
            %
            % Usage:
            % myPlotObj.plotTimeVsCount('tmr','gfp','green')
            
            channelIndex = myPlot.channel.(channel); %IE: type 'tmr' to
            %get the tmr channel
            %Set default parameters for optional inputs
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 3
                error('Function', ...
                    'requires at most 3 optional inputs',...
                    '(transcriptName (for title), figureTitle, and color)');
            end
            
            optargs = {'','','g'};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            
            
            % Place optional args in memorable variable names
            [transcriptName, titleCard, color] = optargs{:};
            
            
            % Set up the figure
            figHandle = figure('Position',[0,0,500,400]);
            axesHandle = axes('Parent', figHandle)
            
            %Draw the plot
            %Note this error bar command should be altered to be more
            %robust
            errorBarSeries = ...
                errorbar(myPlot.myWormDataStruct.spotNum(1:end, myPlot.TIME_IND), ...
                myPlot.myWormDataStruct.spotNum(1:end,channelIndex), ...
                myPlot.myWormDataStruct.L(1:end, channelIndex-3), ...
                myPlot.myWormDataStruct.U(1:end, channelIndex-3),[color,'.'] );
            
            errorbar_tick(errorBarSeries,1000);
            
            %Make sure the errorBar Graphics Object is attached to
            %the current axesHandle
            set(errorBarSeries,'Parent',axesHandle);
            
            %Set axes for specified axesHandle
            axis(axesHandle,[myPlot.xtmin myPlot.xtmax myPlot.ymin myPlot.ymax]);
            %Set fontsize
            set(axesHandle,'fontsize', myPlot.fontsize);
            
            
            %Label title and axes
            title(axesHandle,[titleCard],'FontSize',14,'FontWeight','bold');
            xlabel(axesHandle,'Time in minutes (extrapolated from # of nuclei)','FontSize',14,'FontWeight','bold');
            ylabel(axesHandle,['Number of ' transcriptName ' transcripts'],'FontSize',14,'FontWeight','bold');
            hold(axesHandle,'on');
        end % plotTimeVsCount
        
        
        
        
        
        function [axesHandle] = plotTimeVsCountScatter(myPlot,channel,varargin)
            %PLOT TIME VS COUNT SCATTER
            %Plot developmental time in minutes versus transcript
            %count as a scatterplot
            %
            % Usage:
            % myPlotObj.plotTimeVsCount('tmr','gfp','green')
            
            channelIndex = myPlot.channel.(channel); %IE: type 'tmr' to
            %get the tmr channel
            %Set default parameters for optional inputs
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 3
                error('Function', ...
                    'requires at most 3 optional inputs',...
                    '(transcriptName (for title), figureTitle, and color)');
            end
            
            optargs = {'','','g'};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            
            
            % Place optional args in memorable variable names
            [transcriptName, titleCard, color] = optargs{:};
            
            
            figHandle = figure('Position',[0,0,500,400]);
            axesHandle = axes('Parent', figHandle)
            
            %Draw the plot
            %Note this error bar command should be altered to be more
            %robust
            scatterSeriesObj = ...
                scatter(myPlot.myWormDataStruct.spotNum(1:end, myPlot.TIME_IND), ...
                myPlot.myWormDataStruct.spotNum(1:end,channelIndex), ...
                20,color,'o','filled');
            
            %Make sure the scatterSeriesObj is attached to
            %the current axesHandle
            set(scatterSeriesObj,'Parent',axesHandle);
            
            %Set axes for specified axesHandle
            axis(axesHandle,[myPlot.xtmin myPlot.xtmax myPlot.ymin myPlot.ymax]);
            %Set fontsize
            set(axesHandle,'fontsize', myPlot.fontsize);
            
            
            %Label title and axes
            title(axesHandle,[titleCard],'FontSize',14,'FontWeight','bold');
            xlabel(axesHandle,'Time in minutes (extrapolated from # of nuclei)','FontSize',14,'FontWeight','bold');
            ylabel(axesHandle,['Number of ' transcriptName ' transcripts'],'FontSize',14,'FontWeight','bold');
            hold(axesHandle,'on');
        end % plotTimeVsCountScatter
        
        
        
        
        function [axesHandle] = plotMeanTrajSavitskyGolay(myPlot, axesHandle, xCol, yCol, varargin)
            %PLOT MEAN TRAJ SAVITSKY GOLAY
            %Plot smooth mean trajectory curve overlaying
            %existing axes.
            %
            % Usage:
            % myPlotObj.plotMeanTrajSavitskyGolay(gca,myPlotObj.NUCLEI_IND,myPlotObj.TMR_IND,'green')
            % You enter something like 'myPlotObj.NUCLEI_IND' for xCol
            % or the raw column number.
            
            color = varargin{1};
            countsY = myPlot.meanTrajectorySavitskyGolayFromCol(xCol,yCol);
            countsX = myPlot.myWormDataStruct.spotNum(1:end,xCol);
            plot(axesHandle,countsX, countsY,'Color',color);
            hold(axesHandle,'on');
        end % plotPointsToAxes
        
        
                
        
        
        function [axesHandle] = plotCoeffVar(myPlot, axesHandle, xCol,yCol,varargin)
            %PLOT COEFF VAR
            %Plot Coefficient of Variation curve
            %over existing axes
            %
            % Usage:
            % myPlotObj.plotCoeffVar(gca,myPlotObj.NUCLEI_IND,myPlotObj.TMR_IND,'green')
            % You enter something like 'myPlotObj.NUCLEI_IND' for xCol
            % or the raw column number.
            
            color = varargin{1};
            dataX=myPlot.myWormDataStruct.spotNum(1:end,xCol)
            cvY = myPlot.coeffVarFromCol(xCol,yCol);
            plot(axesHandle, dataX, cvY,'Color',color);
            hold(axesHandle,'on');
            axis(axesHandle,[80 170 0 2])
            set(axesHandle,'FontSize',myPlot.fontsize,'LineWidth',4);
            xlabel(axesHandle,['Time in Minutes (extrapolated from # ' ...
                'of nuclei']);
            ylabel(axesHandle,['Coefficient of Variation']);
        end
        
                
                
   
        
        
        
        function [axesHandle] = plotCountsVsCounts(myPlot,xCol,yCol,varargin)
            %PLOT COUNTS VS COUNTS
            %
            %Plot transcript counts for different channels against each
            %other. Useful for looking for correlations in gene
            %expression.
            %
            % Usage:
            % myPlotObj.plotCountsVsCounts(myPlotObj.A594_IND,myPlotObj.TMR_IND,'green')
            
            axesHandle = axes;
            %Set default parameters for optional inputs
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 1
                error('Function', ...
                    'requires at most 1 optional inputs ',...
                    '(color)');
            end
            
            optargs = {'g'};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            
            %title(axesHandle,[myPlot.strainName]);
            % Place optional args in memorable variable names
            [color] = optargs{:};
            
            xData = myPlot.myWormDataStruct.spotNum(1:end,xCol)
            yData = myPlot.myWormDataStruct.spotNum(1:end,yCol)
            scatter(xData,yData,15,'MarkerEdgeColor',color,'MarkerFaceColor',color);
            axis(axesHandle,[0 1000 0 1000]);
            %Set fontsize
            set(axesHandle,'fontsize', myPlot.fontsize);
            
            %xlabel(axesHandle,['# of ' xCol 'transcripts']);
            %ylabel(axesHandle,['']);
            hold(axesHandle,'on');
        end %plotCountsVsCounts
        
        
        
        
        
        function [axesHandle] = plotEcellsVsCount(myPlot,channel, ...
                varargin)
            %PLOT E CELLS VS COUNT
            %Plot endoderm cells (E cells) versus transcript counts
            %with error bars.
            %
            % Usage:
            % myPlotObj.plotEcellsVsCount('tmr','gfp','green')
            
            channelIndex = myPlot.channel.(channel); %IE: type 'tmr' to
            %get the tmr channel
            %Set default parameters for optional inputs
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 3
                error('Function', ...
                    'requires at most 3 optional inputs',...
                    '(transcriptName (for title), figureTitle, and color)');
            end
            
            optargs = {'','','g'};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            
            
            % Place optional args in memorable variable names
            [transcriptName, titleCard, color] = optargs{:};
            
            
            
            figHandle = figure('Position',[0,0,500,400]);
            axesHandle = axes('Parent', figHandle)
            
            %Draw the plot
            %Note this error bar command should be altered to be more
            %robust
            
            errorBarSeries = ...
                errorbar(myPlot.myWormDataStruct.spotNum(1:end, myPlot.ECELL_IND), ...
                myPlot.myWormDataStruct.spotNum(1:end,channelIndex), ...
                myPlot.myWormDataStruct.L(1:end, channelIndex-3), ...
                myPlot.myWormDataStruct.U(1:end, channelIndex-3),[color,'.'] );
            
            errorbar_tick(errorBarSeries,1000);
            
            %Make sure the errorBar Graphics Object is attached to
            %the current axesHandle
            set(errorBarSeries,'Parent',axesHandle);
            
            %Set axes for specified axesHandle
            axis(axesHandle,[0 20 myPlot.ymin myPlot.ymax]);
            %Set fontsize
            set(axesHandle,'fontsize', myPlot.fontsize);
            
            
            %Label title and axes
            title(axesHandle,[titleCard],'FontSize',14,'FontWeight','bold');
            xlabel(axesHandle,'Number of E cells','FontSize',14,'FontWeight','bold');
            ylabel(axesHandle,['Number of ' transcriptName ' transcripts'],'FontSize',14,'FontWeight','bold');
            hold(axesHandle, 'on');
        end % plotEcellsVsCount
        
        
        
        
        function [axesHandle] = plotEcellsVsCountScatter(myPlot, ...
                channel,...
                varargin)
            
            %PLOT E CELLS VS COUNT SCATTER
            %Plot endoderm cells (E cells) versus transcript counts
            %as scatterplot
            %
            % Usage:
            % myPlotObj.plotEcellsVsCount('tmr','gfp','green')
            
            channelIndex = myPlot.channel.(channel); %IE: type 'tmr' to
            %get the tmr channel
            
            %Set default parameters for optional inputs
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 3
                error('Function', ...
                    'requires at most 3 optional inputs',...
                    '(transcriptName (for title), figureTitle, and color)');
            end
            
            optargs = {'','','g'};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            
            
            % Place optional args in memorable variable names
            [transcriptName, titleCard, color] = optargs{:};

            
            %Set up the plot
            figHandle = figure('Position',[0,0,500,400]);
            axesHandle = axes('Parent', figHandle)
            
            %Draw the plot
            %Note this error bar command should be altered to be more
            %robust
            
            scatterSeriesObj = ...
                scatter(myPlot.myWormDataStruct.spotNum(1:end, myPlot.ECELL_IND), ...
                myPlot.myWormDataStruct.spotNum(1:end,channelIndex), ...
                20, ...
                color,'o','filled' );
            
            
            %Make sure the scatterSeriesObj Graphics Object is attached to
            %the current axesHandle
            set(scatterSeriesObj,'Parent',axesHandle);
            
            %Set axes for specified axesHandle
            axis(axesHandle,[0 20 myPlot.ymin myPlot.ymax]);
            %Set fontsize
            set(axesHandle,'fontsize', myPlot.fontsize);
            
            %Label title and axes
            title(axesHandle,[titleCard],'FontSize',14,'FontWeight','bold');
            xlabel(axesHandle,'Number of E cells','FontSize',14,'FontWeight','bold');
            ylabel(axesHandle,['Number of ' transcriptName ' transcripts'],'FontSize',14,'FontWeight','bold');
            hold(axesHandle, 'on');
        end % plotEcellsVsCountScatter
        
        
        
        
        
        function [axesHandle] = plotNucVsCountScatterSpline(myPlot,channel,varargin)
            %PLOT NUC VS COUNT SCATTER
            %Plot nuclei counts versus transcript
            %counts as a scatter plot with mean trajectory
            %spline.
            % Usage:
            % myPlotObj.plotNucVsCount('tmr', 'gfp', 'g');
            % Optional parameters(in order): 'title','color',smoothingParameter
            
            channelIndex = myPlot.channel.(channel); %IE: type 'tmr' to
            %get the tmr channel
            
            
            %Set default parameters for optional inputs
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 4
                error('plotNucVsCountScatterSpline', ...
                    'requires at most 4 optional inputs: ',...
                    'transcriptName, titleCard, color, smoothingParam');
            end
            
            optargs = {'','', 'blue', 0.0001};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            % or ...
            % [optargs{1:numvarargs}] = varargin{:};
            
            % Place optional args in memorable variable names
            [transcriptName,titleCard,color, smoothingParam] = optargs{:};
            
            
            
            %%Set up the figure
            
            figHandle = figure('Position',[0,0,500,400]);
            axesHandle = axes('Parent', figHandle);
            

            
            %Draw the scatter plot
            
            scatterSeriesObj = ...
                scatter(myPlot.myWormDataStruct.spotNum(1:end, myPlot.NUCLEI_IND), ...
                myPlot.myWormDataStruct.spotNum(1:end,channelIndex), ...
                15, ...
                'MarkerFaceColor',color,'MarkerEdgeColor',color);
    
            hold(axesHandle,'on');
            
            %Calculate spline
            
           
            
            xNucs = myPlot.myWormDataStruct.spotNum(:,myPlot.NUCLEI_IND);
            yCounts = myPlot.myWormDataStruct.spotNum(:,channelIndex);
            
             %%Only calculate spline for values between xmin and xmax
            %withinRange = (xNucs > myPlot.xmin).*(xNucs < myPlot.xmax);
            %%Delete rows outside of range
            %xNucs(not(withinRange)) = [];
            %yCounts(not(withinRange)) = [];
            
            splineY = ...
                myPlot.meanTrajectorySmoothingSpline(xNucs,yCounts,smoothingParam);
            splineX = 1:length(xNucs);
            %Draw spline
            splinePlot = plot(splineX,splineY,'Color',color,'LineWidth',2);
            hold(axesHandle,'off');
            
            
            
            
            %Make sure the scatterSeries Graphics Object is attached to
            %the current axesHandle
            set(scatterSeriesObj,'Parent',axesHandle);
            set(splinePlot, 'Parent', axesHandle);
            %Set axes for specified axesHandle
            axis(axesHandle,[myPlot.xmin myPlot.xmax myPlot.ymin myPlot.ymax]);
            
            
            %Label title and axes
            title(axesHandle,[titleCard],'FontSize',14,'FontWeight','bold');
            xlabel(axesHandle,'Number of nuclei','FontSize',14,'FontWeight','bold');
            ylabel(axesHandle,['Number of ' transcriptName ' transcripts'],...
                'FontSize',14,'FontWeight','bold');
            hold(axesHandle, 'on');
            
            %%Add E-cell annotations
            %DuPlotWormData.addECellAnnotationsNucs(figHandle,axesHandle);
            
        end
        
        
        
        function [axesHandle] = plotNucVsCountErrorBarSpline(myPlot,channel,varargin)
            %PLOT NUC VS COUNT SCATTER
            %Plot nuclei counts versus transcript
            %counts as a scatter plot with mean trajectory
            %spline.
            % Usage:
            % myPlotObj.plotNucVsCount('tmr', 'gfp', 'g');
            % Optional parameters(in order): 'title','color',smoothingParameter
            
            channelIndex = myPlot.channel.(channel); %IE: type 'tmr' to
            %get the tmr channel
            
            
            %Set default parameters for optional inputs
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 4
                error('plotNucVsCountScatterSpline', ...
                    'requires at most 4 optional inputs: ',...
                    'transcriptName, titleCard, color, smoothingParam');
            end
            
            optargs = {'','', 'blue', 0.0001};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            % or ...
            % [optargs{1:numvarargs}] = varargin{:};
            
            % Place optional args in memorable variable names
            [transcriptName,titleCard,color, smoothingParam] = optargs{:};
            
            
            
             %Set up figure
            figHandle = figure('Position',[0,0,500,400]);
            axesHandle = axes('Parent', figHandle)
            
           
            
            %Draw the plot
            
            
            errorBarSeries = ...
                errorbar(myPlot.myWormDataStruct.spotNum(1:end, myPlot.NUCLEI_IND), ...
                myPlot.myWormDataStruct.spotNum(1:end,channelIndex), ...
                myPlot.myWormDataStruct.L(1:end, channelIndex-3), ...
                myPlot.myWormDataStruct.U(1:end, channelIndex-3),['g','.'] );
            
            errorbar_tick(errorBarSeries,1000);
            
            
            %Make sure the errorBar Graphics Object is attached to
            %the current axesHandle
            set(errorBarSeries,'Parent',axesHandle,'Color',color);
            hold(axesHandle,'on'); 
            
            
            
            %Calculate spline
            
           
            
            xNucs = myPlot.myWormDataStruct.spotNum(:,myPlot.NUCLEI_IND);
            yCounts = myPlot.myWormDataStruct.spotNum(:,channelIndex);
            
             %%Only calculate spline for values between xmin and xmax
            %withinRange = (xNucs > myPlot.xmin).*(xNucs < myPlot.xmax);
            %%Delete rows outside of range
            %xNucs(not(withinRange)) = [];
            %yCounts(not(withinRange)) = [];
            
            splineY = ...
                myPlot.meanTrajectorySmoothingSpline(xNucs,yCounts,smoothingParam);
            splineX = 1:length(xNucs);
            %Draw spline
            splinePlot = plot(splineX,splineY,'Color',color,'LineWidth',2);
            hold(axesHandle,'off');
            
            
            
            
            %Make sure the errorBar Graphics Object is attached to
            %the current axesHandle
       
            set(splinePlot, 'Parent', axesHandle);
            %Set axes for specified axesHandle
            axis(axesHandle,[myPlot.xmin myPlot.xmax myPlot.ymin myPlot.ymax]);
            
            
            %Label title and axes
            title(axesHandle,[titleCard],'FontSize',14,'FontWeight','bold');
            xlabel(axesHandle,'Number of nuclei','FontSize',14,'FontWeight','bold');
            ylabel(axesHandle,['Number of ' transcriptName ' transcripts'],...
                'FontSize',14,'FontWeight','bold');
            hold(axesHandle, 'on');
            
            %%Add E-cell annotations
            %DuPlotWormData.addECellAnnotationsNucs(figHandle,axesHandle);
            
        end
        
        
        
        function [axesHandle] = plotNucVsCountSpline(myPlot,channel,varargin)
            %PLOT NUC VS COUNT SCATTER
            %Plot nuclei counts versus transcript
            %counts as a scatter plot with mean trajectory
            %spline.
            % Usage:
            % myPlotObj.plotNucVsCount('tmr', 'gfp', 'g');
            % Optional parameters(in order): 'title','color',smoothingParameter
            
            channelIndex = myPlot.channel.(channel); %IE: type 'tmr' to
            %get the tmr channel
            
            
            %Set default parameters for optional inputs
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 4
                error('plotNucVsCountScatterSpline', ...
                    'requires at most 4 optional inputs: ',...
                    'transcriptName, titleCard, color, smoothingParam');
            end
            
            optargs = {'','', 'blue', 0.0001};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            % or ...
            % [optargs{1:numvarargs}] = varargin{:};
            
            % Place optional args in memorable variable names
            [transcriptName,titleCard,color, smoothingParam] = optargs{:};
            
            
            
            %%Set up the figure
            
            figHandle = figure('Position',[0,0,500,400]);
            axesHandle = axes('Parent', figHandle);


            
          %Calculate spline
            
           
            
            xNucs = myPlot.myWormDataStruct.spotNum(:,myPlot.NUCLEI_IND);
            yCounts = myPlot.myWormDataStruct.spotNum(:,channelIndex);
             %%Only calculate spline for values between xmin and xmax
           % withinRange = (xNucs > myPlot.xmin).*(xNucs < myPlot.xmax);
            %Delete rows outside of range
            %xNucs(not(withinRange)) = [];
            %yCounts(not(withinRange)) = [];
            splineY = ...
                myPlot.meanTrajectorySmoothingSpline(xNucs,yCounts,smoothingParam);
            splineX = 1:length(xNucs);
            %Draw spline
            splinePlot =  plot(splineX,splineY,'Color',color,'LineWidth',2);
            hold(axesHandle,'off');
            
            
            
            
            %Make sure the errorBar Graphics Object is attached to
            %the current axesHandle
            set(splinePlot,'Parent',axesHandle);
            
            %Set axes for specified axesHandle
            axis(axesHandle,[myPlot.xmin myPlot.xmax myPlot.ymin myPlot.ymax]);
            
            
            %Label title and axes
            title(axesHandle,[titleCard],'FontSize',14,'FontWeight','bold');
            xlabel(axesHandle,'Number of nuclei','FontSize',14,'FontWeight','bold');
            ylabel(axesHandle,['Mean Number of ' transcriptName ' transcripts'],...
                'FontSize',14,'FontWeight','bold');
            hold(axesHandle, 'on');
            
            %%Add E-cell annotations
            %DuPlotWormData.addECellAnnotationsNucs(figHandle,axesHandle);
            
           
        end
        
        
        function [axesHandle] = plotEctopicExpEstimates(myPlot,varargin)
           %Calls DuProcessWormData.quantifyEctopicExp
           %and creates a plot
           
           titleCard = strcat(myPlot.strainName,' ectopic gfp estimates');
            numvarargs = length(varargin);
            if numvarargs > 2
                error('plotNucVsCountScatterSpline', ...
                    'requires at most 2 optional inputs: ',...
                    ' color, titleCard');
            end
            
            optargs = {'b','Title'};
            
  
            optargs(1:numvarargs) = varargin;

            [color,titleCard] = optargs{:};
            
           ectopicMat = myPlot.quantifyEctopicExp();
           %figHandle = figure('Position',[0,0,500,400]);
           %axesHandle = axes('Parent', figHandle);
           %axesHandle = axes();
           plotHandle = plot(ectopicMat(:,2),ectopicMat(:,3),'.','Color',color);
           %set(plotHandle, 'Parent',axesHandle)
           axis([0 100 0 100]);
           title(titleCard,'FontSize',14,'FontWeight','bold');
           
        end
        
        
   
        
        
        
        
        %function [] = saveDuPlotWormData(myPlot,filename)
            %DO NOT USE THIS METHOD YET!
         %   save ( filename,getobjname(myPlot));
            
        %end%saveDuPlotWormData(DO NOT USE YET!)
        
        
        
        %%%%%%%%%Get and set functions go here%%%%%%%%%%%%%%%
        %Note: I am not using matlab's standard get set syntax
        
        function []  = setCommonAxes (myPlot,xMin,xMax,yMin,yMax)
            %SET COMMON AXES
            %A quick way to set the common axes for nuclei
            %vs. transcript plots
            % Usage:
            % myPlotObj.setCommonAxes(xmin,xmax,ymin,ymax);
            
            myPlot.xmin = xMin;
            myPlot.xmax = xMax;
            myPlot.ymin = yMin;
            myPlot.ymax = yMax;
            
        end %setCommonAxes
        
        function [] = setTimeAxes (myPlot,tXmin,tXmax)
            %SET TIME AXES
            %A quick way to set the upper and lower x bounds for time
            %vs. transcript plots
            % Usage:
            % myPlotObj.setTimeAxes(txmin,txmax);
            
            myPlot.xtmin = tXmin;
            myPlot.xtmax = tXmax;
        end %setCommonAxes

        
        
        
        
        
    end %non-static methods
    
    
    
       
    
    
    
    
    
    %%Static methods
    methods(Static = true)
        %Invoke static methods by typing DuPlotWormData.myStaticMethod()
        %instead of calling the object!
        
        %%The following methods should be treated as static
        %%They do not necessarily act on the current object's data,
        %%but rather, are placed in this class for the sake of convenience
        
        %Permutation test
        [axesHandle1] = compareTwoDataClouds(myPlot, ...
                DuPlotData1,...
                channel1,...
                DuPlotData2, ...
                channel2, ...
                varargin)
       

            
       %Plot the permutation test as a subplot along with a comparison
       %of mean trajectories
       [hax1, cloudAx] = plotCloudsAndMeans(Data1,channel1,Data2,...
           channel2,varargin);
           
        %Making subplots
       [figHandle] = subplotTwoFigures(axes1,axes2);
       
       %Add different plots to the same axes by copying
       [axesHandle3] = addPlotToAxes(donorAxes,varargin);
       
       %Add Ecell position annotations to nuclei vs. transcript plots
       addECellAnnotationsNucs(figHandle,axesHandle,varargin);

       
       
              

       
       
        
    end% static methods
    
end %end classdef



