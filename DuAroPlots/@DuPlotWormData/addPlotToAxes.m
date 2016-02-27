function [axesHandle] = addPlotToAxes (donorAxes, varargin)
%ADD PLOT TO AXES     This will add the object data
%from the axes handle axDonor to axRecieve
%If axReceive is not specified, this method will use
%the current axes gca instead


%Set default parameters for optional inputs
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 1
                error('addPlot to Axes', ...
                    'requires at most 1 optional inputs: ',...
                    'an axesHandle to plot data onto!');
            end
            
            optargs = {gca};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            % or ...
            % [optargs{1:numvarargs}] = varargin{:};
            
            % Place optional args in memorable variable names
            [receiverAxes] = optargs{:};
            
            copyobj(allchild(donorAxes),receiverAxes)       
end