classdef DuProcessWormData < handle
%DUPROCESSWORMDATA      A class containing functions for sorting, 
%rearrangement, and manual curation of wormData
%structs from AroSpotFindingSuite. 
% Usage:
% MyPlot = DuProcessWormData(wormData)
% Note: I recommend using the subclass DuPlotWormData instead
% DuPlotWormData inherits all methods from this class, and 
% has methods for plotting.     

    properties (Access = public) %usable by subclasses
        myWormDataStruct;
        sgolaySpan; %Set using constructor
        
        blacklisted; %List of rows of wormDataStruct.spotNum
                     %removed by blacklistWorm method
        notes={}; %User notes on manual data parsing.
    end %End private properties
    
    properties ( Constant= true)
        %Note, these are all constant values that have an integer
        %value corresponding to columns in the wormData struct.
        %This is all purely for ease of readability
        %Future changes to the wormData structure (ie: adding
        %another column corresponding to a new channel may result
        %in changes to some of these numbers
        WORM_INDEX = 1;
        POSITION_NUM = 2;
        WORM_NUM = 3;
        A594_IND = 4;
        CY5_IND = 5;
        TMR_IND = 6;
        YFP_IND =7;
        NUCLEI_IND = 8;
        TIME_IND = 9;
        ECELL_IND = 10;
        ERRBAR_A594_IND = 1;
        ERRBAR_CY5_IND = 2;
        ERRBAR_TMR_IND =3;

        %Set up a dictionary
        
        channel = struct('a594',4,'cy5',5,'tmr',6,'yfp',7);
        
        %What color should each channel be displayed as?
        colorTable = containers.Map({'a594','cy5','tmr','dapi'},{'y','r','g','b'});
       
    end
    

    
    methods
        
        
        
        function myObj = DuProcessWormData (wormDataStruct)
        %Constuctor method. Requires a wormData as input 
        %Requires DuNucsTime.m
            
            myObj.myWormDataStruct = wormDataStruct;     
     
            
        
            %Sort rows by nuclei count 
            A = myObj.myWormDataStruct.spotNum;
            A = sortrows(A,myObj.NUCLEI_IND);
            
            %%Remove all rows with uncounted nuclei from the
            %%wormDataStruct.spotNum matrix
            B = A(:,myObj.NUCLEI_IND)<0;  %01 matrix with 1's for
                                         %entries with <0
            %Remove rows with nuclei count <0
            A(B,:)=[];
            
            %Remove corresponding rows from errorbar matrices
            myObj.myWormDataStruct.U(B,:)=[];
            myObj.myWormDataStruct.L(B,:)=[];
             
            %Append nuclei->time (in minutes) conversions to 
            %column spotNum matrix 
            %Requires DuNucTime.m
            
            timeVec = DuNucsTime.convertNucsToTime(A(:,myObj.NUCLEI_IND));
            EcellVec = DuNucsTime.convertNucsToEcells(A(:,myObj.NUCLEI_IND));
           
            %Append a nuclei->time column into spotNum matrix 
            A(:,myObj.TIME_IND)=[timeVec];
            A(:,myObj.ECELL_IND) = [EcellVec];
            myObj.myWormDataStruct.spotNum = A;
                        
            %Set the sgolaySpan value
            myObj.sgolaySpan = ...
               myObj.savitskyGolaySpan(myObj.NUCLEI_IND);
           
            %%%%Subfields in myWormDataStruct
            %Add another subfield in myWormDataStruct for holding xyz coords
            myObj.myWormDataStruct.coordList = {};
            %Add another subfield in myWormDataStruct for holding ectopic
            %expression values. Columns are wormIndex, nuclei count,
            %ectopic transcript counts.
            myObj.myWormDataStruct.ectopicExp = [];
           
           
            
        end %Constructor method

        
        function [nucX, countsY] = sortByNuclei(myObj,col_index)
            %This will sort data by nuclei and remove data with
            %uncounted nuclei (flagged -1)
            
            
            %Make matrix with nuclei on left column, counts on right
            A = [ myObj.myWormDataStruct.spotNum(1:end, myObj.NUCLEI_IND),...
                  myObj.myWormDataStruct.spotNum(1:end, col_index)];
            %Sort by first column
            A = sortrows(A,1);
            
            %Delete uncounted nuclei rows (-1) by logical indexing
            %Ref:http://www.mathworks.com/matlabcentral/answers
            %/13476-delete-rows-that-have-a-negative-number-in-their-first-column
            %Note: I altered the class constructor to do this by
            %defaults so this code may be a bit redundant
            
            A( A(:,1)<0,:)=[]

            nucX = A(:,1);
            countsY = A(:,2);
            
        end %sortByNuclei
            
            
            
        function [scatteredPoints] = scatteredPointsFromCol(myObj, yCol);  
        % Extract data points vector from wormDataStructure 
            scatteredPoints = myObj.myWormDataStruct.spotNum(1:end, yCol);
            
        end %end scatteredPointsFromCol
            
            
        function [meanCountsY] = meanTrajectorySavitskyGolayFromCol(myObj, xCol,yCol) 
            
        %(ie: yInd = 4 for a vector of A594 mRNA spotCounts)
            
            %Generate smooth mean trajectory via Savitsky-Golay
            %smoothing
            xData = myObj.myWormDataStruct.spotNum(1:end,xCol); 
            yData = myObj.myWormDataStruct.spotNum(1:end,yCol); 
            meanCountsY = smooth(xData,yData,myObj.sgolaySpan,'sgolay',1)
            
        end %end meanTrajectorySavitskyGolayFromCol          
            
        function [meanY] = meanTrajectorySmoothingSplineFromCol(myObj, ...
                                                              xCol,yCol, ...
                                                                varargin) 
                                                          
            %Set default parameters for optional inputs
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 1
                error('meanTrajectorySmoothingSplineFromCol', ...
                    'requires at most 1 option input: smoothingParameter');
            end
            
            optargs = {0.0001};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            % or ...
            % [optargs{1:numvarargs}] = varargin{:};
            
            % Place optional args in memorable variable names
            [smoothingParam] = optargs{:};   
       
               
        
        %(ie: yInd = 4 for a vector of A594 mRNA spotCounts)
            
            %Generate smooth mean trajectory via Savitsky-Golay
            %smoothing
            xData = myObj.myWormDataStruct.spotNum(1:end,xCol); 
            yData = myObj.myWormDataStruct.spotNum(1:end,yCol); 
            xVals = 1:length(myObj.myWormDataStruct.spotNum(1:end,xCol));
            fitObject = fit(xData,yData,'smoothingspline','SmoothingParam',...
                            smoothingParam);
            meanY = feval(fitObject,xVals);
            
        end %end meanTrajectorySmoothingSplineFromCol    
        
            
       function [meanY] = meanTrajectorySmoothingSpline(myObj, ...
                                                 xData,yData, varargin)
                                                          
                                             
             smoothingParam = 0.0001;
    
            
            %Generate smooth mean trajectory via Savitsky-Golay
            %smoothing
  
            xVals = 1:length(xData);
            fitObject = fit(xData,yData,'smoothingspline','SmoothingParam',...
                            smoothingParam);
            meanY = feval(fitObject,xVals);
            
        end %end meanTrajectorySmoothingSplineFromCol    
        
        
            
       function [variance] = varianceFromSavitskyCol(myObj, xCol,yCol)
           % Subtract mean trajectory from original data and square
           % residuals. Voila, you get variance
            
           %xCol is needed to get meanTrajectory    
           xData = myObj.myWormDataStruct.spotNum(1:end,xCol);
               
           squaredResiduals = (myObj.myWormDataStruct.spotNum(1:end,yCol)- ...
               myObj.meanTrajectorySmoothingSplineFromCol(xCol,yCol,0.0001)  ).^2;
           
           
           
           variance = smooth(xData, squaredResiduals,myObj.sgolaySpan,'sgolay',1);
           %figure(10)
           %scatter(xData,variance) 
           %figure(20)
           %Set values < 0 to 0
           %variance(variance(:,1)<0,:)=[0]
       end %end varianceFromCol     
        
       
       
       
       
       
       function [variance] = varianceFromSplineCol(myObj, xCol,yCol)
           % Subtract mean trajectory from original data and square
           % residuals. Voila, you get variance
            
     
                smoothingParam = 0.0001
     
           
           %xCol is needed to get meanTrajectory    
           xData = myObj.myWormDataStruct.spotNum(1:end,xCol);
               
           squaredResiduals = (myObj.myWormDataStruct.spotNum(1:end,yCol)- ...
               myObj.meanTrajectorySmoothingSplineFromCol(xCol,yCol)  ).^2;
           
           
           
           variance = smooth(xData, squaredResiduals,myObj.sgolaySpan,'sgolay',1);
           

            
            %Generate smooth mean trajectory via Savitsky-Golay
            %smoothing
  
            xVals = 1:length(xData);
            fitObject = fit(xData,squaredResiduals,'smoothingspline','SmoothingParam',...
                            smoothingParam);
            variance = feval(fitObject,xVals);
           
           
           
           
           %figure(10)
           %scatter(xData,variance) 
           %figure(20)
           %Set values < 0 to 0
           %variance(variance(:,1)<0,:)=[0]
       end %end varianceFromCol     
            
       
       
       
       
       
       function [stdDevTrajectory] = stdDevFromCol(myObj, xCol, yCol)
           
           stdDevTrajectory= sqrt(myObj.varianceFromSplineCol(xCol,yCol));
           %scatter (myObj.myWormDataStruct.spotNum(1:end,xCol),stdDevTrajectory)
           
       end %end stdDevFromCol
           
            
           
       function [CV] = coeffVarFromCol(myObj, xCol, yCol)
       %Calculate coefficient of variation for xCol and yCol
           CV = myObj.stdDevFromCol(xCol, yCol)./ ...
               myObj.meanTrajectorySmoothingSplineFromCol(xCol, yCol);
           
       end %coeffVarFromCol
           

       
           
           
      function [sgolaySpan] = savitskyGolaySpan(myObj, dataCol)
          %Calculate 50% span value to use for
           %savitsky-golay based on size of current dataset
               
           dataSpan =.5*length(myObj.myWormDataStruct.spotNum(1:end,dataCol)); 
           %Find nearest odd number
           sgolaySpan= 2*floor(dataSpan/2)+1;
           
       end %savitskyGolaySpan   
           
  
            
        function [] = populateCoordList(myObj)
            %Looks in myPlot.myWormDataStruct for images and worms
            myObj.myWormDataStruct.coordList = {}; %Clear the coordList field in case not clear already
            %Copy columns 1-3 from spotNum
            myObj.myWormDataStruct.coordList = num2cell(myObj.myWormDataStruct.spotNum(:,1:3));
            %Look up newNucallembryos_Pos* and extract data 
            %Cell array columns are WormNum,PosNum,EmbIndex, coordList for
            %dye1 (matrix), coordList for dye2, etc.
            
            
            for i=1:size(myObj.myWormDataStruct.coordList,1);
                %wormNum = myObj.myWormDataStruct.coordList(i,1)
                for j=1:size(myObj.myWormDataStruct.dye,2);
                    
                    posNum = myObj.myWormDataStruct.coordList{i,2};
                    embIndex = myObj.myWormDataStruct.coordList{i,3};
                    myObj.myWormDataStruct.dye(j);
                    
                    goodCoords = DuProcessWormData.loadSpotStatsCoords(posNum,embIndex,myObj.myWormDataStruct.dye(j));
                    %A reminder as to what dye correspondes to what
                    %coordinates!
                    %dyetag = repmat(myObj.myWormDataStruct.dye{j},size(goodCoords,1),1);
                    %goodCoords = [dyetag,goodCoords];
                    myObj.myWormDataStruct.coordList(i,3+j) = {goodCoords};
                end
                
            end
            
       
        end
        

        function [imgh] = imshowChannel(myObj,channelName,wormIndex)
            %Input is a string like 'tmr'
            %Lets you see the max projection of a given channel
            %Must be called in the correct directory containing
            %segmentation and image files
            posNum = myObj.posNumOfWorm(wormIndex);
            posIndex = myObj.posIndexOfWorm(wormIndex);
            
            segStruct = importdata(strcat(channelName,'_Pos',num2str(posNum),'_SegStacks.mat'));
            
            img = segStruct.segStacks{posIndex}; %3dImage
            img = DuImgLib.imAutoScale16to8(max(img(:,:,:),[],3));
            
            imgh= imshow(img);
            hold on;
            
        end
        
        function [imgh] = imshowComposite(myObj,wormIndex)
            subplot(2,2,1)
            imshowChannel(myObj,'tmr',wormIndex);
            title('gfp');
            subplot(2,2,2);
            imshowChannel(myObj,'cy5',wormIndex);
            title('elt-2');
            subplot(2,2,3)
            imshowChannel(myObj,'a594',wormIndex);
            title('end-1');
            subplot(2,2,4);
            reviewImageCoords(myObj, wormIndex);
        end
        
        
        function [dapiImg] = reviewImageCoords(myObj,wormIndex)
            %Lets you visualize maxProjection of image with specified position
            %number with final smFISH spots
            %You must be in the correct dI think Hillary demonstrated she's more likely to appeal to moderatirectory
            
            %%Handle varargins
            
            %numvarargs = length(varargin)
            %if numvarargs > 3
            %    error ('Method requires at most 3 optional inputs',
            %    'dye1, dye2, dye3');
            %end
            
            %optargs = {'a594','cy5','tmr'};
            
            
            %%
                 
            %Determine posNum corresponding to wormIndex
            %Note: I found matlab to be really stupid when handling methods
            % called from other methods within the same class
            % See this post:http://www.mathworks.com/matlabcentral/newsreader/view_thread/322199
            posNum = myObj.posNumOfWorm(wormIndex);
            posIndex = myObj.posIndexOfWorm(wormIndex);
            %Load dapi image
            %dapiImg = loadtiff(strcat('dapi_Pos',num2str(posNum),'.tif'));
            %Get max projection of dapi image (flatten it)
            %dapiImg = DuImgLib.maxZstack(dapiImg);
            
            %imgh = image(dapiImg);
            %colormap(gray(65536)); 
            %hold on;
            
            segStruct = importdata(strcat('dapi_Pos',num2str(posNum),'_SegStacks.mat'));
            dapiImg = segStruct.segStacks{posIndex}; %3dImage
            dapiImg = DuImgLib.imAutoScale16to8(max(dapiImg(:,:,:),[],3));
            
            imshow(dapiImg);
            %colormap(gray)
            %set(gca,'YDir','normal')
            hold on
            
            
            wormRowIndex = myObj.spotNumIndexOfWorm(wormIndex);
            for i=1:size(myObj.myWormDataStruct.dye,2);
                %myObj.myWormDataStruct.dye{i}     %type is char
                
                spotCoords = myObj.myWormDataStruct.coordList{wormRowIndex,3+i};
                Y = spotCoords(:,1); 
                X = spotCoords(:,2);
                
                plot(X,Y,'MarkerSize',5,'LineStyle','none','Marker','.','Color',myObj.colorTable(myObj.myWormDataStruct.dye{i}));
                hold on;
                
            end
        
            
            
            
            
            
        end
        
        
        function [] = reviewEarlySpotLocation(myObj,nucThresholdMin,nucThresholdMax)
           % reviewImageCoords for all images with nuclei counts less than
           % the nuc threshold. This is to help me quantify low level
           % ectopic expression
           nucVector = myObj.myWormDataStruct.spotNum(:,myObj.NUCLEI_IND);
           wormIndexVector =  myObj.myWormDataStruct.spotNum(:,myObj.WORM_INDEX);
           wormPosVector = myObj.myWormDataStruct.spotNum(:,myObj.POSITION_NUM);
           %Retrieve all wormNums corresponding to nuclei counts less than nucThreshold
           logicalMat = ((nucVector >= nucThresholdMin) & (nucVector <= nucThresholdMax))
           thresholdedWormIndices = wormIndexVector(logicalMat)
           thresholdedPosVals = wormPosVector(logicalMat)
           thresholdedNucCounts = nucVector(logicalMat)
           for i=1:length(thresholdedWormIndices)
               disp(strcat('Position:',num2str(thresholdedPosVals(i,1))));
               disp(strcat('Worm Index:',num2str(thresholdedWormIndices(i,1))));
               disp(strcat('Nuclei count:',num2str(thresholdedNucCounts(i,1))));
               disp(' ');
               
               imshowComposite(myObj,thresholdedWormIndices(i,1));
               waitforbuttonpress;
               clf;
           end
            
            
        end
        
        function [] = imshowCompositeWithAlphaVol(myObj,channelName,wormIndex)
            dapiImg = myObj.reviewImageCoords(wormIndex);
            hold on;
            D = size(dapiImg);
            myObj.plotAlphaShapePoly(channelName,wormIndex,D(1),D(2));
            hold off;
        end
        
       
        function [polyMask] = plotAlphaShapeMask(myObj,channelName,wormIndex,varargin)
    
            radius = 75;
            %%%%%%%%%%%%%%%
          
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Set default parameters for optional inputs
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 2
                error('plotAlpha', ...
                    'requires at most 2 option input: imgW,imgH');
            end
            
            optargs = {500, 500};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            % or ...
            % [optargs{1:numvarargs}] = varargin{:};
            
            % Place optional args in memorable variable names
            [imgW,imgH] = optargs{:}   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure;
            [V,S] = DuAlphaVol(myObj.retrieveXy(channelName,wormIndex),radius,1)
            %x and y are reversed for this data set
            %handle = patch(S.x,S.y,myObj.colorTable(channelName));
            %axis([0 imgW 0 imgH]);
            polyMask = poly2mask(S.x,S.y,imgW,imgH);
            polyMask = imcomplement(polyMask);
            imshow(polyMask)
            alpha(.04);
    
        end
        
        function [patchHandle] = plotAlphaShapePoly(myObj,channelName,wormIndex,varargin)
    
            radius = 75;
            %%%%%%%%%%%%%%%
          
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Set default parameters for optional inputs
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 2
                error('plotAlpha', ...
                    'requires at most 2 option input: imgW,imgH');
            end
            
            optargs = {500, 500};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            % or ...
            % [optargs{1:numvarargs}] = varargin{:};
            
            % Place optional args in memorable variable names
            [imgW,imgH] = optargs{:}   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            [V,S] = DuAlphaVol(myObj.retrieveXy(channelName,wormIndex),radius)
            %x and y are reversed for this data set
            patchHandle = patch(S.x,S.y,myObj.colorTable(channelName),'FaceAlpha',0.3);
       
            hold off;
         
            
        end
              
        function [ectopicCount,alphaVolX,alphaVolY] = numPointsOutsideAlphaVol( myObj, alphaVolChannel,ectopicChannel, wormIndex, varargin)
            %Used for quantifying ectopic expression.
            %Requires script DuAlphaVol
            %Take the points of the alphaVolChannel and construct a 2D alpha volume
            %Then determine how many points of the ectopic channel 
            
            radius = 75;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Set default parameters for optional inputs (radius
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 2
                error('numPointsOutsideAlphaVol', ...
                    'requires at most 1 option input:radius, boolPlotFig');
            end
            
            optargs = {radius,0};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            % or ...
            % [optargs{1:numvarargs}] = varargin{:};
            
            % Place optional args in memorable variable names
            [radius,boolPlotFig] = optargs{:};   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            alphaXy = myObj.retrieveXy(alphaVolChannel, wormIndex);
            ectopicXy = myObj.retrieveXy(ectopicChannel, wormIndex)
            [V,S]= DuAlphaVol(alphaXy,radius);
            alphaVolX = S.x;
            alphaVolY = S.y;
            ectopic = ~inpolygon(ectopicXy(:,1),ectopicXy(:,2),S.x,S.y)
            ectopicCount = sum(ectopic);
            
            if boolPlotFig
                figure;
                myObj.reviewImageCoords(wormIndex);
                patch(S.x,S.y,'b','FaceAlpha',0.3);
            end
            
            
        end
        
        
        
        function [ectopicCount,alphaVolX,alphaVolY] = numPointsOutsideAlphaVol2( myObj, alphaVolChannel1,alphaVolChannel2,ectopicChannel, wormIndex, varargin)
            %Used for quantifying ectopic expression.
            %Requires script DuAlphaVol
            %Take the points of two specified channels and construct a 2D alpha volume
            %Then determine how many points of the ectopic channel lie
            %outside the volume
            
            radius = 75;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Set default parameters for optional inputs (radius
            %Good instructions at
            %http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
            numvarargs = length(varargin);
            if numvarargs > 2
                error('numPointsOutsideAlphaVol', ...
                    'requires at most 2 option input:radius,boolPlotFig');
            end
            
            optargs = {radius,0};
            
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            % or ...
            % [optargs{1:numvarargs}] = varargin{:};
            
            % Place optional args in memorable variable names
            [radius,boolPlotFig] = optargs{:};   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            alphaXy1 = myObj.retrieveXy(alphaVolChannel1, wormIndex);
            alphaXy2 = myObj.retrieveXy(alphaVolChannel2, wormIndex);
            alphaXyComposite = [alphaXy1;alphaXy2];
            ectopicXy = myObj.retrieveXy(ectopicChannel, wormIndex);
            [V,S]= DuAlphaVol(alphaXyComposite,radius);
            alphaVolX = S.x;
            alphaVolY = S.y;
            ectopic = ~inpolygon(ectopicXy(:,1),ectopicXy(:,2),S.x,S.y);
            ectopicCount = sum(ectopic);
            
            
            %Plot composite image with alpha
            if boolPlotFig
                figure;
                myObj.reviewImageCoords(wormIndex);
                patch(S.x,S.y,'b','FaceAlpha',0.3);
            end
            
        end
        
        
        function [ectopicMat] = quantifyEctopicExp(myObj)
            %Determines ectopic GFP/TMR counts for every embryo
            %Produces a mx3 matrix with the first column consisting of
            %WormIndexes, the second column consisting of nuclei counts,
            %and the third column consisting of ectopic expression
            %estimates.
             %This currently uses 'a594' and 'cy5' channels for making alpha 
           %shape, and 'tmr' for making spots.
            
            
            numSkipped =0;
            numEmbryos =size(myObj.myWormDataStruct.coordList(:,1),1)
            ectopicMat = [];
            curEctopic = [0,0,0]; %init value
            for i=1:numEmbryos
                try
                   
                   curWormIndex = cell2mat(myObj.myWormDataStruct.coordList(i,myObj.WORM_INDEX))
                   curEctopic = myObj.numPointsOutsideAlphaVol2('a594','cy5','tmr',curWormIndex,75,0);%Radius 75. Do not show figure
                   curNuc = myObj.nucCountOfWorm(curWormIndex)
                   ectopicMat = cat(1,ectopicMat,[curWormIndex,curNuc,curEctopic]);
                    %ectopicMat = cat(1,ectopicMat,curEctopic);
                catch 
                    disp(strcat('WormIndex ',num2str(i),' missing data. Skipping...'));
                    numSkipped = numSkipped+1;
                end
            end
            
            disp(strcat('Number of embryos sampled: ',num2str(numEmbryos-numSkipped) ) );
            
            %Store this matrix in the object
            myObj.myWormDataStruct.ectopicExp = ectopicMat;
        end
        
                
        function [ectopicMat] = quantifyEctopicElt2Exp(myObj)
            %Determines ectopic GFP/TMR counts for every embryo
            %Produces a mx3 matrix with the first column consisting of
            %WormIndexes, the second column consisting of nuclei counts,
            %and the third column consisting of ectopic elt-2 expression
            %estimates
             %This currently uses 'a594' and 'cy5' channels for making alpha 
           %shape, and 'tmr' for making spots.
            
            
            numSkipped =0;
            numEmbryos =size(myObj.myWormDataStruct.coordList(:,1),1)
            ectopicMat = [];
            curEctopic = [0,0,0]; %init value
            for i=1:numEmbryos
                try
                   
                   curWormIndex = cell2mat(myObj.myWormDataStruct.coordList(i,myObj.WORM_INDEX))
                   curEctopic = myObj.numPointsOutsideAlphaVol('a594','cy5',curWormIndex,75,0);%Radius 75. Do not show figure
                   curNuc = myObj.nucCountOfWorm(curWormIndex)
                   ectopicMat = cat(1,ectopicMat,[curWormIndex,curNuc,curEctopic]);
                    %ectopicMat = cat(1,ectopicMat,curEctopic);
                catch 
                    disp(strcat('WormIndex ',num2str(i),' missing data. Skipping...'));
                    numSkipped = numSkipped+1;
                end
            end
            
            disp(strcat('Number of embryos sampled: ',num2str(numEmbryos-numSkipped) ) );
            ectopicMat
            %Store this matrix in the object
            
        end
        
        
        
        function [nucCount] = nucCountOfWorm(myObj,wormIndex)
            A = myObj.myWormDataStruct.spotNum(:,myObj.WORM_INDEX)==wormIndex;
            B = myObj.myWormDataStruct.spotNum(:,myObj.NUCLEI_IND);
            nucCount = B(A);
        end
        
        
        function [index] = spotNumIndexOfWorm(myObj,wormIndex)
            A = myObj.myWormDataStruct.spotNum(:,myObj.WORM_INDEX)==wormIndex; %Logical matrix
            B = find(A); %Identify nonzero elements
            index = B(1); %Select first nonzero element index
        end
        
          
        
        function [posNum] = posNumOfWorm(myObj,wormIndex)
               A = myObj.myWormDataStruct.spotNum(:,myObj.WORM_INDEX)==wormIndex; %Logical matrix
               B = myObj.myWormDataStruct.spotNum(:,myObj.POSITION_NUM); %All the Position numbers
               posNum = B(A);
        end
        
        function [posIndex] = posIndexOfWorm(myObj,wormIndex)
            A = myObj.myWormDataStruct.spotNum(:,myObj.WORM_INDEX)==wormIndex; %Logical matrix
            B = myObj.myWormDataStruct.spotNum(:,myObj.WORM_NUM); %All the Position numbers
            posIndex = B(A);
        end
        
        
        function [matXyz] = retrieveXyz(myObj,channelName,wormIndex)
            %Retrieves the n x 3 xyz coordinate matrix of a given wormIndex
            rowNum = myObj.spotNumIndexOfWorm(wormIndex);
            matXyz = myObj.myWormDataStruct.coordList{rowNum,myObj.channel.(channelName)};
            
            %The coordinate order is actually y,x,z in Scott's original
            %coordinate list. This line rearranges things to x,y,z
            matXyz = [matXyz(:,2),matXyz(:,1),matXyz(:,3)];
        end
        
        function [matXy] = retrieveXy(myObj,channelName,wormIndex)
            %Retrieves the n x 2 xy coordinate matrix of a given wormIndex
            xyz = myObj.retrieveXyz(channelName,wormIndex);
            matXy = xyz(:,1:2);
        end
        

        
       function [] = blacklistWorm(myObj, embryoIndex)
           %BLACKLISTWORM
           %Removes row from wormDataStruct.spotNum matrix with
           %corresponding embryo index.
           % Usage: obj.blacklistWorm(75)
           % (will blacklist worm with WORM_IND 75)    
           A=myObj.myWormDataStruct.spotNum;
           EU=myObj.myWormDataStruct.U;
           EL = myObj.myWormDataStruct.L;
           C = myObj.myWormDataStruct.coordList;
           B = A(:,myObj.WORM_INDEX)==embryoIndex; %01 matrix with 1's
           %for entries
           %matching to embryoIndex
           
           %Append row to blacklist
           myObj.blacklisted = [myObj.blacklisted; A(B,:)];
           A(B,:)=[];
           EU(B,:)=[];
           EL(B,:)=[];
           C(B,:)=[];
           myObj.myWormDataStruct.spotNum = A;
           myObj.myWormDataStruct.U = EU;
           myObj.myWormDataStruct.L = EL;
           myObj.myWormDataStruct.coordList = C;
       end%blacklistEmbryo
       
       
       
       
       
       function [] = addNote(myObj)
           myObj.notes{end+1,1} = input('Enter your note: ','s');
       end %addNote
   
       end %end methods
       
       
       
       
       methods (Static = true)
        
           function [goodCoords] = loadSpotStatsCoords(posNum,embIndex,dyeName)
            %Loads xyz coordinates from fluor_PosNum_spotStats.mat
            %Only loads the good coordinates
            
            %Go down the list of non blacklisted worms in myObj.myWormDataStruct.coordList
            % -Open 'dyeName'_Pos'posNum'_spotStats.mat file and load cell
            % 'embIndex'
            % -Extract xyz data from spotStats.locAndClass where 4th column has 1 (for good spot)
            
            filename = strcat(dyeName,'_Pos',num2str(posNum),'_spotStats.mat') %This is a cell
            filename = filename{1} %Get string from cell contents
            spotStats = importdata(filename) %Do not use load
                       
            goodCoords = spotStats{embIndex}.locAndClass(spotStats{embIndex}.locAndClass(:,4)==1,1:3);
           end
           
           
   
           
           function [yTranscripts] = extractRangeInclusive(matIn,xCol,yCol,xMin,xMax)
               %I NO LONGER HAVE USE FOR THIS FUNCTION
               
               %Used for extracting xIndexed data from a matrix.
               %For instance, in a matrix where the second column
               %represents nuclei counts, and the third column represents
               %transcripts, extract transcript counts for values between 8
               %and 60 nuclei (inclusive)
               xVals = matIn(:,xCol)
               yVals = matIn(:,yCol)
               
               xLogical = (xVals>=xMin) & (xVals<=xMax)
               
               yTranscripts = yVals(xLogical);
           end
        
           
           
           function [percentile] = simplePermutation(data1,data2,nReps)
               %Perform a permutation test on two sets of numbers
               
               
               meanDiff = mean(data1)-mean(data2);
                          
               %Preallocate array
               randMeanDiffs = zeros(1,nReps);
               
               for i=1:nReps
                   randMeanDiffs(1,i)= permMeanDiff(data1,data2);
               end
               
               
               percentile = sum(randMeanDiffs<=meanDiff)/nReps     
               
               
               %%Nested functions
               function [randMeanDiff] = permMeanDiff(data1,data2)
                   %Scramble your input data into two scrambled data sets,
                   %and find the mean difference of the two scrambled data sets. 
                   
                   %Concatenate data
                   allData = [data1,data2];
                   %Randomize values 
                   randMat = allData(randperm(numel(allData)));
                   
                   data1Length = size(data1,2);
                   %Slice randomized data into groups with the same sizes
                   %as the original data.
                   randData1 = randMat(1:data1Length);
                   randData2 = randMat(data1Length+1:end);
                   
                   randMeanDiff = mean(randData1)-mean(randData2);
               end
               
           end
           
           function [stagePercentileMat] = binnedEctopicObjPermTest(obj1,obj2,nReps)
               %Take two objects that have
               %myObj.myWormDataStruct.ectopicExp filled in,
               %and do a permutation test binned by E cell number/stage
               obj1.quantifyEctopicExp();
               obj2.quantifyEctopicExp();
               stagePercentileMat = DuProcessWormData.binnedEctopicPermTest(obj1.myWormDataStruct.ectopicExp,obj2.myWormDataStruct.ectopicExp,nReps)
           end
           
           
           function [stagePercentileMat] = binnedEctopicPermTest(ectopicExp1,ectopicExp2,nReps)
       
               %Run permutation tests for each Ecell stage. 7,14,44,100 nucs match
               %to 1E,2E,4E, and 8E respectively.
               
               %Requires DuNucsTime.m
               
               %Stages to perform the permutation test on 
               stages = [1,2,4,8,16]';
                          
               
               ecells1 = DuNucsTime.convertNucsToEcells(ectopicExp1(:,2));
               ecells2 = DuNucsTime.convertNucsToEcells(ectopicExp2(:,2));
               
               %Append ecells annotations to fourth column of ectopic matricies
               
               ectopic1 = [ectopicExp1,ecells1];
               ectopic2 = [ectopicExp2,ecells2];
               
               stagePercentileMat = [stages,zeros(length(stages),1)];
               
               for i=1:length(stages)
                   %Extract rows corresponding to each stage
                   
                   curStage1 = ectopic1((stages(i,1) == ectopic1(:,4)),:);
                   curStage2 = ectopic2((stages(i,1) == ectopic2(:,4)),:);
                  
                   curData1 = curStage1(:,3)';
                   curData2 = curStage2(:,3)';
                  
                  
                  stagePercentileMat(i,2) = DuProcessWormData.simplePermutation(curData1,curData2,nReps);
                  
               end
               
           end
           
       end %End Static methods
       
       


       
    end % end classdef
    
    
