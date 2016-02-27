function [axh] = compareTwoDataClouds(DuPlotData1,...
    channel1,...
    DuPlotData2, ...
    channel2, ...
    varargin)
%Permutation test written by Scott Rifkin
%This will perform a permutation and plot out the results.
%I modified this to take DuPlotWormData objects as
%input --Larry
% Usage:
% SomeObject.compareTwoDataClouds(wormDataObj1,'tmr',wormDataObj2,'tmr')
% SomeObject.compareTwoDataClouds(wormDataObj1,'tmr',wormDataObj2,'tmr')
% Additionally, user can add xmin,xmax,smoothingParam as
% arguments.


%Set default parameters for optional inputs
%Good instructions at
%http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
numvarargs = length(varargin);

if numvarargs > 4
    error('compareTwoDataClouds', ...
        'requires at most 3 optional inputs (xmin, xmax,',...
        'smoothingParam, rRepetitions)');
end

%Default arguments
% >>>>>>>>>>>>>>>>>>>>>>>>
% smoothing parameter = 0.0001; numrepetitions = 2000
optargs = {0,160,0.0001,10000};
% <<<<<<<<<<<<<<<<<<<<<<<<

%44 nuclei corresponds to start of 4E stage
% now put these defaults into the valuesToUse cell array,
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;
% or ...
% [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names
[xmin, xmax, smoothingParam,nRepetitions] = optargs{:}

channelIndex1 = DuPlotData1.channel.(channel1);
channelIndex2 = DuPlotData2.channel.(channel2);

xRange = [xmin, xmax];

%Plot options:
colorRandom = [.65,.84,.52]; %light guacamole
colorMeanDifference = [0,0.45,0]; %green
colorP = 'blue';
colorEcell = [0.7 0.7 0.7]; %lightgray

% x1 is nuclei, y1 is spots.  Same with x2 and y2.
x1 = DuPlotData1.myWormDataStruct.spotNum(:,DuPlotData1.NUCLEI_IND);
x2= DuPlotData2.myWormDataStruct.spotNum(:,DuPlotData2.NUCLEI_IND);


y1 = DuPlotData1.myWormDataStruct.spotNum(:, channelIndex1);
y2= DuPlotData2.myWormDataStruct.spotNum(:, channelIndex2);




%Fit smoothing spline to Data1 channel and Data2 channel
%>>>>>>>>>>>>>>>>>>>>>
disp(smoothingParam)
[curve1,~,output1]=fit(x1,y1,'smoothingspline','SmoothingParam',smoothingParam);
[curve2,~,output2]=fit(x2,y2,'smoothingspline','SmoothingParam',smoothingParam);
%<<<<<<<<<<<<<<<<<<<<<<

%Take the difference of mean curves at each x-value
xValsToUse=[xRange(1):xRange(2)];
mean1=feval(curve1,xValsToUse);
mean2=feval(curve2,xValsToUse);
diff=mean1-mean2;

%%First plot of mean values
%closeFigNum(6);
%figure(6)
%subplot(2,1,1)
%plot(xValsToUse,mean1,'c','LineWidth',2);
%hold on
%plot(xValsToUse,mean2,'m','LineWidth',2);
%plot(x1,y1,'c.');plot(x2,y2,'m.');
%hold off
%xlabel('# nuclei');
%ylabel('# transcripts');
%xlim(xRange);

%Now repeat this process except scramble up the data
%Note that this assumes that all the data in your vectors is good.  If it
%isn't good remove it before passing it into the function
len1=length(x1);
len2=length(x2);
nAllData=len1+len2;

allX=[x1;x2];
allY=[y1;y2];

%nRepetitions=500;

%Matrix to hold the random permutations
randomDifferences=zeros(length(diff),nRepetitions); 
%  >>>>>>>>>>>>>>>>>>>>

% Determine the number of comparisons between the curves based on the number of parameters in the smoothing splines
% Then use the Dunn-Sidak method (Ury 1976) to find an adjusted alpha value
%output1 and outpu2 are from the fit() calls above
nparam1=output1.numparam;
nparam2=output2.numparam;
experimentWiseAlpha=0.05;
DSAlpha=1-(1-experimentWiseAlpha)^(1/max(nparam1,nparam2));



% E cell stage annotations from addECellAnnotationsNucs
% 7, 14, 44, 100
stages={0:13; 14:43; 44:99; 100:300};%300 as sufficiently late endpoint
% want to break the data up by stage so that the permutation of strain labels
% only happens within a stage (except group 0E and 1E) 0E is 0-6, 1E is 7-13 nuclei

allX_byStage={};
allY_byStage={};
nData_byStage=[];
n1Data_byStage=[];
n2Data_byStage=[];

for iStage=1:length(stages)
    nData_byStage(iStage)=sum(ismember(allX,stages{iStage}));%how many datapoints are there in this stage?
    n1Data_byStage(iStage)=sum(ismember(x1,stages{iStage}));
    n2Data_byStage(iStage)=sum(ismember(x2,stages{iStage}));
    allX_byStage{iStage}=reshape(allX(ismember(allX,stages{iStage})),1,nData_byStage(iStage));%ensure it is row vector
    allY_byStage{iStage}=reshape(allY(ismember(allX,stages{iStage})),1,nData_byStage(iStage));%ensure it is row vector
end;

parfor i=1:nRepetitions
    if mod(i,500)==1
        disp(i);
    end;
    randX1=[];
    randX2=[];
    randY1=[];
    randY2=[];
    %permute within a stage
    for iStage=1:length(stages)
        p=randperm(nData_byStage(iStage));
        randX=allX_byStage{iStage}(p);
        randY=allY_byStage{iStage}(p);
        randX1=[randX1 randX(1:n1Data_byStage(iStage))];
        randX2=[randX2 randX((n1Data_byStage(iStage)+1):end)];
        randY1=[randY1 randY(1:n1Data_byStage(iStage))];
        randY2=[randY2 randY((n1Data_byStage(iStage)+1):end)];
    end;
        %<<<<<<<<<<<<<<<<<<<<<<<<

    randCurve1=fit(randX1',randY1','smoothingspline','SmoothingParam',smoothingParam);
    randCurve2=fit(randX2',randY2','smoothingspline','SmoothingParam',smoothingParam);
    
    randMean1=feval(randCurve1,xValsToUse);
    randMean2=feval(randCurve2,xValsToUse);
    randDiff=randMean1-randMean2;
    %%Add randDiff to column in random differences matrix
    randomDifferences(:,i)=randDiff;
end;

%%Overlay percentiles on plot
%Plot percentiles

%Find median value at each X(row in randomDifferece) 
md=median(randomDifferences,2);
dpmd=diff-md;



mdRep=repmat(md,1,nRepetitions);
compareMat=randomDifferences-mdRep;
comparisons=compareMat<=repmat(dpmd,1,nRepetitions);
dp=sum(comparisons,2)/nRepetitions;

percentile2Sided=2*(.5-abs(dp-.5));
log10Percentile2Sided=log10(percentile2Sided);





sortedRandomDifferences=randomDifferences;  %size= length(diff),nRepetitions
for i=1:size(sortedRandomDifferences,1)
    sortedRandomDifferences(i,:)=sort(sortedRandomDifferences(i,:));
end;
    
lowerI=round(nRepetitions*DSAlpha/2);
upperI=round(nRepetitions*(1-DSAlpha/2));

%Set up plots
% % % figHandle= figure('Position' ,[0,0,500,400]);
% % % %First axesHandle will plot 

% % % set(axesHandle1,'YAxisLocation','right','YColor',colorMeanDifference,'xlim',[0 160], 'ylim',[-200 200]);
% % % %axis(axesHandle1,[min(xValsToUse) max(xValsToUse) -200 200]);
% % % ylabel(axesHandle1,'Difference in mean # of transcripts','Color',...
% % %                          colorMeanDifference,'FontSize',14,'FontWeight','bold','Fontname','Arial');
% % % xlabel(axesHandle1,'Number of nuclei','FontSize',14,'FontWeight','bold','Fontname','Arial');
% % % 
% % % axesHandle1Pos = get(axesHandle1,'Position');
% % % 
% % % hold(axesHandle1,'on');
% % % %Note, you must set color of axesHandle2 to 'none' for transparency
% % % %The x axis was not lining up.
% % % axesHandle2 = axes('Parent',figHandle, 'Position',axesHandle1Pos, ...
% % %     'YAxisLocation','left','Color','none','YColor',colorP,'xlim',[0 160], 'ylim',[0 1]);
% % % ylabel(axesHandle2,'Percentile','Color',colorP,...
% % %                                      'FontSize',14,'FontWeight','bold' ,'Fontname','Arial');
% % %                                  linkaxes([axesHandle1, axesHandle2],'x');
% % %                                  set(axesHandle2,'XLim',[0 160]);
%xlabel(axesHandle2,'Number of nuclei','FontSize',14,'FontWeight','bold','Fontname','Arial');

% Percentile plot goes from 0 to 1.  Just scale things



%>>>>>>>>>>>>>>>>>>

%Plotting order: (back to front)
%E cell annotations
%significance cutoff lines
%blue percentile line
%null distribution envelope (if plot the full null distribution, then this should be in back of the blue percentile line)
%zero line
%actual difference line


%Add E cell annotations to plot
% Do this here so that they are in the background


nx=length(xValsToUse);
%ensure column vectors

lowerAlphaLine=reshape(DSAlpha/2   +zeros(size(xValsToUse)),nx,1);
upperAlphaLine=reshape(zeros(size(xValsToUse))+     (1-DSAlpha/2),nx,1);
adjustedAlphaLine=reshape(log10(DSAlpha)+zeros(size(xValsToUse)),nx,1);
pvalLine=reshape(log10Percentile2Sided,nx,1);
pvalLine(isinf(pvalLine))=-log10(nRepetitions);

lowerNullLine=reshape(sortedRandomDifferences(:,lowerI),nx,1);
upperNullLine=reshape(sortedRandomDifferences(:,upperI),nx,1);

zeroLine1=reshape(zeros(size(xValsToUse)),nx,1);
zeroLine2=reshape(md,nx,1);

diffLine=reshape(diff,nx,1);

nTranscriptMax=50*ceil(max(abs([abs(lowerNullLine);upperNullLine]))/50);%go to the next 50 above

x=reshape(xValsToUse,nx,1);
figure;
[axh,lh1, lh2]=plotyy([x x],[adjustedAlphaLine pvalLine], [x x x x ], [zeroLine1 lowerNullLine upperNullLine diffLine ]);

set(axh(1),'YLim',[-log10(nRepetitions) 0],'YTick',(-log10(nRepetitions)):0,'FontName','Arial','FontSize',14,'YTickLabel',{'<0.0001','0.001','0.01','0.1','1'});
set(axh(2),'YLim',[-nTranscriptMax nTranscriptMax],'YTick',-nTranscriptMax:100:nTranscriptMax,'FontName','Arial','FontSize',14,'Xtick',[],'XTickmode','manual');
set(lh1(1),'Color',colorP,'LineStyle',':','LineWidth',1);
set(lh1(2),'Color',colorP,'LineStyle','--','LineWidth',1);
%set(lh1(3),'Color',colorP,'LineStyle','-','LineWidth',1);

set(lh2(1),'Color',colorMeanDifference,'LineStyle','-','LineWidth',1);
set(lh2(2),'Color',colorRandom,'LineStyle','-','LineWidth',1);
set(lh2(3),'Color',colorRandom,'LineStyle','-','LineWidth',1);
set(lh2(4),'Color',colorMeanDifference,'LineStyle','-','LineWidth',3);
% % line([0;0], [upperNullLine(1); lowerNullLine(1)],'Color',colorRandom,'LineStyle','-','LineWidth',1,'Parent',axh(2));
% % line([160; 160], [upperNullLine(end); lowerNullLine(end)],'Color',colorRandom,'LineStyle','-','LineWidth',1,'Parent',axh(2));
% % %find parts the are above 200
% % above200=upperNullLine>=200;
% % belown200=lowerNullLine<=-200;
% % boundaryAbove=find(above200(2:end)-above200(1:(end-1)));
% % disp(boundaryAbove);
% % if mod(length(boundaryAbove),2)==1
% %     boundaryAbove(end+1)=160;
% % end;
% % disp(boundaryAbove);
% % for i=1:2:length(boundaryAbove)
% %     line([boundaryAbove(i) boundaryAbove(i+1)],[200 200],'Color',colorRandom,'LineStyle','-','LineWidth',1,'Parent',axh(2));
% % end;
% % 
% % boundaryBelow=find(belown200(2:end)-belown200(1:(end-1)));
% % if mod(length(boundaryBelow),2)==1
% %     boundaryBelow(end+1)=160;
% % end;
% % for i=1:2:length(boundaryBelow)
% %     line([boundaryBelow(i) boundaryBelow(i+1)],[-200 -200],'Color',colorRandom,'LineStyle','-','LineWidth',1,'Parent',axh(2));
% % end;


xlabel(axh(1), 'Number of nuclei');
ylabel(axh(1),'p-value','Color',colorP);
ylabel(axh(2),'Difference in mean # of transcripts','Color',colorMeanDifference);


DuPlotWormData.addECellAnnotationsNucs (axh(1),1,1);


% % % plot(x,lowerNullLine,'Color',colorRandom,'LineStyle','-','LineWidth',1);
% % % axh=gca;
% % % set(axh,'YLim',[-200 200],'YTick',-200:50:200,'FontName','Arial','FontSize',14);
% % % hold on
% % % lh2=line(x,upperNullLine,'Color',colorRandom,'LineStyle','-','LineWidth',1);
% % % lh3=line(x,zeroLine1,'Color',colorMeanDifference,'LineStyle','-','LineWidth',1);
% % % lh4=line(x,diffLine,'Color',colorMeanDifference,'LineStyle','-','LineWidth',3);
% % % xlabel(axh, 'Number of nuclei');
% % % ylabel(axh,'Difference in mean # of transcripts');
% % % DuPlotWormData.addECellAnnotationsNucs (axh,1,1);


% %<<<<<<<<<<<<<<<<<<<<<<
% %Plot alpha annotation line
% % >>>>>>>>>>>>>>>>>>>>>>>>
% line(xValsToUse,DSAlpha/2+zeros(size(xValsToUse)),...
%     'LineWidth',1,...
%      'Color',colorP,...
%      'Parent',axesHandle2,...
%  'LineStyle','-'); %Alpha line set at DSAlpha/2
% 
% line(xValsToUse,ones(size(xValsToUse))-DSAlpha/2,...
%     'LineWidth',1,...
%      'Color',colorP,...
%      'Parent',axesHandle2,...
%  'LineStyle','-'); %Alpha line set at 1-DSAlpha/2
%  % <<<<<<<<<<<<<<<<<<<<<<<<
% %%Plot P-value line                              
% line(xValsToUse,dp,'Color',colorP,...
%     'LineWidth',1,...
%     'Parent',axesHandle2); %Zero line
% 
% 
% 
% %%Second plot
% %subplot(2,1,2)
% 
% %Plot mean (spline) differences for randomized label datasets
% %>>>>>>>>>>
% %change is don't plot all the random splines...just plot the significance cutoff envelope
% %diffLine = line(xValsToUse,randomDifferences,'Color',colorRandom,'Parent',axesHandle1);
% lowerLine=line(xValsToUse,sortedRandomDifferences(:,lowerI),'Color',colorRandom,'Parent',axesHandle1);
% upperLine=line(xValsToUse,sortedRandomDifferences(:,upperI),'Color',colorRandom,'Parent',axesHandle1);
% %<<<<<<<<<<<<<<<<<<<<<<
% %Plot zero line 
% line(xValsToUse,zeros(size(xValsToUse)),'Parent',axesHandle1,...
%     'Color',colorMeanDifference,'LineWidth',1); %Zero line for mean comparison
% 
% %Plot actual mean (spline) differences for both datasets 
% %(diff is mean1-mean2)
% line(xValsToUse,diff,'Color',colorMeanDifference,'LineWidth',3,'Parent',axesHandle1);
% 
% %axesHandle1Pos = axesHandle1.Position









%axesHandle2 = axes('Parent',figHandle,...
%    'YAxisLocation','left','Color','none','YColor',colorP);








% %Figure annotations here
% %Beware, this gets VERY HAIRY
% %%%%%%%%%Set global positions here:
% 
% yPosE = 0.98;
% w=0.005;
% h=0.005;
% yPosAlpha = 0.08;
% xPosAlpha = 0.80;
% 
% yPosAlpha2 = 1.02;
% 
% %Collect axesHandle position [x,y,w,h] in normalized figure units
% xAxOff = axesHandle1Pos(1);
% yAxOff = axesHandle1Pos(2);
% widthAx= axesHandle1Pos(3);
% heightAx = axesHandle1Pos(4)
% %Alpha annotation
% %Normalize user specified annotation positions
% xPosAlpha = xAxOff+xPosAlpha*widthAx;
% 
% yPosAlpha = yAxOff+yPosAlpha*heightAx
% yPosAlpha2 = yAxOff+yPosAlpha2*heightAx
% %>>>>>>>>>>>>>>>>>>>>>>>
% anAlpha = annotation(figHandle,'textbox',[xPosAlpha,yPosAlpha,w,h],...
%            'String', num2str(DSAlpha/2,'%0.04f'),'LineStyle','none',...
%            'FontSize',10,'Color',colorP);
% 
% 
% anAlphaTop = annotation(figHandle,'textbox',[xPosAlpha,yPosAlpha2,w,h],...
%            'String',  num2str(1-DSAlpha/2,'%0.04f'),'LineStyle','none',...
%            'FontSize',10,'Color',colorP);
%      
% 
% % Moved ECell annotations to first so that those lines are behind everything else
%        %<<<<<<<<<<<<<<<<<<<<<<<<  
% 
% %%Eliminate x-ticks and labels on axesHandle1
% %Warning, this can result incorrect axes alignment upon figure rescaling.
% %Do not rescale figure if this is on.
% set(axesHandle1, 'Xtick',[]);

%&export_fig('filename',[DuPlotData1.strainName '_' DuPlotData2.strainName '.pdf'],'-transparent','-pdf');

end

