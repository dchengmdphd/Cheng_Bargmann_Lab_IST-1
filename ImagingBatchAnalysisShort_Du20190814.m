%read experiment info from xls file and populate ExpInfo structure
%Du Version
%---xls file structure---
%gtRep:         genotype replicate - refers to a set of animals of the same
%               genotype from the same experiment
%exp:
%experimentor:
%expFullName:
%gt:
%gtAn:
%PlotGroup:
%pos:
%folder:
%fullPath:
%Function to change stimulation shading is doseRespTim.m

set(0,'defaultfigurecolor',[1 1 1])
%   clear
reReadTxtFiles=1;
YLimTraces=[-.30 5];



plotNC_yn=1;
plotTraces_yn=1;
plotIndDR_yn=1;

[f,p]=uigetfile('E:\trackingArena\ImagingExperimentList_Du.xls','Select .XLS file containing experiment data:');
if ~isempty(f)
    xlsfilename = [p,f];
    ExpInfo = xls2struct(xlsfilename);
else
    error('No file selected');
end


%READ DATA
%collect experiments selected for plotting (PlotGroup>0 in XLS file)
%read to AllSqIntAll and srtAll cells
%shift traces by adding nan if necessary and assemble AllSqInt matrix and
%srt structure

bool_gtRepToPlot=[ExpInfo.PlotGroup]>0;
gtRepToPlot=find(bool_gtRepToPlot);

for i_gtRep=gtRepToPlot;
    disp(['loading data for ',ExpInfo(i_gtRep).expFullName])
    [tmp_AllSqIntAll,tmp_srt]=loadImagingDataBatch({ExpInfo(i_gtRep).fullTracksPath},[],ExpInfo(i_gtRep),reReadTxtFiles);
   % tmp_srt.StimPat=repmat(ExpInfo(i_gtRep).stim,size(tmp_srt.AnimalNr),1);
    ExpInfo(i_gtRep).numFrames=size(tmp_AllSqIntAll,1);
    
    %delete animals not part of current gtRep to save memory
    
    keep=iselement(tmp_srt.AnimalNr,ExpInfo(i_gtRep).gtAn);
    AllSqIntAll{i_gtRep}=tmp_AllSqIntAll(:,keep);
    srtAll{i_gtRep}=trim_strc(tmp_srt,~keep);
    
    disp(['done loading data'])
end

%assemble AllSqInt matrix and srt structure

% closefigs
unq_PlotGroup=unique([ExpInfo.PlotGroup]);
AllSqInt=[];
srt=[];
AnimalNrCont=[];
cols=sum(unq_PlotGroup>0);

for poolGroup=find(unq_PlotGroup>0)
    for i_gtRep=find([ExpInfo.PlotGroup]==unq_PlotGroup(poolGroup))
        %add nan to beginning of traces to shift if offset specified in XLS
        if ExpInfo(i_gtRep).offset>0
            AllSqIntAll{i_gtRep}=[nan(ExpInfo(i_gtRep).offset,size(AllSqIntAll{i_gtRep},2));AllSqIntAll{i_gtRep}];
        end
        AllSqInt=saveCat(AllSqInt,AllSqIntAll{i_gtRep});
        srt=append_strc(srt,srtAll{i_gtRep});

    end
end

%generate unique animal IDs
%necessary because animals from different experiments have same numbers
%generate new number for each combination of animalNr&ExpSet
[~, ~, srt.AnimalNrCont]=unique([srt.AnimalNr;srt.ExpSet]','rows');
srt.AnimalNrCont=srt.AnimalNrCont';

%need to re-sort PlotGroup sorting vector to match data read-order.
[k, l, m]=unique(srt.ExpSet);
n=[ExpInfo.PlotGroup];
r=n(k);
srt.PlotGroup=r(m);
clear k l m r

%consolidate imaging traces
%generate chopped and stitched and normalized traces
%generate matching sorting vectors

%Plot individual traces
%[traces, srt]=consolidateTraces(AllSqInt,srt,0); %0 is for starting with long video

%generate color maps for trace plotting
if length(unq_PlotGroup)<11
    plotColsNorm=balanced();
else
    plotColsNorm=cbrewer('div', 'RdYlGn', length(unq_PlotGroup));
end
plotCols=cbrewer('div', 'PuOr', length(unq_PlotGroup));
% plotCols=cbrewer('div', 'RdYlGn', length(unq_PlotGroup));

%plots the color-coded individual traces, the averaged traces (both start
%only and clamped prestim)
%if plotTraces_yn
%    errTmp=plotTraces(traces,srt); 
% Figure(4)
%   errTmp=plotDoseResponseTrace_01(traces,srt); %go here to change highlighting protocol
% errTmp=plotStaircases_02(traces,srt); 


%function result=plotDoseResponseTrace_01(traces,srt)
YLimTraces=[-.30 2.0];
unq_PlotGroup=unique(srt.PlotGroup);
xLimTraces=[0 33000];
%xLimTraces=[0 15650];
%xLimTraces=[0 5000];

% 
% stim_b=10.^(1/5);
% stim_c=.01;
% stim(1)=0;
% stim(2)=stim_c;
% for i=3:13;stim(i)=stim(i-1).*stim_b;end
% stim(end+1:end+2)=0;
% stimTimes=600.*[0 12 13:24 27];

doseRespTim; %OUTPUT IS doseResponseTiming, ie go here to change the highlighting!!!!!

stimTimes=doseResponseTiming; 


%function hilite_Du(x,y,col,axlist,replace)
%       x is Nx2 matrix of [start end] x-pairs for rectangular highlight
%       y is Nx2 matrix of [start end] y-pairs for rectangular highlight
%           either x or y can be [] to span axis limits
%       color is 1x3 RGB color value (or [] = default)
%       axlist is list of axes ([] is current axis, or 'all' for all axes in current figure)
%       replace if true, clears prior highlighting (default = true)


%generate color maps for trace plotting

if length(unq_PlotGroup)<11
    plotColsNorm=balanced();
else
    plotColsNorm=cbrewer('div', 'RdYlGn', length(unq_PlotGroup));
end
% plotCols=cbrewer('div', 'PuOr', length(unq_PlotGroup));
plotCols=plotColsNorm;

%%
[x, y] = size(AllSqInt);
n = y;

%SPECIFY Group name and color here
strain = 'wt, ist1';
strainnum=8;
col1=1;
col2=1;
col3=1;

%CUT traces
start = 5950
finish = 25949
AllSqIntAllM1X = AllSqInt(start:finish, :);
len = finish - start +1;

%smooth original traces
sm = 50;
AllSqIntAllM1Xorg = [];

for an=1:y
    AllSqIntAllM1Xorg(:,an)=smooth(AllSqIntAllM1X(:,an),sm); %smooths traces
end

a = [1850 2150; 2450 2750; 3050 3350; 4250 4550; 4850 5150; 5450 5750; 6650 6950; 7250 7550; 7850 8150; 9050 9350; 9650 9950; 10250 10550; 11450 11750; 12050 12350; 12650 12950; 13850 14150; 14450 14750; 15050 15350; 16250 16550; 16850 17150; 17450 17750];

AllSqIntAllM1Xcorr=[];
for an=1:y
    regmin=[AllSqIntAllM1Xorg(16850:17750,an); AllSqIntAllM1Xorg(14450:15350,an)];
    regmin=sort(regmin);
    realmin=nanmean(regmin(10:15));
    regmax=[AllSqIntAllM1Xorg(13400:13850,an); AllSqIntAllM1Xorg(15800:16250,an)] ;
    regmax=sort(regmax);
    realmax=nanmean(regmax(885:890));
    for t=1:len
        if AllSqIntAllM1Xorg(t,an)<(realmin-(abs(realmin*0.3)))
            AllSqIntAllM1Xcorr(t,an)=NaN;
        elseif AllSqIntAllM1Xorg(t,an)>(realmax*1.3)
            AllSqIntAllM1Xcorr(t,an)=NaN;
        else
            AllSqIntAllM1Xcorr(t,an)=AllSqIntAllM1Xorg(t,an);
        end
    end 
end

%Normalize traces
AllSqIntAllM1XNormF = NaN(len,y); 
AllSqIntAllM1XNormFA = NaN(len,y); 
normFmax=NaN(1,y);
for an=1:y
    trace = AllSqIntAllM1Xcorr(:,an);
    %establish Fo
    Fo=trace(16850:17150);
    Fo=sort(Fo);
    Fo=Fo(11:60);
    Fo=nanmean(Fo);
    if Fo < 0
        trace = trace+abs(Fo)+1000;
        Fo=1000;
    end
    %dF/F
    normtraceF = ((trace-Fo)/Fo)*100;
    AllSqIntAllM1XNormF(:,an) = normtraceF; 
    %Normalize 0-1
    %normtraceF = trace;
    sorted_normtraceF=sort(normtraceF);
    sorted_normtraceF=sorted_normtraceF(~isnan(sorted_normtraceF));
    r=numel(sorted_normtraceF);
    p5=r*0.05;
    p5=round(p5);
    minF=mean(sorted_normtraceF(1:p5));
    sorted_normtraceF=sorted_normtraceF-minF;
    maxF=mean(sorted_normtraceF((r-p5+1):r));
    NormFmax(1,an)=maxF;
    normtraceFA = normtraceF - minF;    
    normtraceFA = normtraceFA/maxF;   
    AllSqIntAllM1XNormFA(:,an) = normtraceFA;     
end

%Calculate standard error

%Elias
averagetrace=nanmean(AllSqIntAllM1XNormFA,2);
% plot(averagetrace)
std_error = nanste(AllSqIntAllM1XNormFA,2);
data=AllSqIntAllM1XNormFA;
tracetype='norm01F';
% h = shadedErrorBar(1:length(averagetrace),averagetrace,std_error)
% h = shadedErrorBar(1:length(averagetrace),averagetrace,std_error,{'r-o','markerfacecolor','r'})


%% %Get values
% 
% data=AllSqIntAllM1XNormFA;
% tracetype='norm01F';
% xlsname='052016_MagAllGroups_norm01F.xlsx'; %user specified
% 
% %Measure magnitude of baseline and trough and rec
% Fbase = NaN(21,y); %mean of 3 sec before each pulse
% Fon = NaN(21,y); %mean of the last 3 seconds 
% for an = 1:y
%     for p = 1:21
%         fr = a(p,1);
%         fr2 = a(p,2);
%         Fbase(p,an) = nanmean(data((fr-21):(fr-1),an));
%         Fon(p,an) = nanmean(data((fr2-105):(fr2-5),an));              
%     end 
% end
% 
% %calculate ONmag
% ONmag = NaN(21,y);
% ONmag = Fbase-Fon;
% 
% %calculate recT
% smdata = [];
% for an=1:y
%     smdata(:,an)=smooth(data(:,an),10); %smooths traces
% end
% RpercentA = 0.33;
% RpercentB = 0.5;
% RpercentC = 0.75;
% recTA = NaN(21,y);
% recTB = NaN(21,y);
% recTC = NaN(21,y);
% deltaRec = NaN(21,y);
% for an = 1:y 
%     for c = 1:7
%         deltaRec(c*3-2, an) = Fbase(c*3-2,an)-Fon(c*3-2,an);
%         deltaRec(c*3-1, an) = Fbase(c*3-2,an)-Fon(c*3-1,an);
%         deltaRec(c*3, an) = Fbase(c*3-2,an)-Fon(c*3,an);
%     end
% end
%         
% 
% threshB = NaN(21,y);
% 
% for an = 1:y
%     for p = 1:21
%         fr = a(p,1)+305;
%         fr2 = a(p,2)+900; % if looking at pulse 1 or 2, this should be changed to +300
%         if p == 21
%             fr2 = 19999;
%         end
%         threshB(p,an) = (deltaRec(p,an)*RpercentB)+Fon(p,an);
%         for g = fr:fr2
%             if smdata(g,an) > threshB(p,an)
%                 recTB(p,an) = g - a(p,2);
%                 break               
%             end           
%         end   
%     end
% end   
%      
% 
% for an = 1:y
%     for p = 1:21
%         if deltaRec(p,an)<0.075
%             recTB(p,an) = NaN;
%         end
%     end
% end
% 
% 
% %select the ONmag and recT of pulses of interest
% ONmag1 = NaN(7,y);
% recTA3 = NaN(7,y);
% recTB3 = NaN(7,y);
% recTC3 = NaN(7,y);
% for an = 1:y
%     for c = 1:7
%         ONmag1(c,an) = ONmag(c*3-2, an);
%         recTA3(c,an) = recTA(c*3, an);
%         recTB3(c,an) = recTB(c*3, an);
%         recTC3(c,an) = recTC(c*3, an);
%     end
% end
% 
% %Calculate averages
% ONM = nanmean(ONmag1,2);
% ONste = ste(ONmag1,2);
% recTBM = nanmean(recTB3,2);
% recTBste = ste(recTB3,2);

% %   
% %Save data as xls
% sheettitle = {strcat(strain)};
% xlab = {'animal #'};
% %ylab = {'bufferp1';'bufferp2';'bufferp3'; '-9p1';'-9p2';'-9p3'; '-8p1';'-8p2';'-8p3'; '-7p1';'-7p2';'-7p3'; '-6p1';'-6p2';'-6p3'; '-5p1';'-5p2';'-5p3'; '-4p1';'-4p2';'-4p3'};
% ylab = {'bufferp';'-9p';'-8p';'-7p';'-6p';'-5p';'-4p';};
% xlswrite(xlsname,sheettitle,strainnum,'B2');
% xlswrite(xlsname,'ONmag1',strainnum,'D2');
% xlswrite(xlsname,xlab,strainnum,'C3');
% xlswrite(xlsname,ylab,strainnum,'B4');
% xlswrite(xlsname,ONmag1,strainnum,'C4');
% 
% xlswrite(xlsname,sheettitle,strainnum,'B14');
% xlswrite(xlsname,'recTB3',strainnum,'D14');
% xlswrite(xlsname,xlab,strainnum,'C15');
% xlswrite(xlsname,ylab,strainnum,'B16');
% xlswrite(xlsname,recTB3,strainnum,'C16');

        
%% Plot Traces

%AVERAGE traces
Fmean = nanmean(data,2);
Fste = nanste(data,2);

% %Chop traces by pulse & conc
% P123=NaN(2500,y,7);
%  for c=1:7 
%      inc = (c-1)*2400;
%      P123(:,1:y,c)=data(1750+inc:4249+inc,:);
%  end
% FMeanSteC=NaN(2500,2,7);
% for c=1:7  
%     FMeanSteC(:,1,c) = nanmean(P123(:,:,c),2);
%     FMeanSteC(:,2,c) = NaNste(P123(:,:,c),2);
% end
% 
% %ON and OFF magnitudes at each pulse*conc
% OnOff=NaN(9,y,7);
% for an=1:y 
%     for c=1:7
%         cmin1=sort(P123(101:400,an,c));
%         cmin1=nanmean(cmin1(10:19));
%         OnOff(1,an,c)=nanmean(P123(76:95,an,c))-cmin1;
%         cmax1=sort(P123(401:700,an,c));
%         cmax1=nanmean(cmax1(280:289));
%         OnOff(4,an,c)=cmax1-nanmean(P123(376:395,an,c));
%         cmin2=sort(P123(701:1000,an,c));
%         cmin2=nanmean(cmin2(10:19));
%         OnOff(2,an,c)=nanmean(P123(676:695,an,c))-cmin2;
%         cmax2=sort(P123(1001:1300,an,c));
%         cmax2=nanmean(cmax2(280:289));
%         OnOff(5,an,c)=cmax2-nanmean(P123(976:995,an,c));
%         cmin3=sort(P123(1301:1600,an,c));
%         cmin3=nanmean(cmin3(10:19));
%         OnOff(3,an,c)=nanmean(P123(1276:1295,an,c))-cmin3;
%         cmax3=sort(P123(1601:1900,an,c));
%         cmax3=nanmean(cmax3(280:289));
%         OnOff(6,an,c)=cmax3-nanmean(P123(1576:1595,an,c));
%     end
% end
% 
% %Average info by pulse & conc
% for c=1:7
%     ON1=OnOff(1,:,c);
%     OnOff(7,1,c)=nanmean(ON1,2);
%     OnOff(7,2,c)=NaNste(ON1,2);
%     ON23=[OnOff(3,:,c) OnOff(5,:,c)];
%     OnOff(8,1,c)=nanmean(ON23,2);
%     OnOff(8,2,c)=NaNste(ON23,2);
%     OFF123=[OnOff(2,:,c) OnOff(4,:,c) OnOff(6,:,c)];
%     OnOff(9,1,c)=nanmean(OFF123,2); 
%     OnOff(9,2,c)=NaNste(OFF123,2);
% end

% k=2;  

% %Plot average trace only
%figure(k);
figure(1);
% subplot(2,1,1)
% for an=1:y
%     plot((1:len),data(:,an),'Color', [col1,col2,col3],'LineWidth',1)
%     hold on
% end
% xlim([1 len]);
% ylim([-.1 1.1]);
% title(strcat(strain,'..individuals..', tracetype));
% hilite_Du([a]);
% subplot(2,1,2)
%h = shadedErrorBar(1:length(averagetrace),averagetrace,std_error)
%h = shadedErrorBar(1:length(averagetrace),averagetrace,std_error,{'r-o','markerfacecolor','r'})
h = shadedErrorBar(1:length(averagetrace),averagetrace,std_error,{'markerfacecolor','r'});
%  hold on
%  plot((1:len),Fmean,'Color', [col1,col2,col3],'LineWidth',2)
% jbfill((1:len),Fmean'+Fste',Fmean'-Fste',[col1,col2,col3]);
hold on
xlim([1 len]);
ylim([-.1 1.1]);
title(strcat(strain,'..Mean..', tracetype));
hilite_Du([a]);
hold on
k=k+1;

% % %Plot individual and average trace
% figure(k);
% subplot(2,1,1)
% for an=1:y
%     plot((1:len),data(:,an),'Color', [col1,col2,col3],'LineWidth',1)
%     hold on
% end
% xlim([1 len]);
% ylim([-.1 1.1]);
% title(strcat(strain,'..individuals..', tracetype));
% hilite_Du([a]);
% subplot(2,1,2)
% h = shadedErrorBar(1:length(averagetrace),averagetrace,std_error)
% hold on
% plot((1:len),Fmean,'Color', [col1,col2,col3],'LineWidth',2)
% %jbfill((1:len),Fmean'+Fste',Fmean'-Fste',[col1,col2,col3]);
% hold on
% xlim([1 len]);
% ylim([-.1 1.1]);
% title(strcat(strain,'..Mean..', tracetype));
% hilite_Du([a]);
% hold on
% k=k+1;

% %Plot individual traces separately
% figure(k);
% p=1;
% for an=1:y
%     subplot(2,2,p)    
%     plot((1:len),data(:,an),'Color', [col1,col2,col3],'LineWidth',1)
%     hold on
%     xlim([1 len]);
%     ylim([-.1 1.1]);
%     title(strcat(strain,'..individuals..', tracetype));
%     hilite_Du([a]);
%     hold off
%     p=p+1;
%     if p > 4 && an < y
%         k=k+1;
%         figure(k);
%         p=1;
%     end
%     if an==y
%    
%     end
%     
% end
% 
% k=k+1;

% plot heatmaps
% figure(k);
% imagesc(data');
% colormap(jet);
% xlabel('Time (seconds)','FontSize',14);
% title(strcat(strain), 'FontSize', 14);
% colorbar
% k=k+1;
% 

% figure(k);
% ConcData = data(1:1800,:);
% imagesc(ConcData');
% colormap(jet);
% xlabel('Time (seconds)','FontSize',14);
% title(strcat(strain,'..base'), 'FontSize', 14);
% colorbar
% k=k+1;
    
% 
% for y=1:7,
%     figure(k);
%     s1=1200+(y-1)*2400;
%     s2=s1+2700;
%     ConcData = data(s1:s2,:);
%     imagesc(ConcData');
%     colormap(jet);
%     xlabel('Time (seconds)','FontSize',14);
%     title(strcat(strain,'..[-', num2str(11-y),']'), 'FontSize', 14);
%     colorbar
%     k=k+1;
% end
% 
% for c=1:7
%     subplot(4,2,c+1)
%     plot((1:2500),FMeanSteC(:,1,c),'Color', [col1,col2,col3],'LineWidth',3)
%     hold on
%     jbfill((1:2500),FMeanSteC(:,1,c)'+FMeanSteC(:,2,c)',FMeanSteC(:,1,c)'-FMeanSteC(:,2,c)',[col1,col2,col3]);
%     hold on
%     xlim([1,1900])
%     ylim([ys(4,1) ys(4,2)]);
%     hold on
%     hilite_Du([101 400; 701 1000; 1301 1600])
%     hold on
%     title(num2str(-1*(11-c)));
%     hold on
% end
savefig(strcat(tracetype,'_',strain, '_trace'));
% k=k+1;

