%% AWC VGLUT-pH data import and basic analysis script
% Donovan Ventimilglia
% Matlab 2013b


%% Import Data and collect axon and background values

Datafiles = dir('*Data.txt'); 
numfiles = length(Datafiles);
mydata = cell(1, numfiles); % Cell containing imported txt file: data, textdata, and column headers
ROI_BG1=cell(1,numfiles); % integrated Background FL for left ROIs
ROI_BG2=cell(1,numfiles); % integrated Background FL for right ROIs
Total_BG_mean=[]; % Stores background measure for each video (combines all
ROI_axon_mean_FL=cell(1,numfiles);
Whole_axon_BGdiv=[];
Whole_axon_BGcorrected=[];
All_dF_ROIs=[];
num_of_trials=3;
plotlinkedtrials=0;
plotlinkedtrialsRAW=0;
otherBGcorrect=0;
ROIanalysis=0;

% Quick plot command
%  databrowse(DeltaF',[1:size(DeltaF,1)]./5);

Date=date;
prompt = 'Name this Analysis?';
AnalysisGroupName = ['_' input(prompt) '_' Date];

% Set Normalization point  Multiplulse: 10sP = 66:75; 60sP = 316:325; 3minP
% = 916:925; Single pulse: 60sP = 441:450;
prompt= 'Specify Normalization point:';
NP=input(prompt);

% import text files
for k = 1:numfiles
    mydata{1,k} = importdata(Datafiles(k).name);
    filenamesS{k}=Datafiles(k).name;
end

% Assign trial iteration number to each video
for k = 1:numfiles
    TrialPos=strfind(filenamesS{k},'T')+1;
    if isempty(TrialPos)==1
        TrialNum(k)=1;
    else
    TrialNum(k)=str2num(filenamesS{k}(TrialPos));
    clear TrialPos;
    end
end


% Get Whole axon FL measurements
for k = 1:size(mydata,2);
    Whole_axon_FL(:,k)=mydata{1,k}.data(:,3);
end

% Get axon ROI measurements for each video
for k = 1:size(mydata,2);
    counter=1;
    for n = 9:18:size(mydata{1,k}.data,2)
    ROI_axon_mean_FL{1,k}(:,counter)=mydata{1,k}.data(:,n);
    counter=counter+1;
    end
end

% Get Background values for each ROI (left and right of axon - Raw integrate FL)
for k = 1:size(mydata,2);
    counter=1;
    for n = 19:18:size(mydata{1,k}.data,2)
        ROI_BG1{1,k}(:,counter)=mydata{1,k}.data(:,n);
        counter=counter+1;
    end
    counter=1;
    for n = 25:18:size(mydata{1,k}.data,2)
        ROI_BG2{1,k}(:,counter)=mydata{1,k}.data(:,n);
        counter=counter+1;
    end
end

% Calculate overall mean background (sum FL of each BG Box, add left and
% right together and divdie by the total area)
for k=1:numfiles
Total_BG_mean(:,k)=(sum(ROI_BG1{1,k},2)+sum(ROI_BG2{1,k},2))./repmat(64*size(ROI_BG1{1,k},2)*2,size(ROI_BG1{1,k},1),1);
end



%% Background Bleach correction

% Create linear fits of background signal

for k=1:size(Total_BG_mean,2)

    BG_fit(k,:)=polyfit([1:size(Total_BG_mean,1)]',Total_BG_mean(:,k),1);
end

for k=1:size(Total_BG_mean,2)
    BG_fitval(:,k)=polyval([BG_fit(k,1),BG_fit(k,2)],[1:size(Total_BG_mean,1)]');
end

% If background signal decays (linear fit, slope < 0), assume bleaching and
% correct. Correction = BG - FL change of linear fit (linear fit with last
% value of fit subtracted off)

for k=1:size(Total_BG_mean,2)
    if BG_fit(k,1)<0
        BG_corrected(:,k)=Total_BG_mean(:,k)-(BG_fitval(:,k)-repmat(BG_fitval(size(BG_fitval,1),k),size(BG_fitval,1),1));
    else
        BG_corrected(:,k)=Total_BG_mean(:,k);
    end
end

%% signal bleach test

% Create linear fits of signal

for k=1:size(Whole_axon_FL,2)

    S_fit(k,:)=polyfit([1:size(Whole_axon_FL,1)]',Whole_axon_FL(:,k),1);
end

for k=1:size(Whole_axon_FL,2)
    S_fitval(:,k)=polyval([S_fit(k,1),S_fit(k,2)],[1:size(Whole_axon_FL,1)]');
end


%% Background correction and DeltaF

% Assumption: pHluorin signal does not always exhibit significant
% bleaching but the background often does (the pharyngeal tissue autoflourescence). 
% If the background exhibits bleaching and the signal exhibits bleaching, 
% correct the signal by background subtraction. If the background exhibits
% bleaching, but the signal does not, use background bleach corrected for
% background subtraction. If the background does not bleach, do nothing to
% each signal and perform background subtraction.


BG_subtraction_vector=ones(size(Whole_axon_FL,1),size(Whole_axon_FL,2));
Whole_axon_BGcorrected=ones(size(Whole_axon_FL,1),size(Whole_axon_FL,2));
BG_delta=ones(size(Whole_axon_FL,1),size(Whole_axon_FL,2));

for k=1:size(Whole_axon_FL,2)
    
    Bleachrate(k,1)=BG_fit(k,1)*625./mean(Total_BG_mean(1:10,k),1)*100; % background
    Bleachrate(k,2)=S_fit(k,1)*625./mean(Whole_axon_FL(1:10,k),1)*100; % axon
    
    if Bleachrate(k,1)<-0.5 && Bleachrate(k,2)<-0.5
        BG_subtraction_vector(:,k)=(Total_BG_mean(:,k)-repmat(min(Total_BG_mean(:,k)),size(Total_BG_mean,1),1));
        Whole_axon_BGcorrected(:,k)=Whole_axon_FL(:,k)-BG_subtraction_vector(:,k);
        BG_delta(:,k)=NaN;
    else
        BG_subtraction_vector(:,k)=NaN;
        BG_delta(:,k)=(BG_corrected(:,k)-repmat(min(BG_corrected(:,k)),size(BG_corrected(:,k),1),1));
        Whole_axon_BGcorrected(:,k)=Whole_axon_FL(:,k)-BG_delta(:,k);
    end
    
end


% DeltaF calcuation
Inital=repmat(mean(Whole_axon_BGcorrected(NP,:)),size(Whole_axon_BGcorrected,1),1);
DeltaF=(Whole_axon_BGcorrected-Inital)./Inital*100;
% databrowse(DeltaF',[1:size(DeltaF,1)]./5);%hilite([30 90]);


%% Create matrix of linked trials (RAW)
%---------------------------

if plotlinkedtrialsRAW==1

counter=1;
for k=1:3:numfiles
temp=Whole_axon_FL(:,k:k+2);
Linked_trials_raw(:,counter)=temp(:);
counter=counter+1;
clear temp
end

plot
figure;
for k=1:size(Linked_trials_raw,2);
subplot(size(Linked_trials_raw,2),1,k);plot(Linked_trials_raw(:,k));
end

end
%---------------------------

%% Create matrix of linked trials (DeltaF)
%---------------------------
if plotlinkedtrials==1

counter=1;
for k=1:3:numfiles
temp=DeltaF(:,k:k+2);
Linked_trials_sub(:,counter)=temp(:);
counter=counter+1;
clear temp
end

% plot
figure;
for k=1:size(Linked_trials_sub,2);
subplot(size(Linked_trials_sub,2),1,k);plot(Linked_trials_sub(:,k));
end

end
%---------------------------

%% other background correction methods

if otherBGcorrect==1

Background Division
for k = 1:numfiles
    Whole_axon_BGdiv(:,k)=Whole_axon_FL(:,k)./Total_BG_mean(:,k);
end

% Background subtraction
for k = 1:numfiles
    Whole_axon_BGsub(:,k)=Whole_axon_FL(:,k)-Total_BG_mean(:,k);
end

%---------------------------
% DeltaF for Background Division signal
Inital_div=repmat(mean(Whole_axon_BGdiv(NP,:)),size(Whole_axon_BGdiv,1),1);
DeltaF_div=(Whole_axon_BGdiv-Inital_div)./Inital_div*100;
databrowse(DeltaF_div',[1:size(DeltaF_div,1)]./5);

% DeltaF for Background Subtracted signal
Inital_sub=repmat(mean(Whole_axon_BGsub(NP,:)),size(Whole_axon_BGsub,1),1);
DeltaF_sub=(Whole_axon_BGsub-Inital_sub)./Inital_sub*100;
databrowse(DeltaF_sub',[1:size(DeltaF_sub,1)]./5);

end

%---------------------------



%% ROI anaylsis - all data

if ROIanalysis==1


% Background Correction at each ROI

for k = 1:numfiles
    for n =1:size(ROI_axon_mean_FL{1,k},2)        
ROI_BG_combined{1,k}(:,n)=(ROI_BG1{1,k}(:,n)+ROI_BG2{1,k}(:,n))./128;            
ROI_BG_minsub{1,k}(:,n)=ROI_BG_combined{1,k}(:,n)-repmat(min(ROI_BG_combined{1,k}(:,n)),size(ROI_BG_combined{1,k},1),1);
ROI_axon_BGsub{1,k}(:,n)=ROI_axon_mean_FL{1,k}(:,n)-ROI_BG_minsub{1,k}(:,n);

    end
end


%% ROI DeltaF calcuation
%---------------------------


% DeltaF_sub at each ROI for each movie
for k = 1:numfiles
    for n =1:size(ROI_axon_mean_FL{1,k},2)
    ROI_Initail_sub{1,k}(:,n)=repmat(mean(ROI_axon_BGsub{1,k}(NP,n)),size(ROI_axon_BGsub{1,k},1),1);
    ROI_dF{1,k}(:,n)=(ROI_axon_BGsub{1,k}(:,n)-ROI_Initail_sub{1,k}(:,n))./ROI_Initail_sub{1,k}(:,n)*100;
    end
end

% collection of all DeltaF_sub ROIs
for k = 1:numfiles
All_dF_ROIs=[All_dF_ROIs ROI_dF{1,k}];
end
databrowse(All_dF_ROIs',[1:size(All_dF_ROIs,1)]./5);

%---------------------------




%% Mean ROI response per axon - average of repeated trials
% 
% for k=1:2 %num of animals
%     ROIs=cat(3,ROI_dF{1,Animal_num==k});
%     mean_along_axon=mean(ROIs,3);
% end
    
   
end



FileName=['AWCpHluorinImport_Analysis' AnalysisGroupName];
save(FileName);




