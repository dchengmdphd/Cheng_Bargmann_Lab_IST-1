%This is function is developted by Du Cheng, Alejandro Lopez, Navin Pokala,
%and Elias Sheer. This is intended to plot distance to odor vs time, with
%individual traces in the background and average value plus error shade in
%the forground. odor_linked_tracks from Navin's Chemotaxis_tracker is used
%as input, the array of distance to odor is converted into an attribut
%matrix. In the matrix the last non-Nan valuse of each track is used to
%fill the rest of the track till the end of the recording time (lenght of the longest track).

%not ready for prime time NP and DC 5/31/16
%edited ES, DC 8/10/16 count last non-nan value in matrix
%edited AL, DC 8/12/16 plot from matrix instead of array, changed line
%edited DC 01/04/16 accumulation
%edited DC 10/2/2017 optimize accumulation counting and calculating CI
%coloring
function overlay_chemotaxis_distance_vs_time(target_odor_linkedTracks1, color1, target_odor_linkedTracks2, color2, filename, attribute)

if(nargin<6)
    attribute = 'odor_distance';
end

min_tracklength_minutes = 3;
max_delta_dist = 2;
min_distance_to_odor_mm = 15;
max_distance_to_odor_mm = 75;
transperancy = 0.25;
distance_odor_to_edge_mm = 10;
FrameRate = 3;

target_odor_linkedTracks1 = Tracks_minimum_length(target_odor_linkedTracks1, target_odor_linkedTracks1(1).FrameRate*min_tracklength_minutes);
target_odor_linkedTracks2 = Tracks_minimum_length(target_odor_linkedTracks2, target_odor_linkedTracks2(1).FrameRate*min_tracklength_minutes);

input_array1 = target_odor_linkedTracks1;
input_array2 = target_odor_linkedTracks2;
input_array3 = [];
input_array4 = [];

% remove jitter frames ... inter-frame distance > max_delta_dist. Some
% tracks will have a sudden increase of value to very high distance then
% fall back. This is to remove that.
for m = 1:length(target_odor_linkedTracks1)
    dd = abs(diff(target_odor_linkedTracks1(m).(attribute)));
    idx = find(dd>max_delta_dist);
    target_odor_linkedTracks1(m).(attribute)(idx) = NaN;
end
for m = 1:length(target_odor_linkedTracks2)
    dd = abs(diff(target_odor_linkedTracks2(m).(attribute)));
    idx = find(dd>max_delta_dist);
    target_odor_linkedTracks2(m).(attribute)(idx) = NaN;
end

attribute_matrix1 = final_attribute_matrix(target_odor_linkedTracks1, min_distance_to_odor_mm, max_distance_to_odor_mm, attribute);
attribute_matrix2 = final_attribute_matrix(target_odor_linkedTracks2, min_distance_to_odor_mm, max_distance_to_odor_mm, attribute);

mean_dist_1 = nanmean(attribute_matrix1);
mean_dist_1_err = nanstderr(attribute_matrix1);

mean_dist_2 = nanmean(attribute_matrix2);
mean_dist_2_err = nanstderr(attribute_matrix2);

time1 = (1:size(attribute_matrix1,2))/target_odor_linkedTracks1(1).FrameRate;
time2 = (1:size(attribute_matrix2,2))/target_odor_linkedTracks2(1).FrameRate;

%% mean plots: plotting the mean of each condition only without background individual traces
figure(1);
NameArray_fig = {'PaperUnits','PaperSize', 'PaperPosition', 'Units','Position','PaperOrientation'};
ValueArray_fig = {'inches', [8.5000 11], [0 1 8.5000 9],'inches', [0 1 8.5000 9],'portrait'};
set(figure(1),NameArray_fig,ValueArray_fig)

subplot(3,2,1)
title(['{\color{blue}naive} vs {\color{red}trained}']);
errorshade(time1, mean_dist_1+mean_dist_1_err, mean_dist_1-mean_dist_1_err,color1,transperancy);
hold on;
errorshade(time2, mean_dist_2+mean_dist_2_err, mean_dist_2-mean_dist_2_err,color2,transperancy);

Area_under_curve1 = trapz(time1, mean_dist_1) - 3600*distance_odor_to_edge_mm;
Area_under_curve2 = trapz(time2, mean_dist_2) - 3600*distance_odor_to_edge_mm;

label1 = strjoin({'Area under naive curve =',num2str(Area_under_curve1)});
label2 = strjoin({'Area under adapt curve =',num2str(Area_under_curve2)});
plot(time1, mean_dist_1, color1, 'linewidth',3);
plot(time2, mean_dist_2, color2, 'linewidth',3);
text(200,86, label1,'FontSize', 9)
text(200,80, label2,'FontSize', 9)

xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)');
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec)');
orient portrait;

%% % overlay raw plots: plotting the mean plot on top of the indivitual
% % traces. figure 2 is condition 1, figure 3 is condition 2
% figure(2);
% title(fix_title_string(target_odor_linkedTracks1(1).Name));
% for m = 1:length(target_odor_linkedTracks1)
%     %plotting from attribute matrix
%     plot(time1,attribute_matrix1(m,:), 'Color',[0.85,0.85,0.85],'linewidth',0.5); 
%     %plotting from origional array
%     %plot(((target_odor_linkedTracks1(m).Frames)/(target_odor_linkedTracks1(m).FrameRate)),target_odor_linkedTracks1(m).(attribute), 'Color',[0.85,0.85,0.85],'linewidth',0.5); hold on;
%     
%     hold on;
% end
% hold on;
% plot(time1, mean_dist_1, color1, 'linewidth',3);
% xlim([0 max([time1 time2])]);
% if(strcmp(attribute,'odor_distance'))
%     ylim([0 90]);
%     ylabel('Distance to odor (mm)');
% else
%     ylabel(fix_title_string(attribute));
% end
% xlabel('Time (sec)');
% orient landscape;
% 
% figure(3);
% title(fix_title_string(target_odor_linkedTracks2(1).Name));
% for m = 1:length(target_odor_linkedTracks2)
%     %plotting from attribute matrix
%     plot(time2,attribute_matrix2(m,:), 'Color',[0.85,0.85,0.85],'linewidth',0.5);     
%     %plotting from origional array
%     %plot(((target_odor_linkedTracks2(m).Frames)/(target_odor_linkedTracks2(m).FrameRate)),target_odor_linkedTracks2(m).(attribute), 'Color',[0.85,0.85,0.85],'linewidth',0.5); hold on;
%     
%     hold on;
% end
% hold on;
% plot(time2, mean_dist_2, color2, 'linewidth',3);
% xlim([0 max([time1 time2])]);
% if(strcmp(attribute,'odor_distance'))
%     ylim([0 90]);
%     ylabel('Distance to odor (mm)');
% else
%     ylabel(fix_title_string(attribute));
% end
% xlabel('Time (sec)');
% orient landscape;


%% overlay both conditions
subplot(3,2,2)
%title(fix_title_string(target_odor_linkedTracks1(1).Name));
title(['{\color{blue}naive} vs {\color{red}trained}']);


for m = 1:length(target_odor_linkedTracks2)
    
    plot(time2,attribute_matrix2(m,:), 'Color',[1,0.8,0.8],'linewidth',0.5);
    %plot(((target_odor_linkedTracks2(m).Frames)/(target_odor_linkedTracks2(m).FrameRate)),target_odor_linkedTracks2(m).(attribute), 'Color',[0.8,0.8,1],'linewidth',0.5); hold on;
    hold on;
end
hold on;
for m = 1:length(target_odor_linkedTracks1)
    
    plot(time1,attribute_matrix1(m,:), 'Color',[0.8,0.8,1],'linewidth',0.5);
    %plot(((target_odor_linkedTracks1(m).Frames)/(target_odor_linkedTracks1(m).FrameRate)),target_odor_linkedTracks1(m).(attribute), 'Color',[0.8,0.8,1],'linewidth',0.5); hold on;
    hold on;
end
hold on;
errorshade(time1, mean_dist_1+mean_dist_1_err, mean_dist_1-mean_dist_1_err,color1,transperancy);
plot(time1, mean_dist_1, color1, 'linewidth',3); 
hold on;
errorshade(time2, mean_dist_2+mean_dist_2_err, mean_dist_2-mean_dist_2_err,color2,transperancy);
plot(time2, mean_dist_2, color2, 'linewidth',3);

plot(time1, mean_dist_1, color1, 'linewidth',3);
plot(time2, mean_dist_2, color2, 'linewidth',3);

xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)');
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec)');
orient portrait;

%% overlay raw plots with error shades

subplot(3,2,3)
%title(fix_title_string(target_odor_linkedTracks1(1).Name));
title(['{\color{blue}naive}']);
for m = 1:length(target_odor_linkedTracks1)
    plot(time1,attribute_matrix1(m,:), 'Color',[0.85,0.85,0.85],'linewidth',0.5); 
    %plot(((target_odor_linkedTracks1(m).Frames)/(target_odor_linkedTracks1(m).FrameRate)),target_odor_linkedTracks1(m).(attribute), 'Color',[0.85,0.85,0.85],'linewidth',0.5); hold on;

    hold on;
end
hold on;
errorshade(time1, mean_dist_1+mean_dist_1_err, mean_dist_1-mean_dist_1_err,color1,transperancy);
plot(time1, mean_dist_1, color1, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)');
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec)');
orient landscape;

subplot(3,2,4)
%title(fix_title_string(target_odor_linkedTracks2(1).Name));
title(['{\color{red}trained}']);
for m = 1:length(target_odor_linkedTracks2)
    
    plot(time2,attribute_matrix2(m,:), 'Color',[0.85,0.85,0.85],'linewidth',0.5); 
    %plot(((target_odor_linkedTracks2(m).Frames)/(target_odor_linkedTracks2(m).FrameRate)),target_odor_linkedTracks2(m).(attribute), 'Color',[0.85,0.85,0.85],'linewidth',0.5); hold on;
    
    hold on;
end
hold on;
errorshade(time2, mean_dist_2+mean_dist_2_err, mean_dist_2-mean_dist_2_err,color2,transperancy);
plot(time2, mean_dist_2, color2, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)');
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec)');
orient portrait;


%% Speed Average

time_window = 1; %setting time-window to 1min. Can change here.

[events_twindows_1,time_1] = events_per_animal_twindows(input_array1,time_window,'all',FrameRate);
avg_1 = nanmean(events_twindows_1);
stdev_1 = nanstd(events_twindows_1);
sem_1 = stdev_1./sqrt((sum(~isnan(events_twindows_1))));


%for arrays 2-4 will only create calculate these things if the array is not
%empty. If it is empty I will make these variables empty
if isempty(input_array2)==0
    [events_twindows_2,time_2] = events_per_animal_twindows(input_array2,time_window,'all',FrameRate);
    avg_2 = nanmean(events_twindows_2);
    stdev_2 = nanstd(events_twindows_2);
    sem_2 = stdev_2./sqrt((sum(~isnan(events_twindows_2))));
else
    time_2=NaN(1,1);  avg_2 = NaN(1,1);  sem_2 = NaN(1,1);
end

if isempty(input_array3)==0
    [events_twindows_3,time_3] = events_per_animal_twindows(input_array3,time_window,'all',FrameRate);
    avg_3 = nanmean(events_twindows_3);
    stdev_3 = nanstd(events_twindows_3);
    sem_3 = stdev_3./sqrt((sum(~isnan(events_twindows_3))));
else
    time_3=NaN(1,1);  avg_3 = NaN(1,1); sem_3 = NaN(1,1);
end

if isempty(input_array4)==0
    [events_twindows_4,time_4] = events_per_animal_twindows(input_array4,time_window,'all',FrameRate);
    avg_4 = nanmean(events_twindows_4);
    stdev_4 = nanstd(events_twindows_4);
    sem_4 = stdev_4./sqrt((sum(~isnan(events_twindows_4))));

else
    time_4=NaN(1,1);  avg_4 = NaN(1,1); sem_4 = NaN(1,1);
end

[avg_1,sem_1, time_1] = avg_twindow_fwd(input_array1,'Speed',time_window,FrameRate);

maxx = [nanmax(time_1), nanmax(time_2), nanmax(time_3), nanmax(time_4)];
maxx= nanmax(maxx);

%for arrays 2-4 will only create calculate these things if the array is not
%empty. If it is empty I will make these variables empty
if isempty(input_array2)==0
    [avg_2,sem_2, time_2] = avg_twindow_fwd(input_array2,'Speed',time_window,FrameRate);
else
    time_2=NaN(1,1);  avg_2 = NaN(1,1);  sem_2 = NaN(1,1);
end

if isempty(input_array3)==0
    [avg_3,sem_3, time_3] = avg_twindow_fwd(input_array3,'Speed',time_window,FrameRate);
else
    time_3=NaN(1,1);  avg_3 = NaN(1,1); sem_3 = NaN(1,1);
end

if isempty(input_array4)==0
    [avg_4,sem_4, time_4] = avg_twindow_fwd(input_array4,'Speed',time_window,FrameRate);
else
    time_4=NaN(1,1);  avg_4 = NaN(1,1); sem_4 = NaN(1,1);
end

%Plotting the average reorientation rates with the SEM
subplot(3,2,5)
one = errorline(time_1,avg_1,sem_1); hold on; two =  errorline(time_2,avg_2,sem_2);  hold on; three =   errorline(time_3,avg_3,sem_3);  hold on;  four =  errorline(time_4,avg_4,sem_4);%properties for plot of input_array1
%properties for plot of input_array1
NameArray_1 = {'LineWidth','Color'};
ValueArray_1 = {0.75,[0 0 1]};
set(one,NameArray_1,ValueArray_1)
set(one(end),'LineWidth',1.5) %setting the last handle, which is the line to be thicker in width
%properties for plot of input_array
NameArray_2 = {'LineWidth','Color'}; 
ValueArray_2 = {0.75,[1 0 0]};
set(two,NameArray_2,ValueArray_2)
set(two(end),'LineWidth',1.5)
%properties for plot of input_array3
NameArray_3 = {'LineWidth','Color'};
ValueArray_3 = {0.75,[0 0.25 1]};
set(three,NameArray_3,ValueArray_3)
set(three(end),'LineWidth',1.5)
%properties for plot of input_array4
NameArray_4 = {'LineWidth','Color'};
ValueArray_4 = {0.75,[0.25 0.25 0.25]};
set(four,NameArray_4,ValueArray_4)
set(four(end),'LineWidth',1.5)
%gca properties
NameArray_gca = {'FontSize','XLim','YLim','FontName'};
ValueArray_gca = {12, [0 maxx],[0 0.2],'Arial'};
set(gca,NameArray_gca,ValueArray_gca)
Xl = xlabel('Time (min)');
Yl = ylabel('Speed(mm/s)');
set(Xl,'FontSize',14,'FontName','Arial')
set(Yl,'FontSize',14,'FontName','Arial')
% labels = {label_1 label_2 label_3 label_4};
% labels = labels(1:length(varargin));
% handle = [one(end) two(end) three(end) four(end)]; %so that it only incldues plot and nopt the fill in the legend
% handle = handle(1:length(varargin));
% l = legend(handle,labels);
% legend('boxoff');
% set(l,'FontSize',9)
% set(legend,'Location','southeast');
title('Speed','FontSize',12)

%% Mean Square Displacement (MSD)

%MSD- will calculate and plot the mean MSD as a
%function of time since food removal.

smooth_window = 0.2; %total smooth window of 1 min. Can change here
time_window = 2; %setting time-window to 2mins. Can change here.

%Calculates the mean MSD for each animal in a specific time window (see MSD_per_animal_twindows).
%Then calculates the average and SEM for each time window for the
%population
%Does this for each input.  If input is empty it just creates an empty matrix.

%First will plot nonUpsilon event- only thing that I am changing from
%preious code is the reor_type argument ('allnonUpsilon') inputed to events_per_animal_twindows
[MSDavg_twindows_1,time_1] = MSD_per_animal_twindows(input_array1,smooth_window,time_window,FrameRate);
avg_1 = nanmean(MSDavg_twindows_1);
stdev_1 = nanstd(MSDavg_twindows_1);
sem_1 = stdev_1./sqrt((sum(~isnan(MSDavg_twindows_1))));


%for arrays 2-4 will only create calculate these things if the array is not
%empty. If it is empty I will make these variables empty
if isempty(input_array2)==0
    [MSDavg_twindows_2,time_2] = MSD_per_animal_twindows(input_array2,smooth_window,time_window,FrameRate);
    avg_2 = nanmean(MSDavg_twindows_2);
    stdev_2 = nanstd(MSDavg_twindows_2);
    sem_2 = stdev_2./sqrt((sum(~isnan(MSDavg_twindows_2))));
else
    time_2=NaN(1,1);  avg_2 = NaN(1,1);  sem_2 = NaN(1,1);
end

if isempty(input_array3)==0
    [MSDavg_twindows_3,time_3] = MSD_per_animal_twindows(input_array3,smooth_window,time_window,FrameRate);
    avg_3 = nanmean(MSDavg_twindows_3);
    stdev_3 = nanstd(MSDavg_twindows_3);
    sem_3 = stdev_3./sqrt((sum(~isnan(MSDavg_twindows_3))));
else
    time_3=NaN(1,1);  avg_3 = NaN(1,1); sem_3 = NaN(1,1);
end

if isempty(input_array4)==0
    [MSDavg_twindows_4,time_4] = MSD_per_animal_twindows(input_array4,smooth_window,time_window,FrameRate);
    avg_4 = nanmean(MSDavg_twindows_4);
    stdev_4 = nanstd(MSDavg_twindows_4);
    sem_4 = stdev_4./sqrt((sum(~isnan(MSDavg_twindows_4))));
else
    time_4=NaN(1,1);  avg_4 = NaN(1,1); sem_4 = NaN(1,1);
end

%Plotting the average MSD with the SEM
maxx = [nanmax(time_1), nanmax(time_2), nanmax(time_3), nanmax(time_4)];
minx  = nanmin(maxx); %used later on in cumulative reorientation events
maxx= nanmax(maxx); %not being used for now (commented out on 2/19/15)

subplot(3,2,6)
one = errorline(time_1,avg_1,sem_1); hold on; two =  errorline(time_2,avg_2,sem_2);  hold on; three =   errorline(time_3,avg_3,sem_3);  hold on;  four =  errorline(time_4,avg_4,sem_4);%properties for plot of input_array1
%properties for plot of input_array1
NameArray_1 = {'LineWidth','Color'};
ValueArray_1 = {1.5,[0 0 1]};
set(one,NameArray_1,ValueArray_1)
%properties for plot of input_array2
NameArray_2 = {'LineWidth','Color'};
ValueArray_2 = {1.5,[1 0 0]};
set(two,NameArray_2,ValueArray_2)
%properties for plot of input_array3
NameArray_3 = {'LineWidth','Color'};
ValueArray_3 = {1.5,[0 0.25 1]};
set(three,NameArray_3,ValueArray_3)
%properties for plot of input_array4
NameArray_4 = {'LineWidth','Color'};
ValueArray_4 = {1.5,[0.25 0.25 0.25]};
set(four,NameArray_4,ValueArray_4)
%gca properties
NameArray_gca = {'FontSize','XLim','YLim','FontName'};
ValueArray_gca = {12, [0 maxx], [0 inf], 'Arial'};
set(gca,NameArray_gca,ValueArray_gca)
Xl = xlabel('Time (min)');
Yl = ylabel({'Mean Squared';'Displacement(mm^2)'});
set(Xl,'FontSize',14,'FontName','Arial')
set(Yl,'FontSize',14,'FontName','Arial')
% labels = {label_1 label_2 label_3 label_4};
% labels = labels(1:3);
% handle = [one(end) two(end) three(end) four(end)];
% handle = handle(1:length(varargin));
% l = legend(handle,labels);
% legend('boxoff');
% set(l,'FontSize',9)
% set(l,'Location','southeast')
title('MSD','FontSize',12)

%% figure 2
figure(2);
NameArray_fig = {'PaperUnits','PaperSize', 'PaperPosition', 'Units','Position','PaperOrientation'};
ValueArray_fig = {'inches', [8.5000 11], [0 1 8.5000 9],'inches', [0 1 8.5000 9],'portrait'};
set(figure(2),NameArray_fig,ValueArray_fig)
%% all reorientation naive
subplot(3,2,1)
%title(fix_title_string(target_odor_linkedTracks1(1).Name));
title('All Reoreintations');
for m = 1:length(target_odor_linkedTracks1)
    %plotting from attribute matrix
    plot(time1,attribute_matrix1(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    hold on;
    %plotting from origional array
    %plot(((target_odor_linkedTracks1(m).Frames)/(target_odor_linkedTracks1(m).FrameRate)),target_odor_linkedTracks1(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks1)
    stateNo = 1;
    stmatrix1 = track_field_to_matrix(target_odor_linkedTracks1,'State');
    st1 = find((abs(stmatrix1(m,:) - stateNo))>0.2);
    plot(st1/3, attribute_matrix1(m,st1), 'x','MarkerSize',2, 'MarkerEdgeColor','b')
    hold on;
    end
    hold on;
plot(time1, mean_dist_1, color1, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:All Reorientations', 'FontSize', 12);
orient portrait;
%%  %all reorientation adapt
subplot(3,2,2)
%title(fix_title_string(target_odor_linkedTracks2(1).Name));
title('All Reoreintations');
for m = 1:length(target_odor_linkedTracks2)
    %plotting from attribute matrix
    plot(time2,attribute_matrix2(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    %plotting from origional array
    %plot(((target_odor_linkedTracks2(m).Frames)/(target_odor_linkedTracks2(m).FrameRate)),target_odor_linkedTracks2(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks2)
    stateNo = 1;
    stmatrix2 = track_field_to_matrix(target_odor_linkedTracks2,'State');
    st2 = find((abs(stmatrix2(m,:) - stateNo))>0.2);
    plot(st2/3, attribute_matrix2(m,st2), 'x','MarkerSize',2,'MarkerEdgeColor','r')
    hold on;
    end
    hold on;
plot(time2, mean_dist_2, color2, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:All Reorientations', 'FontSize', 12);
orient portrait;
%% pause naive
subplot(3,2,3)
%title(fix_title_string(target_odor_linkedTracks1(1).Name));
title('Pause');
for m = 1:length(target_odor_linkedTracks1)
    %plotting from attribute matrix
    plot(time1,attribute_matrix1(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    hold on;
    %plotting from origional array
    %plot(((target_odor_linkedTracks1(m).Frames)/(target_odor_linkedTracks1(m).FrameRate)),target_odor_linkedTracks1(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks1)
    stateNo = 1.1;
    stmatrix1 = track_field_to_matrix(target_odor_linkedTracks1,'State');
    st1 = find((abs(stmatrix1(m,:) - stateNo))>0.2);
    plot(st1/3, attribute_matrix1(m,st1), 'x','MarkerSize',2, 'MarkerEdgeColor','b')
    hold on;
    end
    hold on;
plot(time1, mean_dist_1, color1, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:Pause', 'FontSize', 12);
orient portrait;
%% pause adapt
subplot(3,2,4)
%title(fix_title_string(target_odor_linkedTracks2(1).Name));
title('Pause');
for m = 1:length(target_odor_linkedTracks2)
    %plotting from attribute matrix
    plot(time2,attribute_matrix2(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    %plotting from origional array
    %plot(((target_odor_linkedTracks2(m).Frames)/(target_odor_linkedTracks2(m).FrameRate)),target_odor_linkedTracks2(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks2)
    stateNo = 1.1;
    stmatrix2 = track_field_to_matrix(target_odor_linkedTracks2,'State');
    st2 = find((abs(stmatrix2(m,:) - stateNo))>0.2);
    plot(st2/3, attribute_matrix2(m,st2), 'x','MarkerSize',2,'MarkerEdgeColor','r')
    hold on;
    end
    hold on;
plot(time2, mean_dist_2, color2, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:Pause', 'FontSize', 12);
orient portrait;
%% All Reorientation Frequencies


%Reorientation frequency- will calculate and plot the reorientation frequency as a
%function of time since food removal. This includes all types of
%reorientations
time_window = 2; %setting time-window to 2mins. Can change here.

%Calculates the events that each animal did in a specific time wondow (see events_per_animal_twindows).
%Calculate the average and SEM for each time window
%Does this for each input.  If input is empty it just creates an empty matrix.


[events_twindows_1,time_1] = events_per_animal_twindows(input_array1,time_window,'all',FrameRate);
avg_1 = nanmean(events_twindows_1);
stdev_1 = nanstd(events_twindows_1);
sem_1 = stdev_1./sqrt((sum(~isnan(events_twindows_1))));


%for arrays 2-4 will only create calculate these things if the array is not
%empty. If it is empty I will make these variables empty
if isempty(input_array2)==0
    [events_twindows_2,time_2] = events_per_animal_twindows(input_array2,time_window,'all',FrameRate);
    avg_2 = nanmean(events_twindows_2);
    stdev_2 = nanstd(events_twindows_2);
    sem_2 = stdev_2./sqrt((sum(~isnan(events_twindows_2))));
else
    time_2=NaN(1,1);  avg_2 = NaN(1,1);  sem_2 = NaN(1,1);
end

if isempty(input_array3)==0
    [events_twindows_3,time_3] = events_per_animal_twindows(input_array3,time_window,'all',FrameRate);
    avg_3 = nanmean(events_twindows_3);
    stdev_3 = nanstd(events_twindows_3);
    sem_3 = stdev_3./sqrt((sum(~isnan(events_twindows_3))));
else
    time_3=NaN(1,1);  avg_3 = NaN(1,1); sem_3 = NaN(1,1);
end

if isempty(input_array4)==0
    [events_twindows_4,time_4] = events_per_animal_twindows(input_array4,time_window,'all',FrameRate);
    avg_4 = nanmean(events_twindows_4);
    stdev_4 = nanstd(events_twindows_4);
    sem_4 = stdev_4./sqrt((sum(~isnan(events_twindows_4))));

else
    time_4=NaN(1,1);  avg_4 = NaN(1,1); sem_4 = NaN(1,1);
end


maxx = [nanmax(time_1), nanmax(time_2), nanmax(time_3), nanmax(time_4)];
maxx= nanmax(maxx);
%maxx = 50; %can type a number here if you want to only go to a certain x(time) value)

sub = subplot(3,2,5);
one = errorline(time_1,avg_1,sem_1); hold on; two =  errorline(time_2,avg_2,sem_2);  hold on; three =   errorline(time_3,avg_3,sem_3);  hold on;  four =  errorline(time_4,avg_4,sem_4);
NameArray_1 = {'LineWidth','Color'};
ValueArray_1 = {0.75,[0 0 1]};
set(one,NameArray_1,ValueArray_1)
set(one(end),'LineWidth',1.5) %setting the last handle, which is the line to be thicker in width
%properties for plot of input_array
NameArray_2 = {'LineWidth','Color'}; 
ValueArray_2 = {0.75,[1 0 0]};
set(two,NameArray_2,ValueArray_2)
set(two(end),'LineWidth',1.5)
%properties for plot of input_array3
NameArray_3 = {'LineWidth','Color'};
ValueArray_3 = {0.75,[0 0.25 1]};
set(three,NameArray_3,ValueArray_3)
set(three(end),'LineWidth',1.5)
%properties for plot of input_array4
NameArray_4 = {'LineWidth','Color'};
ValueArray_4 = {0.75,[0.25 0.25 0.25]};
set(four,NameArray_4,ValueArray_4)
set(four(end),'LineWidth',1.5)
%gca properties
NameArray_gca = {'FontSize','XLim','FontName'};
ValueArray_gca = {12, [0 maxx],'Arial'};
set(sub,NameArray_gca,ValueArray_gca)
Xl = xlabel('Time (min)');
Yl = ylabel('Frenquency');
set(Xl,'FontSize',14,'FontName','Arial')
set(Yl,'FontSize',14,'FontName','Arial')
% labels = {label_1 label_2 label_3 label_4};
% labels = labels(1:length(varargin));
% handle = [one(end) two(end) three(end) four(end)]; %so that it only incldues plot and nopt the fill in the legend
% handle = handle(1:length(varargin));
% l = legend(handle,labels);
% legend('boxoff');
% set(l,'FontSize',9)
title('All Reorientations Frequency','FontSize',12)

% %subplots properties
% ylim = get(sub,'YLim'); %getting the y lim so that I can use it for other reorientation plots.
% ymax = max(ylim);

%% Pause frequency

%Reorientation frequency- will calculate and plot the reorientation frequency as a
%function of time since food removal. This includes lRevOmega types of
%reorientations
time_window = 2; %setting time-window to 2mins. Can change here.
input_array1 = []; input_array2 = []; input_array3 = []; input_array4 = [];
input_array1 = target_odor_linkedTracks1;
input_array2 = target_odor_linkedTracks2;
%Calculates the events that each animal did in a specific time wondow (see events_per_animal_twindows).
%Calculate the average and SEM for each time window
%Does this for each input.  If input is empty it just creates an empty matrix.

[events_twindows_1,time_1] = events_per_animal_twindows(input_array1,time_window,'pause',FrameRate);
avg_1 = nanmean(events_twindows_1);
stdev_1 = nanstd(events_twindows_1);
sem_1 = stdev_1./sqrt((sum(~isnan(events_twindows_1))));


%for arrays 2-4 will only create calculate these things if the array is not
%empty. If it is empty I will make these variables empty
if isempty(input_array2)==0
    [events_twindows_2,time_2] = events_per_animal_twindows(input_array2,time_window,'pause',FrameRate);
    avg_2 = nanmean(events_twindows_2);
    stdev_2 = nanstd(events_twindows_2);
    sem_2 = stdev_2./sqrt((sum(~isnan(events_twindows_2))));
else
    time_2=NaN(1,1);  avg_2 = NaN(1,1);  sem_2 = NaN(1,1);
end

if isempty(input_array3)==0
    [events_twindows_3,time_3] = events_per_animal_twindows(input_array3,time_window,'pause',FrameRate);
    avg_3 = nanmean(events_twindows_3);
    stdev_3 = nanstd(events_twindows_3);
    sem_3 = stdev_3./sqrt((sum(~isnan(events_twindows_3))));
else
    time_3=NaN(1,1);  avg_3 = NaN(1,1); sem_3 = NaN(1,1);
end

if isempty(input_array4)==0
    [events_twindows_4,time_4] = events_per_animal_twindows(input_array4,time_window,'pause',FrameRate);
    avg_4 = nanmean(events_twindows_4);
    stdev_4 = nanstd(events_twindows_4);
    sem_4 = stdev_4./sqrt((sum(~isnan(events_twindows_4))));

else
    time_4=NaN(1,1);  avg_4 = NaN(1,1); sem_4 = NaN(1,1);
end


maxx = [nanmax(time_1), nanmax(time_2), nanmax(time_3), nanmax(time_4)];
maxx= nanmax(maxx);
%maxx = 50; %can type a number here if you want to only go to a certain x(time) value)

sub = subplot(3,2,6);
one = errorline(time_1,avg_1,sem_1); hold on; two =  errorline(time_2,avg_2,sem_2);  hold on; three =   errorline(time_3,avg_3,sem_3);  hold on;  four =  errorline(time_4,avg_4,sem_4);
NameArray_1 = {'LineWidth','Color'};
ValueArray_1 = {0.75,[0 0 1]};
set(one,NameArray_1,ValueArray_1)
set(one(end),'LineWidth',1.5) %setting the last handle, which is the line to be thicker in width
%properties for plot of input_array
NameArray_2 = {'LineWidth','Color'}; 
ValueArray_2 = {0.75,[1 0 0]};
set(two,NameArray_2,ValueArray_2)
set(two(end),'LineWidth',1.5)
%properties for plot of input_array3
NameArray_3 = {'LineWidth','Color'};
ValueArray_3 = {0.75,[0 0.25 1]};
set(three,NameArray_3,ValueArray_3)
set(three(end),'LineWidth',1.5)
%properties for plot of input_array4
NameArray_4 = {'LineWidth','Color'};
ValueArray_4 = {0.75,[0.25 0.25 0.25]};
set(four,NameArray_4,ValueArray_4)
set(four(end),'LineWidth',1.5)
%gca properties
NameArray_gca = {'FontSize','XLim','FontName'};
ValueArray_gca = {12, [0 maxx],'Arial'};
set(sub,NameArray_gca,ValueArray_gca)
Xl = xlabel('Time (min)');
Yl = ylabel('Frequency');
set(Xl,'FontSize',14,'FontName','Arial')
set(Yl,'FontSize',14,'FontName','Arial')
% labels = {label_1 label_2 label_3 label_4};
% labels = labels(1:length(varargin));
% handle = [one(end) two(end) three(end) four(end)]; %so that it only incldues plot and nopt the fill in the legend
% handle = handle(1:length(varargin));
% l = legend(handle,labels);
% legend('boxoff');
% set(l,'FontSize',9)
title('Pause Frequency','FontSize',12)

% %subplots properties
% ylim = get(sub,'YLim'); %getting the y lim so that I can use it for other reorientation plots.
% ymax = max(ylim);

%% figure 3
figure(3);
NameArray_fig = {'PaperUnits','PaperSize', 'PaperPosition', 'Units','Position','PaperOrientation'};
ValueArray_fig = {'inches', [8.5000 11], [0 1 8.5000 9],'inches', [0 1 8.5000 9],'portrait'};
set(figure(3),NameArray_fig,ValueArray_fig)
%% %LRevOmega naive
subplot(3,2,1)
%title(fix_title_string(target_odor_linkedTracks1(1).Name));
for m = 1:length(target_odor_linkedTracks1)
    %plotting from attribute matrix
    plot(time1,attribute_matrix1(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    hold on;
    %plotting from origional array
    %plot(((target_odor_linkedTracks1(m).Frames)/(target_odor_linkedTracks1(m).FrameRate)),target_odor_linkedTracks1(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks1)
    stateNo = 4.7;
    stmatrix1 = track_field_to_matrix(target_odor_linkedTracks1,'State');
    st1 = find((abs(stmatrix1(m,:) - stateNo))<0.01);
    plot(st1/3, attribute_matrix1(m,st1), 'x','MarkerSize',2, 'MarkerEdgeColor','b')
    hold on;
    end
    hold on;
plot(time1, mean_dist_1, color1, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:LRevOmega', 'FontSize', 12);
orient portrait;

%% LRevOmega adapt
subplot(3,2,2);
%title(fix_title_string(target_odor_linkedTracks2(1).Name));
for m = 1:length(target_odor_linkedTracks2)
    %plotting from attribute matrix
    plot(time2,attribute_matrix2(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    %plotting from origional array
    %plot(((target_odor_linkedTracks2(m).Frames)/(target_odor_linkedTracks2(m).FrameRate)),target_odor_linkedTracks2(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks2)
    stateNo = 4.7;
    stmatrix2 = track_field_to_matrix(target_odor_linkedTracks2,'State');
    st2 = find((abs(stmatrix2(m,:) - stateNo))<0.01);
    plot(st2/3, attribute_matrix2(m,st2), 'x','MarkerSize',2,'MarkerEdgeColor','r')
    hold on;
    end
    hold on;
plot(time2, mean_dist_2, color2, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:LRevOmega', 'FontSize', 12);
orient portrait;
%% %LRevUpsilon naive
subplot(3,2,3)
%title(fix_title_string(target_odor_linkedTracks1(1).Name));
for m = 1:length(target_odor_linkedTracks1)
    %plotting from attribute matrix
    plot(time1,attribute_matrix1(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    hold on;
    %plotting from origional array
    %plot(((target_odor_linkedTracks1(m).Frames)/(target_odor_linkedTracks1(m).FrameRate)),target_odor_linkedTracks1(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks1)
    stateNo = 4.3;
    stmatrix1 = track_field_to_matrix(target_odor_linkedTracks1,'State');
    st1 = find((abs(stmatrix1(m,:) - stateNo))<0.01);
    plot(st1/3, attribute_matrix1(m,st1), 'x','MarkerSize',2, 'MarkerEdgeColor','b')
    hold on;
    end
    hold on;
plot(time1, mean_dist_1, color1, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:LRevUpsilon', 'FontSize', 12);
orient portrait;

%% LRevUpsilon adapt
subplot(3,2,4);
%title(fix_title_string(target_odor_linkedTracks2(1).Name));
for m = 1:length(target_odor_linkedTracks2)
    %plotting from attribute matrix
    plot(time2,attribute_matrix2(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    %plotting from origional array
    %plot(((target_odor_linkedTracks2(m).Frames)/(target_odor_linkedTracks2(m).FrameRate)),target_odor_linkedTracks2(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks2)
    stateNo = 4.3;
    stmatrix2 = track_field_to_matrix(target_odor_linkedTracks2,'State');
    st2 = find((abs(stmatrix2(m,:) - stateNo))<0.01);
    plot(st2/3, attribute_matrix2(m,st2), 'x','MarkerSize',2,'MarkerEdgeColor','r')
    hold on;
    end
    hold on;
plot(time2, mean_dist_2, color2, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:LRevUpsilon', 'FontSize', 12);
orient portrait;
%% Long Reversal Omega


%Reorientation frequency- will calculate and plot the reorientation frequency as a
%function of time since food removal. This includes lRevOmega types of
%reorientations
time_window = 2; %setting time-window to 2mins. Can change here.
input_array1 = []; input_array2 = []; input_array3 = []; input_array4 = [];
input_array1 = target_odor_linkedTracks1;
input_array2 = target_odor_linkedTracks2;
%Calculates the events that each animal did in a specific time wondow (see events_per_animal_twindows).
%Calculate the average and SEM for each time window
%Does this for each input.  If input is empty it just creates an empty matrix.

[events_twindows_1,time_1] = events_per_animal_twindows(input_array1,time_window,'lRevOmega',FrameRate);
avg_1 = nanmean(events_twindows_1);
stdev_1 = nanstd(events_twindows_1);
sem_1 = stdev_1./sqrt((sum(~isnan(events_twindows_1))));


%for arrays 2-4 will only create calculate these things if the array is not
%empty. If it is empty I will make these variables empty
if isempty(input_array2)==0
    [events_twindows_2,time_2] = events_per_animal_twindows(input_array2,time_window,'lRevOmega',FrameRate);
    avg_2 = nanmean(events_twindows_2);
    stdev_2 = nanstd(events_twindows_2);
    sem_2 = stdev_2./sqrt((sum(~isnan(events_twindows_2))));
else
    time_2=NaN(1,1);  avg_2 = NaN(1,1);  sem_2 = NaN(1,1);
end

if isempty(input_array3)==0
    [events_twindows_3,time_3] = events_per_animal_twindows(input_array3,time_window,'lRevOmega',FrameRate);
    avg_3 = nanmean(events_twindows_3);
    stdev_3 = nanstd(events_twindows_3);
    sem_3 = stdev_3./sqrt((sum(~isnan(events_twindows_3))));
else
    time_3=NaN(1,1);  avg_3 = NaN(1,1); sem_3 = NaN(1,1);
end

if isempty(input_array4)==0
    [events_twindows_4,time_4] = events_per_animal_twindows(input_array4,time_window,'lRevOmega',FrameRate);
    avg_4 = nanmean(events_twindows_4);
    stdev_4 = nanstd(events_twindows_4);
    sem_4 = stdev_4./sqrt((sum(~isnan(events_twindows_4))));

else
    time_4=NaN(1,1);  avg_4 = NaN(1,1); sem_4 = NaN(1,1);
end


maxx = [nanmax(time_1), nanmax(time_2), nanmax(time_3), nanmax(time_4)];
maxx= nanmax(maxx);
%maxx = 50; %can type a number here if you want to only go to a certain x(time) value)

sub = subplot(3,2,5);
one = errorline(time_1,avg_1,sem_1); hold on; two =  errorline(time_2,avg_2,sem_2);  hold on; three =   errorline(time_3,avg_3,sem_3);  hold on;  four =  errorline(time_4,avg_4,sem_4);
NameArray_1 = {'LineWidth','Color'};
ValueArray_1 = {0.75,[0 0 1]};
set(one,NameArray_1,ValueArray_1)
set(one(end),'LineWidth',1.5) %setting the last handle, which is the line to be thicker in width
%properties for plot of input_array
NameArray_2 = {'LineWidth','Color'}; 
ValueArray_2 = {0.75,[1 0 0]};
set(two,NameArray_2,ValueArray_2)
set(two(end),'LineWidth',1.5)
%properties for plot of input_array3
NameArray_3 = {'LineWidth','Color'};
ValueArray_3 = {0.75,[0 0.25 1]};
set(three,NameArray_3,ValueArray_3)
set(three(end),'LineWidth',1.5)
%properties for plot of input_array4
NameArray_4 = {'LineWidth','Color'};
ValueArray_4 = {0.75,[0.25 0.25 0.25]};
set(four,NameArray_4,ValueArray_4)
set(four(end),'LineWidth',1.5)
%gca properties
NameArray_gca = {'FontSize','XLim','FontName'};
ValueArray_gca = {12, [0 maxx],'Arial'};
set(sub,NameArray_gca,ValueArray_gca)
Xl = xlabel('Time (min)');
Yl = ylabel('Reorientations');
set(Xl,'FontSize',14,'FontName','Arial')
set(Yl,'FontSize',14,'FontName','Arial')
% labels = {label_1 label_2 label_3 label_4};
% labels = labels(1:length(varargin));
% handle = [one(end) two(end) three(end) four(end)]; %so that it only incldues plot and nopt the fill in the legend
% handle = handle(1:length(varargin));
% l = legend(handle,labels);
% legend('boxoff');
% set(l,'FontSize',9)
title('lRevOmega Reorientations Frequency','FontSize',12)

% %subplots properties
% ylim = get(sub,'YLim'); %getting the y lim so that I can use it for other reorientation plots.
% ymax = max(ylim);
%% Long Reversal Upsilon


%Reorientation frequency- will calculate and plot the reorientation frequency as a
%function of time since food removal. This includes lRevOmega types of
%reorientations
time_window = 2; %setting time-window to 2mins. Can change here.
input_array1 = []; input_array2 = []; input_array3 = []; input_array4 = [];
input_array1 = target_odor_linkedTracks1;
input_array2 = target_odor_linkedTracks2;
%Calculates the events that each animal did in a specific time wondow (see events_per_animal_twindows).
%Calculate the average and SEM for each time window
%Does this for each input.  If input is empty it just creates an empty matrix.

[events_twindows_1,time_1] = events_per_animal_twindows(input_array1,time_window,'lRevUpsilon',FrameRate);
avg_1 = nanmean(events_twindows_1);
stdev_1 = nanstd(events_twindows_1);
sem_1 = stdev_1./sqrt((sum(~isnan(events_twindows_1))));


%for arrays 2-4 will only create calculate these things if the array is not
%empty. If it is empty I will make these variables empty
if isempty(input_array2)==0
    [events_twindows_2,time_2] = events_per_animal_twindows(input_array2,time_window,'lRevUpsilon',FrameRate);
    avg_2 = nanmean(events_twindows_2);
    stdev_2 = nanstd(events_twindows_2);
    sem_2 = stdev_2./sqrt((sum(~isnan(events_twindows_2))));
else
    time_2=NaN(1,1);  avg_2 = NaN(1,1);  sem_2 = NaN(1,1);
end

if isempty(input_array3)==0
    [events_twindows_3,time_3] = events_per_animal_twindows(input_array3,time_window,'lRevUpsilon',FrameRate);
    avg_3 = nanmean(events_twindows_3);
    stdev_3 = nanstd(events_twindows_3);
    sem_3 = stdev_3./sqrt((sum(~isnan(events_twindows_3))));
else
    time_3=NaN(1,1);  avg_3 = NaN(1,1); sem_3 = NaN(1,1);
end

if isempty(input_array4)==0
    [events_twindows_4,time_4] = events_per_animal_twindows(input_array4,time_window,'lRevUpsilon',FrameRate);
    avg_4 = nanmean(events_twindows_4);
    stdev_4 = nanstd(events_twindows_4);
    sem_4 = stdev_4./sqrt((sum(~isnan(events_twindows_4))));

else
    time_4=NaN(1,1);  avg_4 = NaN(1,1); sem_4 = NaN(1,1);
end


maxx = [nanmax(time_1), nanmax(time_2), nanmax(time_3), nanmax(time_4)];
maxx= nanmax(maxx);
%maxx = 50; %can type a number here if you want to only go to a certain x(time) value)

sub = subplot(3,2,6);
one = errorline(time_1,avg_1,sem_1); hold on; two =  errorline(time_2,avg_2,sem_2);  hold on; three =   errorline(time_3,avg_3,sem_3);  hold on;  four =  errorline(time_4,avg_4,sem_4);
NameArray_1 = {'LineWidth','Color'};
ValueArray_1 = {0.75,[0 0 1]};
set(one,NameArray_1,ValueArray_1)
set(one(end),'LineWidth',1.5) %setting the last handle, which is the line to be thicker in width
%properties for plot of input_array
NameArray_2 = {'LineWidth','Color'}; 
ValueArray_2 = {0.75,[1 0 0]};
set(two,NameArray_2,ValueArray_2)
set(two(end),'LineWidth',1.5)
%properties for plot of input_array3
NameArray_3 = {'LineWidth','Color'};
ValueArray_3 = {0.75,[0 0.25 1]};
set(three,NameArray_3,ValueArray_3)
set(three(end),'LineWidth',1.5)
%properties for plot of input_array4
NameArray_4 = {'LineWidth','Color'};
ValueArray_4 = {0.75,[0.25 0.25 0.25]};
set(four,NameArray_4,ValueArray_4)
set(four(end),'LineWidth',1.5)
%gca properties
NameArray_gca = {'FontSize','XLim','FontName'};
ValueArray_gca = {12, [0 maxx],'Arial'};
set(sub,NameArray_gca,ValueArray_gca)
Xl = xlabel('Time (min)');
Yl = ylabel('Reorientations');
set(Xl,'FontSize',14,'FontName','Arial')
set(Yl,'FontSize',14,'FontName','Arial')
% labels = {label_1 label_2 label_3 label_4};
% labels = labels(1:length(varargin));
% handle = [one(end) two(end) three(end) four(end)]; %so that it only incldues plot and nopt the fill in the legend
% handle = handle(1:length(varargin));
% l = legend(handle,labels);
% legend('boxoff');
% set(l,'FontSize',9)
title('lRevUpsilon Reorientations Frequency','FontSize',12)

% %subplots properties
% ylim = get(sub,'YLim'); %getting the y lim so that I can use it for other reorientation plots.
% ymax = max(ylim);




%% figure 4
figure(4);
NameArray_fig = {'PaperUnits','PaperSize', 'PaperPosition', 'Units','Position','PaperOrientation'};
ValueArray_fig = {'inches', [8.5000 11], [0 1 8.5000 9],'inches', [0 1 8.5000 9],'portrait'};
set(figure(4),NameArray_fig,ValueArray_fig)

%% sRevOmega naive
subplot(3,2,1)
%title(fix_title_string(target_odor_linkedTracks1(1).Name));
for m = 1:length(target_odor_linkedTracks1)
    %plotting from attribute matrix
    plot(time1,attribute_matrix1(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    hold on;
    %plotting from origional array
    %plot(((target_odor_linkedTracks1(m).Frames)/(target_odor_linkedTracks1(m).FrameRate)),target_odor_linkedTracks1(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks1)
    stateNo = 5.7;
    stmatrix1 = track_field_to_matrix(target_odor_linkedTracks1,'State');
    st1 = find((abs(stmatrix1(m,:) - stateNo))<0.01);
    plot(st1/3, attribute_matrix1(m,st1), 'x','MarkerSize',2, 'MarkerEdgeColor','b')
    hold on;
    end
    hold on;
plot(time1, mean_dist_1, color1, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:sRev Omega', 'FontSize', 12);
orient portrait;

%% sRevOmega adapt
subplot(3,2,2)
%title(fix_title_string(target_odor_linkedTracks2(1).Name));
for m = 1:length(target_odor_linkedTracks2)
    %plotting from attribute matrix
    plot(time2,attribute_matrix2(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    %plotting from origional array
    %plot(((target_odor_linkedTracks2(m).Frames)/(target_odor_linkedTracks2(m).FrameRate)),target_odor_linkedTracks2(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks2)
    stateNo = 5.7;
    stmatrix2 = track_field_to_matrix(target_odor_linkedTracks2,'State');
    st2 = find((abs(stmatrix2(m,:) - stateNo))<0.01);
    plot(st2/3, attribute_matrix2(m,st2), 'x','MarkerSize',2,'MarkerEdgeColor','r')
    hold on;
    end
    hold on;
plot(time2, mean_dist_2, color2, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:sRev Omega', 'FontSize', 12);
orient portrait;

%% sRevUpsilon naive
subplot(3,2,3)
%title(fix_title_string(target_odor_linkedTracks1(1).Name));
for m = 1:length(target_odor_linkedTracks1)
    %plotting from attribute matrix
    plot(time1,attribute_matrix1(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    hold on;
    %plotting from origional array
    %plot(((target_odor_linkedTracks1(m).Frames)/(target_odor_linkedTracks1(m).FrameRate)),target_odor_linkedTracks1(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks1)
    stateNo = 5.3;
    stmatrix1 = track_field_to_matrix(target_odor_linkedTracks1,'State');
    st1 = find((abs(stmatrix1(m,:) - stateNo))<0.01);
    plot(st1/3, attribute_matrix1(m,st1), 'x','MarkerSize',2, 'MarkerEdgeColor','b')
    hold on;
    end
    hold on;
plot(time1, mean_dist_1, color1, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:sRev Omega', 'FontSize', 12);
orient portrait;

%% sRevUpsilon adapt
subplot(3,2,4)
%title(fix_title_string(target_odor_linkedTracks2(1).Name));
for m = 1:length(target_odor_linkedTracks2)
    %plotting from attribute matrix
    plot(time2,attribute_matrix2(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    %plotting from origional array
    %plot(((target_odor_linkedTracks2(m).Frames)/(target_odor_linkedTracks2(m).FrameRate)),target_odor_linkedTracks2(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks2)
    stateNo = 5.3;
    stmatrix2 = track_field_to_matrix(target_odor_linkedTracks2,'State');
    st2 = find((abs(stmatrix2(m,:) - stateNo))<0.01);
    plot(st2/3, attribute_matrix2(m,st2), 'x','MarkerSize',2,'MarkerEdgeColor','r')
    hold on;
    end
    hold on;
plot(time2, mean_dist_2, color2, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:sRev Omega', 'FontSize', 12);
orient portrait;

%% sRev Omega Frequencies

%Same as above. Only thing that is different is that in
%events_per_animal_twindows the input for reor_type is 'pure_omega')

time_window = 2; %setting time-window to 2mins. Can change here.


[events_twindows_1,time_1] = events_per_animal_twindows(input_array1,time_window,'sRevOmega',FrameRate);
avg_1 = nanmean(events_twindows_1);
stdev_1 = nanstd(events_twindows_1);
sem_1 = stdev_1./sqrt((sum(~isnan(events_twindows_1))));


%for arrays 2-4 will only create calculate these things if the array is not
%empty. If it is empty I will make these variables empty
if isempty(input_array2)==0
    [events_twindows_2,time_2] = events_per_animal_twindows(input_array2,time_window,'sRevOmega',FrameRate);
    avg_2 = nanmean(events_twindows_2);
    stdev_2 = nanstd(events_twindows_2);
    sem_2 = stdev_2./sqrt((sum(~isnan(events_twindows_2))));
else
    time_2=NaN(1,1);  avg_2 = NaN(1,1);  sem_2 = NaN(1,1);
end

if isempty(input_array3)==0
    [events_twindows_3,time_3] = events_per_animal_twindows(input_array3,time_window,'sRevOmega',FrameRate);
    avg_3 = nanmean(events_twindows_3);
    stdev_3 = nanstd(events_twindows_3);
    sem_3 = stdev_3./sqrt((sum(~isnan(events_twindows_3))));
else
    time_3=NaN(1,1);  avg_3 = NaN(1,1); sem_3 = NaN(1,1);
end

if isempty(input_array4)==0
    [events_twindows_4,time_4] = events_per_animal_twindows(input_array4,time_window,'sRevOmega',FrameRate);
    avg_4 = nanmean(events_twindows_4);
    stdev_4 = nanstd(events_twindows_4);
    sem_4 = stdev_4./sqrt((sum(~isnan(events_twindows_4))));

else
    time_4=NaN(1,1);  avg_4 = NaN(1,1); sem_4 = NaN(1,1);
end

%Plotting the average reorientation rates with the SEM
subplot(3,2,5)
one = errorline(time_1,avg_1,sem_1); hold on; two =  errorline(time_2,avg_2,sem_2);  hold on; three =   errorline(time_3,avg_3,sem_3);  hold on;  four =  errorline(time_4,avg_4,sem_4);
%properties for plot of input_array1
NameArray_1 = {'LineWidth','Color'};
ValueArray_1 = {0.75,[0 0 1]};
set(one,NameArray_1,ValueArray_1)
set(one(end),'LineWidth',1.5) %setting the last handle, which is the line to be thicker in width
%properties for plot of input_array
NameArray_2 = {'LineWidth','Color'}; 
ValueArray_2 = {0.75,[1 0 0]};
set(two,NameArray_2,ValueArray_2)
set(two(end),'LineWidth',1.5)
%properties for plot of input_array3
NameArray_3 = {'LineWidth','Color'};
ValueArray_3 = {0.75,[0 0.25 1]};
set(three,NameArray_3,ValueArray_3)
set(three(end),'LineWidth',1.5)
%properties for plot of input_array4
NameArray_4 = {'LineWidth','Color'};
ValueArray_4 = {0.75,[0.25 0.25 0.25]};
set(four,NameArray_4,ValueArray_4)
set(four(end),'LineWidth',1.5)
%gca properties
NameArray_gca = {'FontSize','XLim','YLim','FontName'};
ValueArray_gca = {12, [0 maxx],[0 1],'Arial'};
set(gca,NameArray_gca,ValueArray_gca)
Xl = xlabel('Time (min)');
Yl = ylabel('Reorientations');
set(Xl,'FontSize',14,'FontName','Arial')
set(Yl,'FontSize',14,'FontName','Arial')
% labels = {label_1 label_2 label_3 label_4};
% labels = labels(1:length(varargin));
% handle = [one(end) two(end) three(end) four(end)]; %so that it only incldues plot and nopt the fill in the legend
% handle = handle(1:length(varargin));
% l = legend(handle,labels);
% legend('boxoff');
% set(l,'FontSize',9)
title('sRev Omega Frequency','FontSize',12)

%% sRevUpsilon Frequencies

%Same as above. Only thing that is different is that in
%events_per_animal_twindows the input for reor_type is 'pure_omega')

time_window = 2; %setting time-window to 2mins. Can change here.


[events_twindows_1,time_1] = events_per_animal_twindows(input_array1,time_window,'sRevUpsilon',FrameRate);
avg_1 = nanmean(events_twindows_1);
stdev_1 = nanstd(events_twindows_1);
sem_1 = stdev_1./sqrt((sum(~isnan(events_twindows_1))));


%for arrays 2-4 will only create calculate these things if the array is not
%empty. If it is empty I will make these variables empty
if isempty(input_array2)==0
    [events_twindows_2,time_2] = events_per_animal_twindows(input_array2,time_window,'sRevUpsilon',FrameRate);
    avg_2 = nanmean(events_twindows_2);
    stdev_2 = nanstd(events_twindows_2);
    sem_2 = stdev_2./sqrt((sum(~isnan(events_twindows_2))));
else
    time_2=NaN(1,1);  avg_2 = NaN(1,1);  sem_2 = NaN(1,1);
end

if isempty(input_array3)==0
    [events_twindows_3,time_3] = events_per_animal_twindows(input_array3,time_window,'sRevUpsilon',FrameRate);
    avg_3 = nanmean(events_twindows_3);
    stdev_3 = nanstd(events_twindows_3);
    sem_3 = stdev_3./sqrt((sum(~isnan(events_twindows_3))));
else
    time_3=NaN(1,1);  avg_3 = NaN(1,1); sem_3 = NaN(1,1);
end

if isempty(input_array4)==0
    [events_twindows_4,time_4] = events_per_animal_twindows(input_array4,time_window,'sRevUpsilon',FrameRate);
    avg_4 = nanmean(events_twindows_4);
    stdev_4 = nanstd(events_twindows_4);
    sem_4 = stdev_4./sqrt((sum(~isnan(events_twindows_4))));

else
    time_4=NaN(1,1);  avg_4 = NaN(1,1); sem_4 = NaN(1,1);
end

%Plotting the average reorientation rates with the SEM
subplot(3,2,6)
one = errorline(time_1,avg_1,sem_1); hold on; two =  errorline(time_2,avg_2,sem_2);  hold on; three =   errorline(time_3,avg_3,sem_3);  hold on;  four =  errorline(time_4,avg_4,sem_4);
%properties for plot of input_array1
NameArray_1 = {'LineWidth','Color'};
ValueArray_1 = {0.75,[0 0 1]};
set(one,NameArray_1,ValueArray_1)
set(one(end),'LineWidth',1.5) %setting the last handle, which is the line to be thicker in width
%properties for plot of input_array
NameArray_2 = {'LineWidth','Color'}; 
ValueArray_2 = {0.75,[1 0 0]};
set(two,NameArray_2,ValueArray_2)
set(two(end),'LineWidth',1.5)
%properties for plot of input_array3
NameArray_3 = {'LineWidth','Color'};
ValueArray_3 = {0.75,[0 0.25 1]};
set(three,NameArray_3,ValueArray_3)
set(three(end),'LineWidth',1.5)
%properties for plot of input_array4
NameArray_4 = {'LineWidth','Color'};
ValueArray_4 = {0.75,[0.25 0.25 0.25]};
set(four,NameArray_4,ValueArray_4)
set(four(end),'LineWidth',1.5)
%gca properties
NameArray_gca = {'FontSize','XLim','YLim','FontName'};
ValueArray_gca = {12, [0 maxx],[0 1],'Arial'};
set(gca,NameArray_gca,ValueArray_gca)
Xl = xlabel('Time (min)');
Yl = ylabel('Reorientations');
set(Xl,'FontSize',14,'FontName','Arial')
set(Yl,'FontSize',14,'FontName','Arial')
% labels = {label_1 label_2 label_3 label_4};
% labels = labels(1:length(varargin));
% handle = [one(end) two(end) three(end) four(end)]; %so that it only incldues plot and nopt the fill in the legend
% handle = handle(1:length(varargin));
% l = legend(handle,labels);
% legend('boxoff');
% set(l,'FontSize',9)
title('sRev Upsilon Frequency','FontSize',12)

%% figure 5
figure(5);
NameArray_fig = {'PaperUnits','PaperSize', 'PaperPosition', 'Units','Position','PaperOrientation'};
ValueArray_fig = {'inches', [8.5000 11], [0 1 8.5000 9],'inches', [0 1 8.5000 9],'portrait'};
set(figure(5),NameArray_fig,ValueArray_fig)

%% pure_Omega naive
subplot(3,2,1)
%title(fix_title_string(target_odor_linkedTracks1(1).Name));
for m = 1:length(target_odor_linkedTracks1)
    %plotting from attribute matrix
    plot(time1,attribute_matrix1(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    hold on;
    %plotting from origional array
    %plot(((target_odor_linkedTracks1(m).Frames)/(target_odor_linkedTracks1(m).FrameRate)),target_odor_linkedTracks1(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks1)
    stateNo = 7.1;
    stmatrix1 = track_field_to_matrix(target_odor_linkedTracks1,'State');
    st1 = find((abs(stmatrix1(m,:) - stateNo))<0.01);
    plot(st1/3, attribute_matrix1(m,st1), 'x','MarkerSize',2, 'MarkerEdgeColor','b')
    hold on;
    end
    hold on;
plot(time1, mean_dist_1, color1, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:pure Omega', 'FontSize', 12);
orient portrait;

%% pure_Omega adapt
subplot(3,2,2)
%title(fix_title_string(target_odor_linkedTracks2(1).Name));
for m = 1:length(target_odor_linkedTracks2)
    %plotting from attribute matrix
    plot(time2,attribute_matrix2(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    %plotting from origional array
    %plot(((target_odor_linkedTracks2(m).Frames)/(target_odor_linkedTracks2(m).FrameRate)),target_odor_linkedTracks2(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks2)
    stateNo = 7.1;
    stmatrix2 = track_field_to_matrix(target_odor_linkedTracks2,'State');
    st2 = find((abs(stmatrix2(m,:) - stateNo))<0.01);
    plot(st2/3, attribute_matrix2(m,st2), 'x','MarkerSize',2,'MarkerEdgeColor','r')
    hold on;
    end
    hold on;
plot(time2, mean_dist_2, color2, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:pure Omega', 'FontSize', 12);
orient portrait;

%% pure_upsilon naive
subplot(3,2,3)
%title(fix_title_string(target_odor_linkedTracks1(1).Name));
for m = 1:length(target_odor_linkedTracks1)
    %plotting from attribute matrix
    plot(time1,attribute_matrix1(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    hold on;
    %plotting from origional array
    %plot(((target_odor_linkedTracks1(m).Frames)/(target_odor_linkedTracks1(m).FrameRate)),target_odor_linkedTracks1(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks1)
    stateNo = 3.1;
    stmatrix1 = track_field_to_matrix(target_odor_linkedTracks1,'State');
    st1 = find((abs(stmatrix1(m,:) - stateNo))<0.01);
    plot(st1/3, attribute_matrix1(m,st1), 'x','MarkerSize',2, 'MarkerEdgeColor','b')
    hold on;
    end
    hold on;
plot(time1, mean_dist_1, color1, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:pure Upsilon(time independent)', 'FontSize', 12);
orient portrait;

%% pure_upsilon adapt
subplot(3,2,4)
%title(fix_title_string(target_odor_linkedTracks2(1).Name));
for m = 1:length(target_odor_linkedTracks2)
    %plotting from attribute matrix
    plot(time2,attribute_matrix2(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    %plotting from origional array
    %plot(((target_odor_linkedTracks2(m).Frames)/(target_odor_linkedTracks2(m).FrameRate)),target_odor_linkedTracks2(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks2)
    stateNo = 3.1;
    stmatrix2 = track_field_to_matrix(target_odor_linkedTracks2,'State');
    st2 = find((abs(stmatrix2(m,:) - stateNo))<0.01);
    plot(st2/3, attribute_matrix2(m,st2), 'x','MarkerSize',2,'MarkerEdgeColor','r')
    hold on;
    end
    hold on;
plot(time2, mean_dist_2, color2, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:pure Upsilon(time independent)', 'FontSize', 12);
orient portrait;
%% Pure Omega Frequencies

%Same as above. Only thing that is different is that in
%events_per_animal_twindows the input for reor_type is 'pure_omega')

time_window = 2; %setting time-window to 2mins. Can change here.


[events_twindows_1,time_1] = events_per_animal_twindows(input_array1,time_window,'pure_omega',FrameRate);
avg_1 = nanmean(events_twindows_1);
stdev_1 = nanstd(events_twindows_1);
sem_1 = stdev_1./sqrt((sum(~isnan(events_twindows_1))));


%for arrays 2-4 will only create calculate these things if the array is not
%empty. If it is empty I will make these variables empty
if isempty(input_array2)==0
    [events_twindows_2,time_2] = events_per_animal_twindows(input_array2,time_window,'pure_omega',FrameRate);
    avg_2 = nanmean(events_twindows_2);
    stdev_2 = nanstd(events_twindows_2);
    sem_2 = stdev_2./sqrt((sum(~isnan(events_twindows_2))));
else
    time_2=NaN(1,1);  avg_2 = NaN(1,1);  sem_2 = NaN(1,1);
end

if isempty(input_array3)==0
    [events_twindows_3,time_3] = events_per_animal_twindows(input_array3,time_window,'pure_omega',FrameRate);
    avg_3 = nanmean(events_twindows_3);
    stdev_3 = nanstd(events_twindows_3);
    sem_3 = stdev_3./sqrt((sum(~isnan(events_twindows_3))));
else
    time_3=NaN(1,1);  avg_3 = NaN(1,1); sem_3 = NaN(1,1);
end

if isempty(input_array4)==0
    [events_twindows_4,time_4] = events_per_animal_twindows(input_array4,time_window,'pure_omega',FrameRate);
    avg_4 = nanmean(events_twindows_4);
    stdev_4 = nanstd(events_twindows_4);
    sem_4 = stdev_4./sqrt((sum(~isnan(events_twindows_4))));

else
    time_4=NaN(1,1);  avg_4 = NaN(1,1); sem_4 = NaN(1,1);
end

%Plotting the average reorientation rates with the SEM
subplot(3,2,5)
one = errorline(time_1,avg_1,sem_1); hold on; two =  errorline(time_2,avg_2,sem_2);  hold on; three =   errorline(time_3,avg_3,sem_3);  hold on;  four =  errorline(time_4,avg_4,sem_4);
%properties for plot of input_array1
NameArray_1 = {'LineWidth','Color'};
ValueArray_1 = {0.75,[0 0 1]};
set(one,NameArray_1,ValueArray_1)
set(one(end),'LineWidth',1.5) %setting the last handle, which is the line to be thicker in width
%properties for plot of input_array
NameArray_2 = {'LineWidth','Color'}; 
ValueArray_2 = {0.75,[1 0 0]};
set(two,NameArray_2,ValueArray_2)
set(two(end),'LineWidth',1.5)
%properties for plot of input_array3
NameArray_3 = {'LineWidth','Color'};
ValueArray_3 = {0.75,[0 0.25 1]};
set(three,NameArray_3,ValueArray_3)
set(three(end),'LineWidth',1.5)
%properties for plot of input_array4
NameArray_4 = {'LineWidth','Color'};
ValueArray_4 = {0.75,[0.25 0.25 0.25]};
set(four,NameArray_4,ValueArray_4)
set(four(end),'LineWidth',1.5)
%gca properties
NameArray_gca = {'FontSize','XLim','YLim','FontName'};
ValueArray_gca = {12, [0 maxx],[0 1],'Arial'};
set(gca,NameArray_gca,ValueArray_gca)
Xl = xlabel('Time (min)');
Yl = ylabel('Reorientations');
set(Xl,'FontSize',14,'FontName','Arial')
set(Yl,'FontSize',14,'FontName','Arial')
% labels = {label_1 label_2 label_3 label_4};
% labels = labels(1:length(varargin));
% handle = [one(end) two(end) three(end) four(end)]; %so that it only incldues plot and nopt the fill in the legend
% handle = handle(1:length(varargin));
% l = legend(handle,labels);
% legend('boxoff');
% set(l,'FontSize',9)
title('pure Omega Frequency','FontSize',12)
%% Pure Upsilon Frequencies

%Same as above. Only thing that is different is that in
%events_per_animal_twindows the input for reor_type is 'pure_omega')

time_window = 2; %setting time-window to 2mins. Can change here.


[events_twindows_1,time_1] = events_per_animal_twindows(input_array1,time_window,'pureUpsilon',FrameRate);
avg_1 = nanmean(events_twindows_1);
stdev_1 = nanstd(events_twindows_1);
sem_1 = stdev_1./sqrt((sum(~isnan(events_twindows_1))));


%for arrays 2-4 will only create calculate these things if the array is not
%empty. If it is empty I will make these variables empty
if isempty(input_array2)==0
    [events_twindows_2,time_2] = events_per_animal_twindows(input_array2,time_window,'pureUpsilon',FrameRate);
    avg_2 = nanmean(events_twindows_2);
    stdev_2 = nanstd(events_twindows_2);
    sem_2 = stdev_2./sqrt((sum(~isnan(events_twindows_2))));
else
    time_2=NaN(1,1);  avg_2 = NaN(1,1);  sem_2 = NaN(1,1);
end

if isempty(input_array3)==0
    [events_twindows_3,time_3] = events_per_animal_twindows(input_array3,time_window,'pureUpsilon',FrameRate);
    avg_3 = nanmean(events_twindows_3);
    stdev_3 = nanstd(events_twindows_3);
    sem_3 = stdev_3./sqrt((sum(~isnan(events_twindows_3))));
else
    time_3=NaN(1,1);  avg_3 = NaN(1,1); sem_3 = NaN(1,1);
end

if isempty(input_array4)==0
    [events_twindows_4,time_4] = events_per_animal_twindows(input_array4,time_window,'pureUpsilon',FrameRate);
    avg_4 = nanmean(events_twindows_4);
    stdev_4 = nanstd(events_twindows_4);
    sem_4 = stdev_4./sqrt((sum(~isnan(events_twindows_4))));

else
    time_4=NaN(1,1);  avg_4 = NaN(1,1); sem_4 = NaN(1,1);
end

%Plotting the average reorientation rates with the SEM
subplot(3,2,6)
one = errorline(time_1,avg_1,sem_1); hold on; two =  errorline(time_2,avg_2,sem_2);  hold on; three =   errorline(time_3,avg_3,sem_3);  hold on;  four =  errorline(time_4,avg_4,sem_4);
%properties for plot of input_array1
NameArray_1 = {'LineWidth','Color'};
ValueArray_1 = {0.75,[0 0 1]};
set(one,NameArray_1,ValueArray_1)
set(one(end),'LineWidth',1.5) %setting the last handle, which is the line to be thicker in width
%properties for plot of input_array
NameArray_2 = {'LineWidth','Color'}; 
ValueArray_2 = {0.75,[1 0 0]};
set(two,NameArray_2,ValueArray_2)
set(two(end),'LineWidth',1.5)
%properties for plot of input_array3
NameArray_3 = {'LineWidth','Color'};
ValueArray_3 = {0.75,[0 0.25 1]};
set(three,NameArray_3,ValueArray_3)
set(three(end),'LineWidth',1.5)
%properties for plot of input_array4
NameArray_4 = {'LineWidth','Color'};
ValueArray_4 = {0.75,[0.25 0.25 0.25]};
set(four,NameArray_4,ValueArray_4)
set(four(end),'LineWidth',1.5)
%gca properties
NameArray_gca = {'FontSize','XLim','YLim','FontName'};
ValueArray_gca = {12, [0 maxx],[0 4],'Arial'};
set(gca,NameArray_gca,ValueArray_gca)
Xl = xlabel('Time (min)');
Yl = ylabel('Reorientations');
set(Xl,'FontSize',14,'FontName','Arial')
set(Yl,'FontSize',14,'FontName','Arial')
% labels = {label_1 label_2 label_3 label_4};
% labels = labels(1:length(varargin));
% handle = [one(end) two(end) three(end) four(end)]; %so that it only incldues plot and nopt the fill in the legend
% handle = handle(1:length(varargin));
% l = legend(handle,labels);
% legend('boxoff');
% set(l,'FontSize',9)
title('pure Upsilon Frequency','FontSize',12)
%% figure 6
figure(6);
NameArray_fig = {'PaperUnits','PaperSize', 'PaperPosition', 'Units','Position','PaperOrientation'};
ValueArray_fig = {'inches', [8.5000 11], [0 1 8.5000 9],'inches', [0 1 8.5000 9],'portrait'};
set(figure(6),NameArray_fig,ValueArray_fig)

%% pure_lRev naive
subplot(3,2,1)
%title(fix_title_string(target_odor_linkedTracks1(1).Name));
for m = 1:length(target_odor_linkedTracks1)
    %plotting from attribute matrix
    plot(time1,attribute_matrix1(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    hold on;
    %plotting from origional array
    %plot(((target_odor_linkedTracks1(m).Frames)/(target_odor_linkedTracks1(m).FrameRate)),target_odor_linkedTracks1(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks1)
    stateNo = 4.1;
    stmatrix1 = track_field_to_matrix(target_odor_linkedTracks1,'State');
    st1 = find((abs(stmatrix1(m,:) - stateNo))<0.01);
    plot(st1/3, attribute_matrix1(m,st1), 'x','MarkerSize',2, 'MarkerEdgeColor','b')
    hold on;
    end
    hold on;
plot(time1, mean_dist_1, color1, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:pure sRev(time independent)', 'FontSize', 12);
orient portrait;

%% pure_lRev adapt
subplot(3,2,2)
%title(fix_title_string(target_odor_linkedTracks2(1).Name));
for m = 1:length(target_odor_linkedTracks2)
    %plotting from attribute matrix
    plot(time2,attribute_matrix2(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    %plotting from origional array
    %plot(((target_odor_linkedTracks2(m).Frames)/(target_odor_linkedTracks2(m).FrameRate)),target_odor_linkedTracks2(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks2)
    stateNo = 4.1;
    stmatrix2 = track_field_to_matrix(target_odor_linkedTracks2,'State');
    st2 = find((abs(stmatrix2(m,:) - stateNo))<0.01);
    plot(st2/3, attribute_matrix2(m,st2), 'x','MarkerSize',2,'MarkerEdgeColor','r')
    hold on;
    end
    hold on;
plot(time2, mean_dist_2, color2, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:pure sRev(time independent)', 'FontSize', 12);
orient portrait;
%% pure_sRev naive
subplot(3,2,3)
%title(fix_title_string(target_odor_linkedTracks1(1).Name));
for m = 1:length(target_odor_linkedTracks1)
    %plotting from attribute matrix
    plot(time1,attribute_matrix1(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    hold on;
    %plotting from origional array
    %plot(((target_odor_linkedTracks1(m).Frames)/(target_odor_linkedTracks1(m).FrameRate)),target_odor_linkedTracks1(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks1)
    stateNo = 5.1;
    stmatrix1 = track_field_to_matrix(target_odor_linkedTracks1,'State');
    st1 = find((abs(stmatrix1(m,:) - stateNo))<0.01);
    plot(st1/3, attribute_matrix1(m,st1), 'x','MarkerSize',2, 'MarkerEdgeColor','b')
    hold on;
    end
    hold on;
plot(time1, mean_dist_1, color1, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:pure sRev(time independent)', 'FontSize', 12);
orient portrait;

%% pure_sRev adapt
subplot(3,2,4)
%title(fix_title_string(target_odor_linkedTracks2(1).Name));
for m = 1:length(target_odor_linkedTracks2)
    %plotting from attribute matrix
    plot(time2,attribute_matrix2(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
    set(gca, 'FontSize', 12)
    %plotting from origional array
    %plot(((target_odor_linkedTracks2(m).Frames)/(target_odor_linkedTracks2(m).FrameRate)),target_odor_linkedTracks2(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
    hold on;
    end
    hold on;
for m = 1:length(target_odor_linkedTracks2)
    stateNo = 5.1;
    stmatrix2 = track_field_to_matrix(target_odor_linkedTracks2,'State');
    st2 = find((abs(stmatrix2(m,:) - stateNo))<0.01);
    plot(st2/3, attribute_matrix2(m,st2), 'x','MarkerSize',2,'MarkerEdgeColor','r')
    hold on;
    end
    hold on;
plot(time2, mean_dist_2, color2, 'linewidth',3);
xlim([0 max([time1 time2])]);
if(strcmp(attribute,'odor_distance'))
    ylim([0 90]);
    ylabel('Distance to odor (mm)', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:pure sRev(time independent)', 'FontSize', 12);
orient portrait;
%% Pure lRev Frequencies

%Same as above. Only thing that is different is that in
%events_per_animal_twindows the input for reor_type is 'pure_omega')

time_window = 2; %setting time-window to 2mins. Can change here.


[events_twindows_1,time_1] = events_per_animal_twindows(input_array1,time_window,'pure_lRev',FrameRate);
avg_1 = nanmean(events_twindows_1);
stdev_1 = nanstd(events_twindows_1);
sem_1 = stdev_1./sqrt((sum(~isnan(events_twindows_1))));


%for arrays 2-4 will only create calculate these things if the array is not
%empty. If it is empty I will make these variables empty
if isempty(input_array2)==0
    [events_twindows_2,time_2] = events_per_animal_twindows(input_array2,time_window,'pure_lRev',FrameRate);
    avg_2 = nanmean(events_twindows_2);
    stdev_2 = nanstd(events_twindows_2);
    sem_2 = stdev_2./sqrt((sum(~isnan(events_twindows_2))));
else
    time_2=NaN(1,1);  avg_2 = NaN(1,1);  sem_2 = NaN(1,1);
end

if isempty(input_array3)==0
    [events_twindows_3,time_3] = events_per_animal_twindows(input_array3,time_window,'pure_lRev',FrameRate);
    avg_3 = nanmean(events_twindows_3);
    stdev_3 = nanstd(events_twindows_3);
    sem_3 = stdev_3./sqrt((sum(~isnan(events_twindows_3))));
else
    time_3=NaN(1,1);  avg_3 = NaN(1,1); sem_3 = NaN(1,1);
end

if isempty(input_array4)==0
    [events_twindows_4,time_4] = events_per_animal_twindows(input_array4,time_window,'pure_lRev',FrameRate);
    avg_4 = nanmean(events_twindows_4);
    stdev_4 = nanstd(events_twindows_4);
    sem_4 = stdev_4./sqrt((sum(~isnan(events_twindows_4))));

else
    time_4=NaN(1,1);  avg_4 = NaN(1,1); sem_4 = NaN(1,1);
end

%Plotting the average reorientation rates with the SEM
subplot(3,2,5)
one = errorline(time_1,avg_1,sem_1); hold on; two =  errorline(time_2,avg_2,sem_2);  hold on; three =   errorline(time_3,avg_3,sem_3);  hold on;  four =  errorline(time_4,avg_4,sem_4);
%properties for plot of input_array1
NameArray_1 = {'LineWidth','Color'};
ValueArray_1 = {0.75,[0 0 1]};
set(one,NameArray_1,ValueArray_1)
set(one(end),'LineWidth',1.5) %setting the last handle, which is the line to be thicker in width
%properties for plot of input_array
NameArray_2 = {'LineWidth','Color'}; 
ValueArray_2 = {0.75,[1 0 0]};
set(two,NameArray_2,ValueArray_2)
set(two(end),'LineWidth',1.5)
%properties for plot of input_array3
NameArray_3 = {'LineWidth','Color'};
ValueArray_3 = {0.75,[0 0.25 1]};
set(three,NameArray_3,ValueArray_3)
set(three(end),'LineWidth',1.5)
%properties for plot of input_array4
NameArray_4 = {'LineWidth','Color'};
ValueArray_4 = {0.75,[0.25 0.25 0.25]};
set(four,NameArray_4,ValueArray_4)
set(four(end),'LineWidth',1.5)
%gca properties
NameArray_gca = {'FontSize','XLim','YLim','FontName'};
ValueArray_gca = {12, [0 maxx],[0 3],'Arial'};
set(gca,NameArray_gca,ValueArray_gca)
Xl = xlabel('Time (min)');
Yl = ylabel('Reorientations');
set(Xl,'FontSize',14,'FontName','Arial')
set(Yl,'FontSize',14,'FontName','Arial')
% labels = {label_1 label_2 label_3 label_4};
% labels = labels(1:length(varargin));
% handle = [one(end) two(end) three(end) four(end)]; %so that it only incldues plot and nopt the fill in the legend
% handle = handle(1:length(varargin));
% l = legend(handle,labels);
% legend('boxoff');
% set(l,'FontSize',9)
title('Pure lRev Frequency','FontSize',12)

%% Pure sRev Frequencies

%Same as above. Only thing that is different is that in
%events_per_animal_twindows the input for reor_type is 'pure_omega')

time_window = 2; %setting time-window to 2mins. Can change here.


[events_twindows_1,time_1] = events_per_animal_twindows(input_array1,time_window,'pure_sRev',FrameRate);
avg_1 = nanmean(events_twindows_1);
stdev_1 = nanstd(events_twindows_1);
sem_1 = stdev_1./sqrt((sum(~isnan(events_twindows_1))));


%for arrays 2-4 will only create calculate these things if the array is not
%empty. If it is empty I will make these variables empty
if isempty(input_array2)==0
    [events_twindows_2,time_2] = events_per_animal_twindows(input_array2,time_window,'pure_sRev',FrameRate);
    avg_2 = nanmean(events_twindows_2);
    stdev_2 = nanstd(events_twindows_2);
    sem_2 = stdev_2./sqrt((sum(~isnan(events_twindows_2))));
else
    time_2=NaN(1,1);  avg_2 = NaN(1,1);  sem_2 = NaN(1,1);
end

if isempty(input_array3)==0
    [events_twindows_3,time_3] = events_per_animal_twindows(input_array3,time_window,'pure_sRev',FrameRate);
    avg_3 = nanmean(events_twindows_3);
    stdev_3 = nanstd(events_twindows_3);
    sem_3 = stdev_3./sqrt((sum(~isnan(events_twindows_3))));
else
    time_3=NaN(1,1);  avg_3 = NaN(1,1); sem_3 = NaN(1,1);
end

if isempty(input_array4)==0
    [events_twindows_4,time_4] = events_per_animal_twindows(input_array4,time_window,'pure_sRev',FrameRate);
    avg_4 = nanmean(events_twindows_4);
    stdev_4 = nanstd(events_twindows_4);
    sem_4 = stdev_4./sqrt((sum(~isnan(events_twindows_4))));

else
    time_4=NaN(1,1);  avg_4 = NaN(1,1); sem_4 = NaN(1,1);
end

%Plotting the average reorientation rates with the SEM
subplot(3,2,6)
one = errorline(time_1,avg_1,sem_1); hold on; two =  errorline(time_2,avg_2,sem_2);  hold on; three =   errorline(time_3,avg_3,sem_3);  hold on;  four =  errorline(time_4,avg_4,sem_4);
%properties for plot of input_array1
NameArray_1 = {'LineWidth','Color'};
ValueArray_1 = {0.75,[0 0 1]};
set(one,NameArray_1,ValueArray_1)
set(one(end),'LineWidth',1.5) %setting the last handle, which is the line to be thicker in width
%properties for plot of input_array
NameArray_2 = {'LineWidth','Color'}; 
ValueArray_2 = {0.75,[1 0 0]};
set(two,NameArray_2,ValueArray_2)
set(two(end),'LineWidth',1.5)
%properties for plot of input_array3
NameArray_3 = {'LineWidth','Color'};
ValueArray_3 = {0.75,[0 0.25 1]};
set(three,NameArray_3,ValueArray_3)
set(three(end),'LineWidth',1.5)
%properties for plot of input_array4
NameArray_4 = {'LineWidth','Color'};
ValueArray_4 = {0.75,[0.25 0.25 0.25]};
set(four,NameArray_4,ValueArray_4)
set(four(end),'LineWidth',1.5)
%gca properties
NameArray_gca = {'FontSize','XLim','YLim','FontName'};
ValueArray_gca = {12, [0 maxx],[0 3],'Arial'};
set(gca,NameArray_gca,ValueArray_gca)
Xl = xlabel('Time (min)');
Yl = ylabel('Reorientations');
set(Xl,'FontSize',14,'FontName','Arial')
set(Yl,'FontSize',14,'FontName','Arial')
% labels = {label_1 label_2 label_3 label_4};
% labels = labels(1:length(varargin));
% handle = [one(end) two(end) three(end) four(end)]; %so that it only incldues plot and nopt the fill in the legend
% handle = handle(1:length(varargin));
% l = legend(handle,labels);
% legend('boxoff');
% set(l,'FontSize',9)
title('Pure sRev Frequency','FontSize',12)

%% figure 7
figure(7);
NameArray_fig = {'PaperUnits','PaperSize', 'PaperPosition', 'Units','Position','PaperOrientation'};
ValueArray_fig = {'inches', [8.5000 11], [0 1 8.5000 9],'inches', [0 1 8.5000 9],'portrait'};
set(figure(7),NameArray_fig,ValueArray_fig)

%% Accumulation Naive number
subplot(3,2,1)
%title(fix_title_string(target_odor_linkedTracks1(1).Name));

%load('C:\Users\DuCheng\Documents\MATLAB\Du\standard1cmsemicircle.chemotaxis_regions.mat');%load a standard chemotaxis region
%load('C:\Users\DuCheng\Documents\MATLAB\Du\standard1.2cmsemicircle.chemotaxis_regions.mat');
%load('C:\Users\DuCheng\Documents\MATLAB\Du\standard2cmsemicircle.chemotaxis_regions.mat');
load('C:\Users\DuCheng\Documents\MATLAB\Du\standard3sections.chemotaxis_regions.mat');

SmoothX1 = track_field_to_matrix(target_odor_linkedTracks1, 'SmoothX');
SmoothY1 = track_field_to_matrix(target_odor_linkedTracks1, 'SmoothY');
odor_distance1 = track_field_to_matrix(target_odor_linkedTracks1, 'odor_distance');
odor_angle1 = track_field_to_matrix(target_odor_linkedTracks1, 'odor_angle');
State1 = track_field_to_matrix(target_odor_linkedTracks1, 'State');
Speed1 = track_field_to_matrix(target_odor_linkedTracks1, 'Speed');

[num_targ1, num_control1, targ_finished_tracks, control_finished_tracks,SmoothX1_updated, SmoothY1_updated, avg_dist_target1, avg_dist_control1, std_dev_dist_target1, std_dev_dist_control1] = ...
    number_dist_at_target_twochoice(SmoothX1, SmoothY1, target_point, target_verticies, control_point, control_verticies); 

plot(time1, num_targ1,'LineWidth',2,'Color','b','LineStyle','-')
hold on
plot(time1, num_control1,'LineWidth',2,'Color','b','LineStyle','--')

arrivalcounts1 = [0 diff(num_targ1)];
arrivalframeinds1 = find(arrivalscounts1>0);
arrivalframes1 = [];

for i = 1:length(arrivalframeinds1)
arrivalframes1 = [arrivalframes1; repmat(arrivalframeinds1(i),arrivalcounts1(arrivalframeinds1(i)),1)];
plot(arrivalframes1, 1, 'x','MarkerSize',2, 'MarkerEdgeColor','b')
end

arrivalcounts2 = [0 diff(num_control1)];
arrivalframeinds2 = find(arrivalscounts2>0);
arrivalframes2 = [];
for i = 1:length(arrivalframeinds2)
arrivalframes2 = [arrivalframes2; repmat(arrivalframeinds2(i),arrivalcounts1(arrivalframeinds2(i)),1)];
end

xlim([0 max([time1 time2])]);
ylim([0 max(num_targ1+num_control1)]);
if(strcmp(attribute,'odor_distance'))
    ylabel('Number of Animals', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:Accumulation', 'FontSize', 12);
orient portrait;
%% Accumlation adapt number
subplot(3,2,2)
%title(fix_title_string(target_odor_linkedTracks2(1).Name));

%load('C:\Users\ChengDu\Documents\MATLAB\Du\standard2cmsemicircle.chemotaxis_regions.mat');%load standard chemotaxis region file
%load('C:\Users\ChengDu\Documents\MATLAB\Du\standard1.5cmsemicircle.chemotaxis_regions.mat');

SmoothX2 = track_field_to_matrix(target_odor_linkedTracks2, 'SmoothX');
SmoothY2 = track_field_to_matrix(target_odor_linkedTracks2, 'SmoothY');
odor_distance2 = track_field_to_matrix(target_odor_linkedTracks2, 'odor_distance');
odor_angle2 = track_field_to_matrix(target_odor_linkedTracks2, 'odor_angle');
State2 = track_field_to_matrix(target_odor_linkedTracks2, 'State');
Speed2 = track_field_to_matrix(target_odor_linkedTracks2, 'Speed');

[num_targ2, num_control2, targ_finished_tracks, control_finished_tracks, SmoothX2_updated, SmoothY2_updated, avg_dist_target2, avg_dist_control2, std_dev_dist_target2, std_dev_dist_control2] = ...
    number_dist_at_target_twochoice(SmoothX2, SmoothY2, target_point, target_verticies, control_point, control_verticies); 

plot(time2, num_targ2,'LineWidth',2,'Color','r','LineStyle','-')
hold on
plot(time2, num_control2,'LineWidth',2,'Color','r','LineStyle','--')

xlim([0 max([time1 time2])]);
ylim([0 max(num_targ2+num_control2)]);
if(strcmp(attribute,'odor_distance'))
    ylabel('Number of Animals', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:Accumulation', 'FontSize', 12);
orient portrait;
%% Accumulation Naive Fraction
subplot(3,2,3)
%title(fix_title_string(target_odor_linkedTracks1(1).Name));

numactivetracks1 = sum(~isnan(SmoothX1),1);
if max(numactivetracks1)>(num_targ1(end)+num_control1(end)+numactivetracks1(end));
    totalworm1 = max(numactivetracks1);
else
    totalworm1 = (num_targ1(end)+num_control1(end)+numactivetracks1(end));
end

perctarget1 = num_targ1/totalworm1;
perccontrol1 = num_control1/totalworm1;
%end point chemotaxis index
Chemotaxis_index1 = (num_targ1(end)-num_control1(end))/totalworm1;
CIlabel1 = strjoin({'CI =',num2str(Chemotaxis_index1,'%.2f')});

plot(time1, perctarget1,'LineWidth',2,'Color','b','LineStyle','-')
hold on
plot(time1, perccontrol1,'LineWidth',2,'Color','b','LineStyle','--')
text(100,0.5, CIlabel1,'FontSize', 9)

xlim([0 max([time1 time2])]);
ylim([0 1]);
if(strcmp(attribute,'odor_distance'))
    ylabel('Fraction of Animals', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:Accumulation', 'FontSize', 12);
orient portrait;
%% Accumulation Adapt Fraction
subplot(3,2,4)
%title(fix_title_string(target_odor_linkedTracks2(1).Name));

numactivetracks2 = sum(~isnan(SmoothX2),1);

if max(numactivetracks2)>(num_targ2(end)+num_control2(end)+numactivetracks2(end));
    totalworm2 = max(numactivetracks2);
else
    totalworm2 = (num_targ2(end)+num_control2(end)+numactivetracks2(end));
end

perctarget2 = num_targ2/totalworm2;
perccontrol2 = num_control2/totalworm2;

%end point chemotaxis index
Chemotaxis_index2 = (num_targ2(end)-num_control2(end))/totalworm2;
CIlabel2 = strjoin({'CI =',num2str(Chemotaxis_index2,'%.2f')});

plot(time2, perctarget2,'LineWidth',2,'Color','r','LineStyle','-')
hold on
plot(time2, perccontrol2,'LineWidth',2,'Color','r','LineStyle','--')
text(100,0.5, CIlabel2,'FontSize', 9)

xlim([0 max([time1 time2])]);
ylim([0 1]);
if(strcmp(attribute,'odor_distance'))
    ylabel('Fraction of Animals', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end
xlabel('Time (sec), x:Accumulation', 'FontSize', 12);
orient portrait;

%% Accumulation Naive Fraction inverted
subplot(3,2,5)
%title(fix_title_string(target_odor_linkedTracks1(1).Name));
% 
% numactivetracks1 = sum(~isnan(SmoothX1),1);
% 
% if max(numactivetracks1)>(num_targ1(end)+num_control1(end));
%     totalworm1 = max(numactivetracks1);
% else
%     totalworm1 = (num_targ1(end)+num_control1(end));
% end
% 
% perctarget1 = num_targ1/totalworm1;
% perccontrol1 = num_control1/totalworm1;
perccontrolinv1 = 1 - perccontrol1;

% plot(time1, perctarget1,'LineWidth',2,'Color','b','LineStyle','-')
% hold on
% plot(time1, perccontrol1,'LineWidth',2,'Color','b','LineStyle','--')

yyaxis left
plot(time1, perctarget1,'LineWidth',2,'Color','b','LineStyle','-')
text(100,0.5, CIlabel1,'FontSize', 9)

xlim([0 max([time1 time2])]);
ylim([0 1]);
if(strcmp(attribute,'odor_distance'))
    ylabel('Fraction of Animals', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end

yyaxis right
plot(time1, perccontrolinv1,'LineWidth',2,'Color','b','LineStyle','--')

yticks([0 0.2 0.4 0.6 0.8 1])
yticklabels({'1' '0.8' '0.6' '0.4' '0.2' '0'})
xlim([0 max([time1 time2])]);
ylim([0 1]);
if(strcmp(attribute,'odor_distance'))
    ylabel('Fraction of Animals', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end

xlabel('Time (sec), x:Accumulation', 'FontSize', 12);
orient portrait;
%% Accumulation Adapt Fraction inverted
subplot(3,2,6)
%title(fix_title_string(target_odor_linkedTracks2(1).Name));

numactivetracks2 = sum(~isnan(SmoothX2),1);

% perctarget2 = num_targ2/max(numactivetracks2);
% perccontrol2 = num_control2/max(numactivetracks2);
perccontrolinv2 = 1 - perccontrol2;

% plot(time2, perctarget2,'LineWidth',2,'Color','r','LineStyle','-')
% hold on
% plot(time2, perccontrol2,'LineWidth',2,'Color','r','LineStyle','--')

yyaxis left
plot(time2, perctarget2,'LineWidth',2,'Color','r','LineStyle','-')
text(100,0.5, CIlabel2,'FontSize', 9)

xlim([0 max([time1 time2])]);
ylim([0 1]);
if(strcmp(attribute,'odor_distance'))
    ylabel('Fraction of Animals', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end

yyaxis right
plot(time2, perccontrolinv2,'LineWidth',2,'Color','r','LineStyle','--')

yticks([0 0.2 0.4 0.6 0.8 1])
yticklabels({'1' '0.8' '0.6' '0.4' '0.2' '0'})
xlim([0 max([time1 time2])]);
ylim([0 1]);
if(strcmp(attribute,'odor_distance'))
    ylabel('Fraction of Animals', 'FontSize', 12);
else
    ylabel(fix_title_string(attribute));
end

xlabel('Time (sec), x:Accumulation', 'FontSize', 12);
orient portrait;


% %% figure 6
% figure(6);
% NameArray_fig = {'PaperUnits','PaperSize', 'PaperPosition', 'Units','Position','PaperOrientation'};
% ValueArray_fig = {'inches', [8.5000 11], [0 1 8.5000 9],'inches', [0 1 8.5000 9],'portrait'};
% set(figure(6),NameArray_fig,ValueArray_fig)
% %% pause naive
% 
% subplot(3,2,3)
% title(fix_title_string(target_odor_linkedTracks1(1).Name));
% for m = 1:length(target_odor_linkedTracks1)
%     %plotting from attribute matrix
%     plot(time1,attribute_matrix1(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
%     set(gca, 'FontSize', 12)
%     hold on;
%     %plotting from origional array
%     %plot(((target_odor_linkedTracks1(m).Frames)/(target_odor_linkedTracks1(m).FrameRate)),target_odor_linkedTracks1(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
%     end
%     hold on;
% for m = 1:length(target_odor_linkedTracks1)
%     stateNo = 1.1;
%     stmatrix1 = track_field_to_matrix(target_odor_linkedTracks1,'State');
%     st1 = find((abs(stmatrix1(m,:) - stateNo))<0.01);
%     plot(st1/3, attribute_matrix1(m,st1), 'x','MarkerSize',2, 'MarkerEdgeColor','b')
%     hold on;
%     end
%     hold on;
% plot(time1, mean_dist_1, color1, 'linewidth',3);
% xlim([0 max([time1 time2])]);
% if(strcmp(attribute,'odor_distance'))
%     ylim([0 90]);
%     ylabel('Distance to odor (mm)', 'FontSize', 12);
% else
%     ylabel(fix_title_string(attribute));
% end
% xlabel('Time (sec), x:Pause', 'FontSize', 12);
% orient portrait;
% 
% %% pause adapt
% subplot(3,2,4)
% title(fix_title_string(target_odor_linkedTracks2(1).Name));
% for m = 1:length(target_odor_linkedTracks2)
%     %plotting from attribute matrix
%     plot(time2,attribute_matrix2(m,:), 'Color',[0.7,0.7,0.7],'linewidth',0.5);
%     set(gca, 'FontSize', 12)
%     %plotting from origional array
%     %plot(((target_odor_linkedTracks2(m).Frames)/(target_odor_linkedTracks2(m).FrameRate)),target_odor_linkedTracks2(m).(attribute), 'Color',[0.7,0.7,0.7],'linewidth',0.5); hold on;
%     hold on;
%     end
%     hold on;
% for m = 1:length(target_odor_linkedTracks2)
%     stateNo = 1.1;
%     stmatrix2 = track_field_to_matrix(target_odor_linkedTracks2,'State');
%     st2 = find((abs(stmatrix2(m,:) - stateNo))<0.01);
%     plot(st2/3, attribute_matrix2(m,st2), 'x','MarkerSize',2,'MarkerEdgeColor','r')
%     hold on;
%     end
%     hold on;
% plot(time2, mean_dist_2, color2, 'linewidth',3);
% xlim([0 max([time1 time2])]);
% if(strcmp(attribute,'odor_distance'))
%     ylim([0 90]);
%     ylabel('Distance to odor (mm)', 'FontSize', 12);
% else
%     ylabel(fix_title_string(attribute));
% end
% xlabel('Time (sec), x:Pause', 'FontSize', 12);
% orient portrait;

%%
figure (8)
NameArray_fig = {'PaperUnits','PaperSize', 'PaperPosition', 'Units','Position','PaperOrientation'};
ValueArray_fig = {'inches', [8.5000 11], [0 1 8.5000 9],'inches', [0 1 8.5000 9],'portrait'};
set(figure(8),NameArray_fig,ValueArray_fig)

%%
subplot(4,3,1)

title('LRevOmega naive');
[event_matrix,time_matrix] = mark_reor_events(target_odor_linkedTracks1,'lRevOmega',3); %this puts a 1 in a new matrix everywhere where the behavior occurs
event_matrix = event_matrix';

% find the indices in the event_matrix where the behavior occurred
[row,col] = find(event_matrix == 1);

% get the frames to check
newframes = NaN(size(row,1),2);

for idx = 1:size(row,1)
    if col(idx)<10
        newidx = 1;
    else
        newidx = col(idx)-9;
    end
    newframes(idx,:) = [row(idx) newidx];
end

% get the odor angles at these frames

odor_angles1 = track_field_to_matrix(target_odor_linkedTracks1,'odor_angle');
anglesofinterest1 = diag(odor_angles1(newframes(:,1),newframes(:,2)));
hist(anglesofinterest1,6);
ylabel('Events');
xlabel('Angles Before LRevOmega');
xlim([0 180]);

%%
subplot(4,3,2)

title('LRevOmega adapt');
[event_matrix,time_matrix] = mark_reor_events(target_odor_linkedTracks2,'lRevOmega',3); %this puts a 1 in a new matrix everywhere where the behavior occurs
event_matrix = event_matrix';

% find the indices in the event_matrix where the behavior occurred
[row,col] = find(event_matrix == 1);

% get the frames to check
newframes = NaN(size(row,1),2);

for idx = 1:size(row,1)
    if col(idx)<10
        newidx = 1;
    else
        newidx = col(idx)-9;
    end
    newframes(idx,:) = [row(idx) newidx];
end

% get the odor angles at these frames

odor_angles2 = track_field_to_matrix(target_odor_linkedTracks2,'odor_angle');
anglesofinterest2 = diag(odor_angles2(newframes(:,1),newframes(:,2)));
hist(anglesofinterest2,6);
ylabel('Events');
xlabel('Angles Before LRevOmega');
xlim([0 180]);

%%
subplot(4,3,3)

angleslessthan901 = anglesofinterest1(anglesofinterest1<90);
anglesmorethan901 = anglesofinterest1(anglesofinterest1>90);
angleslessthan902 = anglesofinterest2(anglesofinterest2<90);
anglesmorethan902 = anglesofinterest2(anglesofinterest2>90);

countangleslessthan901 = numel(angleslessthan901);
countanglesmorethan901 = numel(anglesmorethan901);
countangleslessthan902 = numel(angleslessthan902);
countanglesmorethan902 = numel(anglesmorethan902);

freqangleslessthan901 = countangleslessthan901/totalworm1;
freqanglesmorethan901 = countanglesmorethan901/totalworm1;
freqangleslessthan902 = countangleslessthan902/totalworm2;
freqanglesmorethan902 = countanglesmorethan902/totalworm2;

anglecount = [freqangleslessthan901 freqanglesmorethan901; freqangleslessthan902 freqanglesmorethan902];

anglebarplot = bar(anglecount);
ylabel('turnfrequency');
xlabel('LRevOmega');
%xlim([0 180]);
anglebarplot(2).FaceColor = 'red';

%%
subplot(4,3,4)

title('LRevUpsilon naive');
[event_matrix,time_matrix] = mark_reor_events(target_odor_linkedTracks1,'lRevUpsilon',3); %this puts a 1 in a new matrix everywhere where the behavior occurs
event_matrix = event_matrix';

% find the indices in the event_matrix where the behavior occurred
[row,col] = find(event_matrix == 1);

% get the frames to check
newframes = NaN(size(row,1),2);

for idx = 1:size(row,1)
    if col(idx)<10
        newidx = 1;
    else
        newidx = col(idx)-9;
    end
    newframes(idx,:) = [row(idx) newidx];
end

% get the odor angles at these frames

odor_angles1 = track_field_to_matrix(target_odor_linkedTracks1,'odor_angle');
anglesofinterest1 = diag(odor_angles1(newframes(:,1),newframes(:,2)));
hist(anglesofinterest1,6);
ylabel('Events');
xlabel('Angles Before LRevUpsilon');
xlim([0 180]);

%%
subplot(4,3,5)

title('LRevUpsilon adapt');
[event_matrix,time_matrix] = mark_reor_events(target_odor_linkedTracks2,'lRevUpsilon',3); %this puts a 1 in a new matrix everywhere where the behavior occurs
event_matrix = event_matrix';

% find the indices in the event_matrix where the behavior occurred
[row,col] = find(event_matrix == 1);

% get the frames to check
newframes = NaN(size(row,1),2);

for idx = 1:size(row,1)
    if col(idx)<10
        newidx = 1;
    else
        newidx = col(idx)-9;
    end
    newframes(idx,:) = [row(idx) newidx];
end

% get the odor angles at these frames

odor_angles2 = track_field_to_matrix(target_odor_linkedTracks2,'odor_angle');
anglesofinterest2 = diag(odor_angles2(newframes(:,1),newframes(:,2)));
hist(anglesofinterest2,6);
ylabel('Events');
xlabel('Angles Before LRevUpsilon');
xlim([0 180]);

%%
subplot(4,3,6)

angleslessthan901 = anglesofinterest1(anglesofinterest1<90);
anglesmorethan901 = anglesofinterest1(anglesofinterest1>90);
angleslessthan902 = anglesofinterest2(anglesofinterest2<90);
anglesmorethan902 = anglesofinterest2(anglesofinterest2>90);

countangleslessthan901 = numel(angleslessthan901);
countanglesmorethan901 = numel(anglesmorethan901);
countangleslessthan902 = numel(angleslessthan902);
countanglesmorethan902 = numel(anglesmorethan902);

freqangleslessthan901 = countangleslessthan901/totalworm1;
freqanglesmorethan901 = countanglesmorethan901/totalworm1;
freqangleslessthan902 = countangleslessthan902/totalworm2;
freqanglesmorethan902 = countanglesmorethan902/totalworm2;

anglecount = [freqangleslessthan901 freqanglesmorethan901; freqangleslessthan902 freqanglesmorethan902];

anglebarplot = bar(anglecount);
ylabel('turnfrequency');
xlabel('LRevUpsilon');
%xlim([0 180]);
anglebarplot(2).FaceColor = 'red';

%%
subplot(4,3,7)

title('Pure Omega naive');
[event_matrix,time_matrix] = mark_reor_events(target_odor_linkedTracks1,'pure_omega',3); %this puts a 1 in a new matrix everywhere where the behavior occurs
event_matrix = event_matrix';

% find the indices in the event_matrix where the behavior occurred
[row,col] = find(event_matrix == 1);

% get the frames to check
newframes = NaN(size(row,1),2);

for idx = 1:size(row,1)
    if col(idx)<10
        newidx = 1;
    else
        newidx = col(idx)-9;
    end
    newframes(idx,:) = [row(idx) newidx];
end

% get the odor angles at these frames

odor_angles1 = track_field_to_matrix(target_odor_linkedTracks1,'odor_angle');
anglesofinterest1 = diag(odor_angles1(newframes(:,1),newframes(:,2)));
hist(anglesofinterest1,6);
ylabel('Events');
xlabel('Angles Before Pure Omega');
xlim([0 180]);

%%
subplot(4,3,8)

title('Pure Omega adapt');
[event_matrix,time_matrix] = mark_reor_events(target_odor_linkedTracks2,'pure_omega',3); %this puts a 1 in a new matrix everywhere where the behavior occurs
event_matrix = event_matrix';

% find the indices in the event_matrix where the behavior occurred
[row,col] = find(event_matrix == 1);

% get the frames to check
newframes = NaN(size(row,1),2);

for idx = 1:size(row,1)
    if col(idx)<10
        newidx = 1;
    else
        newidx = col(idx)-9;
    end
    newframes(idx,:) = [row(idx) newidx];
end

% get the odor angles at these frames

odor_angles2 = track_field_to_matrix(target_odor_linkedTracks2,'odor_angle');
anglesofinterest2 = diag(odor_angles2(newframes(:,1),newframes(:,2)));
hist(anglesofinterest2,6);
ylabel('Events');
xlabel('Angles Before Pure Omega');
xlim([0 180]);

%%
subplot(4,3,9)

angleslessthan901 = anglesofinterest1(anglesofinterest1<90);
anglesmorethan901 = anglesofinterest1(anglesofinterest1>90);
angleslessthan902 = anglesofinterest2(anglesofinterest2<90);
anglesmorethan902 = anglesofinterest2(anglesofinterest2>90);

countangleslessthan901 = numel(angleslessthan901);
countanglesmorethan901 = numel(anglesmorethan901);
countangleslessthan902 = numel(angleslessthan902);
countanglesmorethan902 = numel(anglesmorethan902);

freqangleslessthan901 = countangleslessthan901/totalworm1;
freqanglesmorethan901 = countanglesmorethan901/totalworm1;
freqangleslessthan902 = countangleslessthan902/totalworm2;
freqanglesmorethan902 = countanglesmorethan902/totalworm2;

anglecount = [freqangleslessthan901 freqanglesmorethan901; freqangleslessthan902 freqanglesmorethan902];

anglebarplot = bar(anglecount);
ylabel('turnfrequency');
xlabel('PureOmega');
%xlim([0 180]);
anglebarplot(2).FaceColor = 'red';


%%
subplot(4,3,10)

title('Pure Upsilon naive');
[event_matrix,time_matrix] = mark_reor_events(target_odor_linkedTracks1,'pure_Upsilon',3); %this puts a 1 in a new matrix everywhere where the behavior occurs
event_matrix = event_matrix';

% find the indices in the event_matrix where the behavior occurred
[row,col] = find(event_matrix == 1);

% get the frames to check
newframes = NaN(size(row,1),2);

for idx = 1:size(row,1)
    if col(idx)<10
        newidx = 1;
    else
        newidx = col(idx)-9;
    end
    newframes(idx,:) = [row(idx) newidx];
end

% get the odor angles at these frames

odor_angles1 = track_field_to_matrix(target_odor_linkedTracks1,'odor_angle');
anglesofinterest1 = diag(odor_angles1(newframes(:,1),newframes(:,2)));
hist(anglesofinterest1,6);
ylabel('Events');
xlabel('Angles Before Pure Upsilon');
xlim([0 180]);

%%
subplot(4,3,11)

title('Pure Upsilon adapt');
[event_matrix,time_matrix] = mark_reor_events(target_odor_linkedTracks2,'pure_Upsilon',3); %this puts a 1 in a new matrix everywhere where the behavior occurs
event_matrix = event_matrix';

% find the indices in the event_matrix where the behavior occurred
[row,col] = find(event_matrix == 1);

% get the frames to check
newframes = NaN(size(row,1),2);

for idx = 1:size(row,1)
    if col(idx)<10
        newidx = 1;
    else
        newidx = col(idx)-9;
    end
    newframes(idx,:) = [row(idx) newidx];
end

% get the odor angles at these frames

odor_angles2 = track_field_to_matrix(target_odor_linkedTracks2,'odor_angle');
anglesofinterest2 = diag(odor_angles2(newframes(:,1),newframes(:,2)));
hist(anglesofinterest2,6);
ylabel('Events');
xlabel('Angles Before Pure Upsilon');
xlim([0 180]);

%%
subplot(4,3,12)

angleslessthan901 = anglesofinterest1(anglesofinterest1<90);
anglesmorethan901 = anglesofinterest1(anglesofinterest1>90);
angleslessthan902 = anglesofinterest2(anglesofinterest2<90);
anglesmorethan902 = anglesofinterest2(anglesofinterest2>90);

countangleslessthan901 = numel(angleslessthan901);
countanglesmorethan901 = numel(anglesmorethan901);
countangleslessthan902 = numel(angleslessthan902);
countanglesmorethan902 = numel(anglesmorethan902);

freqangleslessthan901 = countangleslessthan901/totalworm1;
freqanglesmorethan901 = countanglesmorethan901/totalworm1;
freqangleslessthan902 = countangleslessthan902/totalworm2;
freqanglesmorethan902 = countanglesmorethan902/totalworm2;

anglecount = [freqangleslessthan901 freqanglesmorethan901; freqangleslessthan902 freqanglesmorethan902];

anglebarplot = bar(anglecount);
ylabel('turnfrequency');
xlabel('PureUpsilon');
%xlim([0 180]);
anglebarplot(2).FaceColor = 'red';


%%
figure (9)
NameArray_fig = {'PaperUnits','PaperSize', 'PaperPosition', 'Units','Position','PaperOrientation'};
ValueArray_fig = {'inches', [8.5000 11], [0 1 8.5000 9],'inches', [0 1 8.5000 9],'portrait'};
set(figure(9),NameArray_fig,ValueArray_fig)

%%
subplot(4,3,1)

title('sRevOmega naive');
[event_matrix,time_matrix] = mark_reor_events(target_odor_linkedTracks1,'sRevOmega',3); %this puts a 1 in a new matrix everywhere where the behavior occurs
event_matrix = event_matrix';

% find the indices in the event_matrix where the behavior occurred
[row,col] = find(event_matrix == 1);

% get the frames to check
newframes = NaN(size(row,1),2);

for idx = 1:size(row,1)
    if col(idx)<10
        newidx = 1;
    else
        newidx = col(idx)-9;
    end
    newframes(idx,:) = [row(idx) newidx];
end

% get the odor angles at these frames

odor_angles1 = track_field_to_matrix(target_odor_linkedTracks1,'odor_angle');
anglesofinterest1 = diag(odor_angles1(newframes(:,1),newframes(:,2)));
hist(anglesofinterest1,6);
ylabel('Events');
xlabel('Angles Before sRevOmega');
xlim([0 180]);

%%
subplot(4,3,2)

title('sRevOmega adapt');
[event_matrix,time_matrix] = mark_reor_events(target_odor_linkedTracks2,'sRevOmega',3); %this puts a 1 in a new matrix everywhere where the behavior occurs
event_matrix = event_matrix';

% find the indices in the event_matrix where the behavior occurred
[row,col] = find(event_matrix == 1);

% get the frames to check
newframes = NaN(size(row,1),2);

for idx = 1:size(row,1)
    if col(idx)<10
        newidx = 1;
    else
        newidx = col(idx)-9;
    end
    newframes(idx,:) = [row(idx) newidx];
end

% get the odor angles at these frames

odor_angles2 = track_field_to_matrix(target_odor_linkedTracks2,'odor_angle');
anglesofinterest2 = diag(odor_angles2(newframes(:,1),newframes(:,2)));
hist(anglesofinterest2,6);
ylabel('Events');
xlabel('Angles Before sRevOmega');
xlim([0 180]);

%%
subplot(4,3,3)

angleslessthan901 = anglesofinterest1(anglesofinterest1<90);
anglesmorethan901 = anglesofinterest1(anglesofinterest1>90);
angleslessthan902 = anglesofinterest2(anglesofinterest2<90);
anglesmorethan902 = anglesofinterest2(anglesofinterest2>90);

countangleslessthan901 = numel(angleslessthan901);
countanglesmorethan901 = numel(anglesmorethan901);
countangleslessthan902 = numel(angleslessthan902);
countanglesmorethan902 = numel(anglesmorethan902);

freqangleslessthan901 = countangleslessthan901/totalworm1;
freqanglesmorethan901 = countanglesmorethan901/totalworm1;
freqangleslessthan902 = countangleslessthan902/totalworm2;
freqanglesmorethan902 = countanglesmorethan902/totalworm2;

anglecount = [freqangleslessthan901 freqanglesmorethan901; freqangleslessthan902 freqanglesmorethan902];

anglebarplot = bar(anglecount);
ylabel('turnfrequency');
xlabel('sRevOmega');
%xlim([0 180]);
anglebarplot(2).FaceColor = 'red';


%%
subplot(4,3,4)

title('sRevUpsilon naive');
[event_matrix,time_matrix] = mark_reor_events(target_odor_linkedTracks1,'sRevUpsilon',3); %this puts a 1 in a new matrix everywhere where the behavior occurs
event_matrix = event_matrix';

% find the indices in the event_matrix where the behavior occurred
[row,col] = find(event_matrix == 1);

% get the frames to check
newframes = NaN(size(row,1),2);

for idx = 1:size(row,1)
    if col(idx)<10
        newidx = 1;
    else
        newidx = col(idx)-9;
    end
    newframes(idx,:) = [row(idx) newidx];
end

% get the odor angles at these frames

odor_angles1 = track_field_to_matrix(target_odor_linkedTracks1,'odor_angle');
anglesofinterest1 = diag(odor_angles1(newframes(:,1),newframes(:,2)));
hist(anglesofinterest1,6);
ylabel('Events');
xlabel('Angles Before sRevOmega');
xlim([0 180]);

%%
subplot(4,3,5)

title('sRevUpsilon adapt');
[event_matrix,time_matrix] = mark_reor_events(target_odor_linkedTracks2,'sRevUpsilon',3); %this puts a 1 in a new matrix everywhere where the behavior occurs
event_matrix = event_matrix';

% find the indices in the event_matrix where the behavior occurred
[row,col] = find(event_matrix == 1);

% get the frames to check
newframes = NaN(size(row,1),2);

for idx = 1:size(row,1)
    if col(idx)<10
        newidx = 1;
    else
        newidx = col(idx)-9;
    end
    newframes(idx,:) = [row(idx) newidx];
end

% get the odor angles at these frames

odor_angles2 = track_field_to_matrix(target_odor_linkedTracks2,'odor_angle');
anglesofinterest2 = diag(odor_angles2(newframes(:,1),newframes(:,2)));
hist(anglesofinterest2,6);
ylabel('Events');
xlabel('Angles Before sRevUpsilon');
xlim([0 180]);

%%
subplot(4,3,6)

angleslessthan901 = anglesofinterest1(anglesofinterest1<90);
anglesmorethan901 = anglesofinterest1(anglesofinterest1>90);
angleslessthan902 = anglesofinterest2(anglesofinterest2<90);
anglesmorethan902 = anglesofinterest2(anglesofinterest2>90);

countangleslessthan901 = numel(angleslessthan901);
countanglesmorethan901 = numel(anglesmorethan901);
countangleslessthan902 = numel(angleslessthan902);
countanglesmorethan902 = numel(anglesmorethan902);

freqangleslessthan901 = countangleslessthan901/totalworm1;
freqanglesmorethan901 = countanglesmorethan901/totalworm1;
freqangleslessthan902 = countangleslessthan902/totalworm2;
freqanglesmorethan902 = countanglesmorethan902/totalworm2;

anglecount = [freqangleslessthan901 freqanglesmorethan901; freqangleslessthan902 freqanglesmorethan902];

anglebarplot = bar(anglecount);
ylabel('turnfrequency');
xlabel('sRevUpsilon');
%xlim([0 180]);
anglebarplot(2).FaceColor = 'red';


%%
subplot(4,3,7)

title('pause naive');
[event_matrix,time_matrix] = mark_reor_events(target_odor_linkedTracks1,'pause',3); %this puts a 1 in a new matrix everywhere where the behavior occurs
event_matrix = event_matrix';

% find the indices in the event_matrix where the behavior occurred
[row,col] = find(event_matrix == 1);

% get the frames to check
newframes = NaN(size(row,1),2);

for idx = 1:size(row,1)
    if col(idx)<10
        newidx = 1;
    else
        newidx = col(idx)-9;
    end
    newframes(idx,:) = [row(idx) newidx];
end

% get the odor angles at these frames

odor_angles1 = track_field_to_matrix(target_odor_linkedTracks1,'odor_angle');
anglesofinterest1 = diag(odor_angles1(newframes(:,1),newframes(:,2)));
hist(anglesofinterest1,6);
ylabel('Events');
xlabel('Angles Before pause');
xlim([0 180]);

%%
subplot(4,3,8)

title('pause adapt');
[event_matrix,time_matrix] = mark_reor_events(target_odor_linkedTracks2,'pause',3); %this puts a 1 in a new matrix everywhere where the behavior occurs
event_matrix = event_matrix';

% find the indices in the event_matrix where the behavior occurred
[row,col] = find(event_matrix == 1);

% get the frames to check
newframes = NaN(size(row,1),2);

for idx = 1:size(row,1)
    if col(idx)<10
        newidx = 1;
    else
        newidx = col(idx)-9;
    end
    newframes(idx,:) = [row(idx) newidx];
end

% get the odor angles at these frames

odor_angles2 = track_field_to_matrix(target_odor_linkedTracks2,'odor_angle');
anglesofinterest2 = diag(odor_angles2(newframes(:,1),newframes(:,2)));
hist(anglesofinterest2,6);
ylabel('Events');
xlabel('Angles Before pause');
xlim([0 180]);

%%
subplot(4,3,9)

angleslessthan901 = anglesofinterest1(anglesofinterest1<90);
anglesmorethan901 = anglesofinterest1(anglesofinterest1>90);
angleslessthan902 = anglesofinterest2(anglesofinterest2<90);
anglesmorethan902 = anglesofinterest2(anglesofinterest2>90);

countangleslessthan901 = numel(angleslessthan901);
countanglesmorethan901 = numel(anglesmorethan901);
countangleslessthan902 = numel(angleslessthan902);
countanglesmorethan902 = numel(anglesmorethan902);

freqangleslessthan901 = countangleslessthan901/totalworm1;
freqanglesmorethan901 = countanglesmorethan901/totalworm1;
freqangleslessthan902 = countangleslessthan902/totalworm2;
freqanglesmorethan902 = countanglesmorethan902/totalworm2;

anglecount = [freqangleslessthan901 freqanglesmorethan901; freqangleslessthan902 freqanglesmorethan902];

anglebarplot = bar(anglecount);
ylabel('turnfrequency');
xlabel('pause');
%xlim([0 180]);
anglebarplot(2).FaceColor = 'red';

%%
subplot(4,3,10)

title('fwd naive');
[event_matrix,time_matrix] = mark_reor_events(target_odor_linkedTracks1,'fwd',3); %this puts a 1 in a new matrix everywhere where the behavior occurs
event_matrix = event_matrix';

% find the indices in the event_matrix where the behavior occurred
[row,col] = find(event_matrix == 1);

% get the frames to check
newframes = NaN(size(row,1),2);

for idx = 1:size(row,1)
    if col(idx)<10
        newidx = 1;
    else
        newidx = col(idx)-9;
    end
    newframes(idx,:) = [row(idx) newidx];
end

% get the odor angles at these frames

odor_angles1 = track_field_to_matrix(target_odor_linkedTracks1,'odor_angle');
anglesofinterest1 = diag(odor_angles1(newframes(:,1),newframes(:,2)));
hist(anglesofinterest1,6);
ylabel('Events');
xlabel('Angles Before fwd');
xlim([0 180]);

%%
subplot(4,3,11)

title('fwd adapt');
[event_matrix,time_matrix] = mark_reor_events(target_odor_linkedTracks2,'fwd',3); %this puts a 1 in a new matrix everywhere where the behavior occurs
event_matrix = event_matrix';

% find the indices in the event_matrix where the behavior occurred
[row,col] = find(event_matrix == 1);

% get the frames to check
newframes = NaN(size(row,1),2);

for idx = 1:size(row,1)
    if col(idx)<10
        newidx = 1;
    else
        newidx = col(idx)-9;
    end
    newframes(idx,:) = [row(idx) newidx];
end

% get the odor angles at these frames

odor_angles2 = track_field_to_matrix(target_odor_linkedTracks2,'odor_angle');
anglesofinterest2 = diag(odor_angles2(newframes(:,1),newframes(:,2)));
hist(anglesofinterest2,6);
ylabel('Events');
xlabel('Angles Before fwd');
xlim([0 180]);
%%
subplot(4,3,12)

angleslessthan901 = anglesofinterest1(anglesofinterest1<90);
anglesmorethan901 = anglesofinterest1(anglesofinterest1>90);
angleslessthan902 = anglesofinterest2(anglesofinterest2<90);
anglesmorethan902 = anglesofinterest2(anglesofinterest2>90);

countangleslessthan901 = numel(angleslessthan901);
countanglesmorethan901 = numel(anglesmorethan901);
countangleslessthan902 = numel(angleslessthan902);
countanglesmorethan902 = numel(anglesmorethan902);

freqangleslessthan901 = countangleslessthan901/totalworm1;
freqanglesmorethan901 = countanglesmorethan901/totalworm1;
freqangleslessthan902 = countangleslessthan902/totalworm2;
freqanglesmorethan902 = countanglesmorethan902/totalworm2;

anglecount = [freqangleslessthan901 freqanglesmorethan901; freqangleslessthan902 freqanglesmorethan902];

anglebarplot = bar(anglecount);
ylabel('turnfrequency');
xlabel('fwd');
%xlim([0 180]);
anglebarplot(2).FaceColor = 'red';

%%
save_pdf([1 2 3 4 5 6 7 8 9],filename);

close all

return;
end

%function on converting the array to matrix, fill the matrix with the last
%non-nan value to the end.
function output_matrix = final_attribute_matrix(Tracks, threshold_low, threshold_high, attribute)

output_matrix = track_field_to_matrix(Tracks,attribute); % odor_distance odor_angle

%"~" invert whatever the commond you put in. "isnan" find the first value
%that is nan
for(i=1:length(Tracks))
    track_attribute = Tracks(i).(attribute);
    non_nan_att_ind = find(~isnan(track_attribute), 1, 'last' );
    last_non_nan_att = track_attribute(non_nan_att_ind);
    if(last_non_nan_att > threshold_high)
        output_matrix(i,Tracks(i).Frames(non_nan_att_ind):size(output_matrix,2)) = last_non_nan_att;
    elseif(last_non_nan_att < threshold_low)
        output_matrix(i,Tracks(i).Frames(non_nan_att_ind):size(output_matrix,2)) = last_non_nan_att;
    end
end

%this part is what didn't work in the earlier versions. We thought simply
%take the end value of each track and fill it to the rest would work,
%however, most track ends with Nan instead of a value.

%output_matrix = matrix_replace(output_matrix, '>', threshold_high,threshold_high);
%output_matrix = matrix_replace(output_matrix, '<', threshold_low,threshold_low);

% for m = 1:size(output_matrix,2)
%     indices = find(isnan(output_matrix(:,m)) == 0);
%     if isempty(indices) == 1
%         
%     elseif indices(end) < size(output_matrix,1)
%         if output_matrix(indices(end),m) < threshold_low  ||  output_matrix(indices(end),m) > threshold_high
%             output_matrix(indices(end)+1:end,m) = output_matrix(indices(end),m);
%         end
%     end
% end 

return;    
end