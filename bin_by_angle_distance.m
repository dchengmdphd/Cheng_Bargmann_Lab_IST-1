%script to make binned reorientations by odor bearing angle and odor
%distance.

%% 1. open one of the merged "GOOD.mat" linkedTracks files that I make
[filename, pathname, ~] = uigetfile({'*.mat'});
tracks = load([pathname filename ]);
tmp = fieldnames(tracks);
tmp = tmp{1};
tracks = eval(['tracks.' tmp]);
%% select reor type
% reor_type = 'pure_upsilon';
reor_type = 'allnonUpsilon';
% reor_type = 'lRevOmega';
%% 2. Extract worm frames at different bearing angles
% extract tracks and frame numbers for different ranges of bearing angles
angle_ranges = [0 44;45 89;90 134;135 180];
% angle_ranges = [0 25.7;25.8 51.4;51.5 77.1;77.2 102.9;103 128.6;128.7 154.3;154.4 180];


reor_matrix = mark_reor_events_mod032317(tracks,reor_type,3);
outlawn_matrix = ~mark_inlawn_events_mod032317(tracks,3); %inverse of in lawn

frac_behav_angle = NaN(size(reor_matrix,1),size(angle_ranges,1));
for i = 1:size(angle_ranges,1)
    binary_angle_matrix = mark_angle_events_mod032317(tracks,reor_matrix,angle_ranges(i,1),angle_ranges(i,2),3);
    binary_angle_matrix_outside_lawn = binary_angle_matrix&outlawn_matrix; %only accept those distances where worms are outside the lawn
    % loop over the rows (each one corresponds to a track) and extract the
    % number of frames that satisfy angle condition --> denominator.
    % Then go to corresponding frames in reor_matrix and sum that -->
    % numerator. fraction is the "frequency" of particular reorientation
    % type in a given bearing angle.
    for j = 1:size(reor_matrix,1)
        %         if j == 14
        %             disp('debug');
        %         end
        total_angle_time = nansum(binary_angle_matrix_outside_lawn(j,:),2)/180;% per minute
        angle_inds = find(binary_angle_matrix_outside_lawn(j,:));
        total_behav_time = nansum(reor_matrix(j,angle_inds),2);
        if total_angle_time == 0
            frac_behav_angle(j,i) = NaN;
        else
            frac_behav_angle(j,i) = total_behav_time/total_angle_time;
        end
    end
end

avg_frac_behav_angle = nanmean(frac_behav_angle,1);
ste_frac_behav_angle = ste(frac_behav_angle,1);

%% 3. plot out reorientations versus bearing angle
figure(1);
errorbar(avg_frac_behav_angle,ste_frac_behav_angle);
set(gca,'XTick',1:size(angle_ranges,1));
set(gca,'XTickLabel',{'0-44','45-89','90-134','135-180'});
% set(gca,'XTickLabel',{'0-25.7','25.8-51.4','51.5-77.1','77.2-102.9','103-128.6','128.7-154.3','154.4-180'});
xlim([0.5 size(angle_ranges,1)+0.5])
% ylim([0.02 0.05]);
xlabel('odor bearing angle (deg)');
ylabel('frequency (event/min)');
% suplabel(filename,'t');
title([filename ' : ' reor_type ' reori vs. odor bearing angle_outsidelawn_1before'],'Interpreter','none');
% savefig([filename(1:end-4) '_' reor_type '_reoriVsodorangle_outsidelawn_1before.fig']);

%% 4. Extract worm frames at different distances to the odor source
% may want to exclude frames inside the lawn (don't do that yet though)

distance_ranges = [1 375;376 750;751 1125;1126 1500]; %pixels
distance_ranges_mm = distance_ranges.*0.03;

reor_matrix = mark_reor_events_mod032317(tracks,reor_type,3);
outlawn_matrix = ~mark_inlawn_events_mod032317(tracks,3); %inverse of in lawn

frac_behav_dist = NaN(size(reor_matrix,1),size(distance_ranges,1));
for i = 1:size(distance_ranges,1)
    binary_distance_matrix = mark_distance_events_mod032317(tracks,reor_matrix,distance_ranges(i,1),distance_ranges(i,2),3);
    binary_distance_matrix_outside_lawn = binary_distance_matrix&outlawn_matrix; %only accept those distances where worms are outside the lawn
    for j = 1:size(reor_matrix,1)
        total_distance_time = nansum(binary_distance_matrix_outside_lawn(j,:),2)/180; %per minute
        distance_inds = find(binary_distance_matrix_outside_lawn(j,:));
        total_behav_time = nansum(reor_matrix(j,distance_inds),2);
        if total_distance_time == 0
            frac_behav_dist(j,i) = NaN;
        else
%             if total_behav_time ~= 0
%                 disp('debug');
%             end
            frac_behav_dist(j,i) = total_behav_time/total_distance_time;
        end
    end
end

avg_frac_behav_dist = nanmean(frac_behav_dist,1);
ste_frac_behav_dist = ste(frac_behav_dist,1);

%% 5. plot out reorientations versus bearing distance
figure(2);
errorbar(avg_frac_behav_dist,ste_frac_behav_dist);
set(gca,'XTick',1:size(distance_ranges,1));
set(gca,'XTickLabel',{'0-11','12-22','23-34','35-45'});
% set(gca,'XTickLabel',{'0-25.7','25.8-51.4','51.5-77.1','77.2-102.9','103-128.6','128.7-154.3','154.4-180'});
xlim([0.5 size(distance_ranges,1)+0.5])
% ylim([0.02 0.05]);
xlabel('distance from odor (mm)');
ylabel('frequency (event/min)');
% suplabel(filename,'t');
title([filename ' : ' reor_type ' reori vs. odor distance_outsidelawn_1before'],'Interpreter','none');
% savefig([filename(1:end-4) '_' reor_type '_reoriVsodordist_outsidelawn_1before.fig']);
