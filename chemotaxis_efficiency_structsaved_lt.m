% % % % CHEMOTAXIS EFFICIENCY ANALYSIS V.6 03/23/16
% % This function analyzes chemotaxis efficiency and saves a struct
% % containing extracted features.
% % It now reads in linkedTracks.mat files, asks user to input chemotaxis
% % region as a circle as a way of avoiding the chemotaxis tracker.
function chemotaxis_efficiency_structsaved_lt( filename, min_track_length, last_edited_frame_global, figure_flag )
% TO DO: put in varargin statement for min_track_length and good_radius
% (defaults are 1000, and 400, respectively)
close all;

v = VideoReader(filename);
% [filename, pathname, ~] = uigetfile({'*.avi'}); %%COMMENTED OUT!!!!
% v = VideoReader([pathname filename]);
% circle the lawn or chemotaxis destination region
figure, imshow(imadjust(rgb2gray(read(v,last_edited_frame_global)))); hold on;
h2 = imellipse;
position = wait(h2);
hmask2 = (h2.createMask());
api = iptgetapi(h2);
vert = api.getVertices();
tgt_x = vert(:,1);
tgt_y = vert(:,2);
dest_pos = getPosition(h2);
% close all;
bottom_left = [dest_pos(1) dest_pos(2)];
top_right = [dest_pos(1)+dest_pos(3) dest_pos(2)+dest_pos(4)];
dest_center = [(bottom_left(1)+top_right(1))/2 (bottom_left(2)+top_right(2))/2]; % midpoint of opposite diagonal vertices of ellipse bounding rectangle
close all;

% linkedTracks_file = [pathname filename(1:end-4) '.linkedTracks.mat'];%%COMMENTED OUT!!!!
linkedTracks_file = [filename(1:end-4) '.linkedTracks.mat'];
tracks = load_Tracks(linkedTracks_file);
total_numtracks = length(tracks);
if total_numtracks ~= 0
%     last_edited_frame_global = 8100;
else
    error('THERE ARE NO TRACKS IN THIS LINKEDTRACKS FILE!');
end

% DOUBLE CHECK THAT EVERYTHING IS A DOUBLE!!!!!
%% grab tracks which terminate in the odor target region
% this definition of good tracks works because we are using azide to
% paralyze worms once they reach the odor.


good_tracks = [];
for i = 1:length(tracks)
    if length(tracks(i).Frames)>min_track_length
        good_tracks = [good_tracks; tracks(i)];
    end
end
good_numtracks = length(good_tracks);
good_bool = 0;
if good_numtracks == 0
    disp('There are no good tracks. Will not output ...good.mat file');
else
    good_bool = 1;
     %change to double precision for future calculations
    good_state = double(track_field_to_matrix(good_tracks,'State')); good_state = pad_matrix(good_state,size(good_state,1),last_edited_frame_global,NaN);
    good_frames = double(track_field_to_matrix(good_tracks,'Frames')); good_frames = pad_matrix(good_frames,size(good_frames,1),last_edited_frame_global,NaN);
    good_smoothx = double(track_field_to_matrix(good_tracks,'SmoothX')); good_smoothx = pad_matrix(good_smoothx,size(good_smoothx,1),last_edited_frame_global,NaN);
    good_smoothy = double(track_field_to_matrix(good_tracks,'SmoothY')); good_smoothy = pad_matrix(good_smoothy,size(good_smoothy,1),last_edited_frame_global,NaN);
    good_head_angle = double(track_field_to_matrix(good_tracks,'head_angle')); good_head_angle = pad_matrix(good_head_angle,size(good_head_angle,1),last_edited_frame_global,NaN);
    good_eccentricity = double(track_field_to_matrix(good_tracks,'Eccentricity')); good_eccentricity = pad_matrix(good_eccentricity,size(good_eccentricity,1),last_edited_frame_global,NaN);
end
%% PLOT good TRACKS BEHAVIORAL SEGMENTATION
if good_bool
    e_a = sqrt(good_numtracks);
    if rem(e_a,1)==0 %figure out the dimensions of the multi-panel plot
        e_a = fix(e_a);
        e_b = e_a;
    else
        e_a = fix(e_a)+1;
        tmp = rem(good_numtracks/e_a,1);
        if tmp==0
            e_b = good_numtracks/e_a;
        else
            e_b = ceil(good_numtracks/e_a);
        end
    end
    
    if figure_flag
        [ax1, h1] = make_panel_figure_xy(e_a,e_b,2); suplabel('good TRACKS: ALL BEHAVIORS','t');
        
        for i = 1:good_numtracks
            firstframe = good_tracks(i).Frames(1);
            lastframe = good_tracks(i).Frames(end);
            frames = firstframe:lastframe;
            curr_x = good_smoothx(i,frames)';
            curr_y = good_smoothy(i,frames)';
            
            state_vect = good_state(i,:);
            
            %extract behavioral state code indices
            pure_fwd = find(abs(state_vect-1.0)<.001);
            pure_pause = find(abs(state_vect-1.1)<.001);
            lRevOmega = find(abs(state_vect-4.7)<.001);
            OmegalRev = find(abs(state_vect-7.4)<.001);
            sRevOmega = find(abs(state_vect-5.7)<.001);
            OmegasRev = find(abs(state_vect-7.5)<.001);
            lRevUpsilon = find(abs(state_vect-4.3)<.001);
            UpsilonlRev = find(abs(state_vect-3.4)<.001);
            sRevUpsilon = find(abs(state_vect-5.3)<.001);
            UpsilonsRev = find(abs(state_vect-3.5)<.001);
            pure_Upsilon = find(abs(state_vect-3.1)<.001);
            pure_lRev = find(abs(state_vect-4.1)<.001);
            pure_sRev = find(abs(state_vect-5.1)<.001);
            pure_omega = find(abs(state_vect-7.1)<.001);
            ring = find(abs(state_vect-99)<.001);
            missing = find(abs(state_vect-100)<.001);
            num_behav = 16;
            % plot color according to behavior
            axes(ax1(i)); hold on;
            cc = hsv(num_behav);
            colormap(cc);
            
            scatter(curr_x(pure_fwd-firstframe+1),curr_y(pure_fwd-firstframe+1),5,cc(1,:),'filled','o');
            scatter(curr_x(pure_pause-firstframe+1),curr_y(pure_pause-firstframe+1),5,cc(2,:),'filled','o');
            scatter(curr_x(lRevOmega-firstframe+1),curr_y(lRevOmega-firstframe+1),5,cc(3,:),'filled','o');
            scatter(curr_x(OmegalRev-firstframe+1),curr_y(OmegalRev-firstframe+1),5,cc(4,:),'filled','o');
            scatter(curr_x(sRevOmega-firstframe+1),curr_y(sRevOmega-firstframe+1),5,cc(5,:),'filled','o');
            scatter(curr_x(OmegasRev-firstframe+1),curr_y(OmegasRev-firstframe+1),5,cc(6,:),'filled','o');
            scatter(curr_x(lRevUpsilon-firstframe+1),curr_y(lRevUpsilon-firstframe+1),5,cc(7,:),'filled','o');
            scatter(curr_x(UpsilonlRev-firstframe+1),curr_y(UpsilonlRev-firstframe+1),5,cc(8,:),'filled','o');
            scatter(curr_x(sRevUpsilon-firstframe+1),curr_y(sRevUpsilon-firstframe+1),5,cc(9,:),'filled','o');
            scatter(curr_x(UpsilonsRev-firstframe+1),curr_y(UpsilonsRev-firstframe+1),5,cc(10,:),'filled','o');
            scatter(curr_x(pure_Upsilon-firstframe+1),curr_y(pure_Upsilon-firstframe+1),5,cc(11,:),'filled','o');
            scatter(curr_x(pure_lRev-firstframe+1),curr_y(pure_lRev-firstframe+1),5,cc(12,:),'filled','o');
            scatter(curr_x(pure_sRev-firstframe+1),curr_y(pure_sRev-firstframe+1),5,cc(13,:),'filled','o');
            scatter(curr_x(pure_omega-firstframe+1),curr_y(pure_omega-firstframe+1),5,cc(14,:),'filled','o');
            scatter(curr_x(ring-firstframe+1),curr_y(ring-firstframe+1),5,cc(15,:),'filled','o');
            scatter(curr_x(missing-firstframe+1),curr_y(missing-firstframe+1),5,cc(16,:),'filled','o');
            if i ==1
                legend('pure fwd','pure pause','lRevOmega','OmegalRev','sRevOmega','OmegasRev','lRevUpsilon','UpsilonlRev','sRevUpsilon','UpsilonsRev','pure Upsilon','pure lRev','pure sRev','pure omega','ring','missing');
            end
            text(double(curr_x(1)),double(curr_y(1)),'START','FontSize',8);
            text(double(curr_x(end)),double(curr_y(end)),'END','FontSize',8);
            scatter(dest_center(1),dest_center(2),50,'g','filled');
            axis equal;
        end
    end
end



%% Extract parameters needed for defining the end of an good track (also load them for other tracks)
% calculate the instantaneous acceleration of the centroid path
% define the time window to examine track (before paralysis)
% using speed and distance from target (hopefully this will be enough)
if good_bool
    good_speed = double(track_field_to_matrix(good_tracks,'Speed'));
    good_odor_distance = NaN(size(good_smoothx));
    for i = 1:good_numtracks
        firstframe = good_tracks(i).Frames(1);
        lastframe = good_tracks(i).Frames(end);
        frames = firstframe:lastframe;
        G = [good_smoothx(i,frames)' good_smoothy(i,frames)'];
        G2 = [dest_center(1).*ones(size(G,1),1) dest_center(2).*ones(size(G,1),1)];
        V = G - G2;
        D = sqrt(sum(abs(V).^2,2));
        good_odor_distance(i,frames) = D;
    end
    good_angular_speed_navin = double(track_field_to_matrix(good_tracks,'AngSpeed'));
end

%% EXTRACT BODY CURVATURE, ANGULAR SPEED, TURNING EVENTS FOR good TRACKS
if good_bool
    if exist('good_body_contour_curve_all','var') %this is here for technical reasons of repeatedly calling this code block leading to problems
        clear('good_body_contour_curve_all');
    end
    %these vectors are the same size as the ones above but their last entry should be NaN
    good_path_angles = NaN(good_numtracks,size(good_smoothx,2));
    good_path_angles_differentiable = NaN(good_numtracks,size(good_smoothx,2));
    
    good_angular_velocity = NaN(good_numtracks,size(good_smoothx,2)); %first angular velocity will always be NaN
    good_turning_events_matrix = zeros(good_numtracks,size(good_smoothx,2));
    good_turning_events = cell(good_numtracks,1);
    
    good_odor_angles = NaN(good_numtracks,size(good_smoothx,2));
    
    good_path_curvature = NaN(good_numtracks,size(good_smoothx,2));
    
    good_bigcurveMat = NaN(good_numtracks*last_edited_frame_global,55);%last two numbers (52+3=55) are worm index and frame index (starting at 1 is beginning of track)
    good_body_contour_curve_all(good_numtracks) = struct();
    % odor vector is in previous version of this file
    for i = 1:good_numtracks
        curr_track = good_tracks(i);
        frames = good_frames(i,:);
        frames = frames(~isnan(frames));
        curr_x = good_smoothx(i,frames)';
        curr_y = good_smoothy(i,frames)';
        x_3rdpt = curr_x+1; %this is needed for angle calculations using law of cosines
        next_x = curr_x(2:end);
        next_y = curr_y(2:end);
        
        ind = 1;
        prev_theta_sign = 0;
        curr_theta_sign = 0;
        prev_theta = 0;
        curr_diffable_theta = 0;
        prev_diffable_theta = 0;
        if exist('body_contour_curve_frame','var')
            clear('body_contour_curve_frame');
        end
        body_contour_curve_frame(length(frames)) = struct();
        for j = 1:length(next_x)
            %this function uses the law of cosines to get the angle between three
            %points A, B, C (angle B is returned)
            
            % odor_theta: these angles will all be positive and range from 0-180
%             if j == 98
%                 disp('debug')
%             end
            [~, odor_theta_deg] = angle_between_three_points_raddeg([dest_center(1) dest_center(2)],[curr_x(j) curr_y(j)],[next_x(j) next_y(j)]);
            good_odor_angles(i,frames(ind)) = odor_theta_deg;
            % theta: these angles are direction of movement in arena coordinates
            [~, theta_deg] = angle_between_three_points_raddeg([x_3rdpt(j) curr_y(j)],[curr_x(j) curr_y(j)],[next_x(j) next_y(j)]);
            if (next_y(j) - curr_y(j)) < 0
                theta = -1*theta_deg; %make angles negative only if y difference is less than 0
            else
                theta = theta_deg;
            end
            curr_theta_sign = sign(theta);
            good_path_angles(i,frames(ind)) = theta;
            if j == 1
                curr_diffable_theta = theta;
                good_path_angles_differentiable(i,frames(ind)) = curr_diffable_theta;
                
            else %if j>=2
                ang_diff = min(360-(abs(theta - prev_theta)), abs(theta - prev_theta)); % choose the lower angular difference
                option_sub = prev_diffable_theta - ang_diff; option_add = prev_diffable_theta + ang_diff;
                if sign(option_sub) == sign(prev_diffable_theta) % if subtraction gives the same sign, use it
                    curr_diffable_theta = option_sub;
                elseif sign(option_add) == sign(prev_diffable_theta) % if addition gives the same sign, use it
                    curr_diffable_theta = option_add;
                elseif sign(option_sub == 0) % we know we can't match the last sign. if subtraction gives zero, use it.
                    curr_diffable_theta = option_sub;
                elseif sign(option_add == 0) % we know we can't match the last sign. if addition gives zero, use it.
                    curr_diffable_theta = option_add;
                elseif 0 == sign(prev_diffable_theta) % if previous sign was 0, add angle difference by convention
                    curr_diffable_theta = option_add;
                end
                good_path_angles_differentiable(i,frames(ind)) = curr_diffable_theta;
                
            end
            prev_theta = theta;
            prev_theta_sign = curr_theta_sign;
            prev_diffable_theta = curr_diffable_theta;
            
            % for extracting body curvature
            j_track_ind = frames(j)-frames(1)+1;
            x_pos = single(curr_track.body_contour(j_track_ind).x);
            y_pos = single(curr_track.body_contour(j_track_ind).y);
            if length(x_pos)==52 && length(y_pos)==52
                image = curr_track.Image{j_track_ind};
                curvature_vs_body = getcurvature(x_pos', y_pos');
                curvature_vs_body = curvature_vs_body - mean(curvature_vs_body); %reorients curvature vector to worm-centered coordinates
                body_contour_curve_frame(j).x = x_pos;
                body_contour_curve_frame(j).y = y_pos;
                body_contour_curve_frame(j).curvature = curvature_vs_body;
                body_contour_curve_frame(j).image = image;
                good_bigcurveMat(ind,:) = [curvature_vs_body' good_tracks(i).State(j_track_ind) i j];
            end
            
            ind = ind+1; % this ind is a running tally of all worms in all frames
        end
        
% 03/03/16 not sure why i had it ending in end-1 before. seems to work fine
% ending at end of frames.
%         abs_ang_diff = abs([NaN diff(good_path_angles_differentiable(i,frames(1:end-1)))])';
%         good_angular_velocity(i,frames(1:end-1)) = abs([NaN diff(good_path_angles_differentiable(i,frames(1:end-1)))]');
        
        abs_ang_diff = abs([NaN diff(good_path_angles_differentiable(i,frames(1:end)))])';
        good_angular_velocity(i,frames) = abs([NaN diff(good_path_angles_differentiable(i,frames))]');
        
        % for the Hongkyun Kim Levy Flight Turning Event definition
        te_thresh = 10;
        good_turning_events_matrix(i,frames(abs_ang_diff >= te_thresh)') = frames(abs_ang_diff >= te_thresh)';
        good_turning_events_matrix(i,frames(1)) = frames(1);%frame numbers -- frame 1 is always the first turning event.
        good_turning_events{i} = find(good_turning_events_matrix);
        
        % for getting path curvature
        good_path_curvature(i,frames) = getcurvature(curr_x,curr_y);
%         highcurv_inds = find(curvature>0.2);%these inds are where path curvature is highest -- often only for a single or few frames
        
        good_body_contour_curve_all(i).body_contour_curve = body_contour_curve_frame;
    end
    good_bigcurveMat(isnan(good_bigcurveMat(:,1)),:) = [];
end

%% NOW SAVE ALL good TRACKS INFO INTO ONE .MAT FILE AND ALL OTHER TRACKS INTO ANOTHER
% make sure that you include the good_radius and the coordinates of the
% target point in the good_tracks.mat file

to_save_pfx = linkedTracks_file(1:end-4);
if good_bool
    track_features = struct();
    track_features(good_numtracks).frames = []; %initialize the struct
    for i = 1:good_numtracks
        track_features(i).frames = good_frames(i,:);
        track_features(i).smoothx = good_smoothx(i,:);
        track_features(i).smoothy = good_smoothy(i,:);
        track_features(i).speed = good_speed(i,:);
        track_features(i).state = good_state(i,:);
        track_features(i).head_angle = good_head_angle(i,:);
        track_features(i).eccentricity = good_eccentricity(i,:);
        track_features(i).body_contour_curve_all = good_body_contour_curve_all(i); %each is a struct
        track_features(i).path_angle = good_path_angles(i,:);
        track_features(i).path_angle_differentiable = good_path_angles_differentiable(i,:);
        track_features(i).path_curvature = good_path_curvature(i,:);
        track_features(i).turning_events = good_turning_events{i}; % each is a cell
        track_features(i).turning_eventsMat = good_turning_events_matrix(i,:);
        track_features(i).target_center = [dest_center(1)*ones(good_numtracks,1) dest_center(2)*ones(good_numtracks,1)]; %same every time
        track_features(i).target_verts = vert; %same every time
        track_features(i).odor_distance = good_odor_distance(i,:);
        track_features(i).odor_angle = good_odor_angles(i,:);

    end
    %what to do about bigcurveMat? leave it out for now 03/23/2017
    
%     track_features.angular_velocity = good_angular_velocity;
% %     track_features.frames = good_frames;
% %     track_features.odor_angles = good_odor_angles;
% %     track_features.odor_distance = good_odor_distance;
% %     track_features.path_angles = good_path_angles;
% %     track_features.path_angles_differentiable = good_path_angles_differentiable;
% %     track_features.path_curvature = good_path_curvature;
% %     track_features.smoothx = good_smoothx;
% %     track_features.smoothy = good_smoothy;
% %     track_features.speed = good_speed;
% %     track_features.state = good_state;
% %     track_features.turning_eventsMat = good_turning_events_matrix;
% %     track_features.turning_events = good_turning_events;
% %     track_features.target_verts = vert;
% %     track_features.target_center = [dest_center(1)*ones(good_numtracks,1) dest_center(2)*ones(good_numtracks,1)];
% %     track_features.body_contour_curve_all = good_body_contour_curve_all;
% ?    track_features.bigcurveMat = good_bigcurveMat;
% %     track_features.head_angle = good_head_angle;
% %     track_features.eccentricity = good_eccentricity;
    save([to_save_pfx '.GOOD.mat'],'track_features');
end

end