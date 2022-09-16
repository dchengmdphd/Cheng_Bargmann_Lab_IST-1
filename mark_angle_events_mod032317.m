function [event_matrix] = mark_angle_events_mod032317(Tracks,reor_matrix, angle_min, angle_max,FrameRate)

%This type of analysis is to look at the behavior as consisting of reorientation events
%(like Bernoulli trials, 0-fail, 1-success)- in this case 0-fwd
%1-reorientation event

%Inputs
%Tracks- an array that contains the output of .linkedTracks, generated
%using join_multiple_sansfood_track_arrays
%reor_type- the reorientation type that you want to flag as an event.
%There are various options here (i)'all'- will take all reorientation as an event
%(ii) 'allnonUpsilon'- will take everything except Upsilon (anything but 3.1) (iii) enter a
%specific pirroutte string. See num_state_convert. Here you should enter the part of the reorientation that occurs
%first- for example if you're lookign for lRevOmega, you should
%enter 'lRevOmega' as opposed to 'OmegalRev'.
%Options for (iii) are:
%'lRevOmega'
%'sRevOmega'
%'lRevUpsilon'
%'sRevUpsilon'
%'pure_Upsilon'
%'pure_lRev'
%'pure_sRev'
%'pure_omega'

%Output
%event_matrix- a matrix (each column is a worm, each row is a frame) 
%that contains '1' where an event occurs for the first time and '0' for the
%rest. The 1 will only be at the frame where the event starts. frames
%where the track doesn't exist will be filled with NaNs.


%Converts the array into matrices where each column is an animal and each
%row is a frame
frames = track_field_to_matrix_mod032317(Tracks, 'frames');
Time = frames ./ FrameRate;
odor_angle_matrix = track_field_to_matrix_mod032317(Tracks, 'odor_angle');

reor_matrix = reor_matrix';
reor_start_inds = find(reor_matrix); %these are starts of behavior of interest
track_start_reor_inds = reor_start_inds(mod(reor_start_inds,size(odor_angle_matrix,1))==1);
non_track_start_reor_inds = reor_start_inds(mod(reor_start_inds,size(odor_angle_matrix,1))~=1);

odor_angle_matrix(track_start_reor_inds) = NaN;%just in case any reorientations are the beginning of a track, don't replace with numbers from a previous track
odor_angle_matrix(non_track_start_reor_inds) = odor_angle_matrix(non_track_start_reor_inds - 1);%look to the angle right before the reorientation


%Creates empty matrix with NaNs
% time_matrix = NaN(size(odor_angle_matrix)); 
event_matrix = zeros(size(odor_angle_matrix));
% event_matrix(isnan(odor_angle_matrix_precise)==1) = NaN;

% Will mark a 1 if the odor_angle is not NaN and falls within the range
% specified between angle_min and angle_max

event_matrix(odor_angle_matrix>=angle_min & odor_angle_matrix<=angle_max) = 1;


end





