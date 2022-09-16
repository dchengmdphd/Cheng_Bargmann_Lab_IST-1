function [event_matrix] = mark_inlawn_events_mod032317(Tracks,FrameRate)

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


smoothx_matrix = track_field_to_matrix_mod032317(Tracks, 'smoothx');
smoothy_matrix = track_field_to_matrix_mod032317(Tracks, 'smoothy');

target_verts = Tracks(1).target_verts;
tgt_x = target_verts(:,1);
tgt_y = target_verts(:,2);

event_matrix = inpolygon(smoothx_matrix,smoothy_matrix,tgt_x,tgt_y);
end





