function [event_matrix,time_matrix] = mark_reor_events_mod032317(Tracks,reor_type,FrameRate)

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
state_matrix_precise = track_field_to_matrix_mod032317(Tracks, 'state');
state_matrix = floor(state_matrix_precise); %this eliminates decimal points in the state matrix so that fwd and pause are 1


%Creates empty matrix with NaNs
time_matrix = NaN(size(state_matrix)); 
binary_reor = zeros(size(state_matrix));
event_matrix = zeros(size(state_matrix));
event_matrix(isnan(state_matrix_precise)==1) = NaN;

%Will first mark reorientations with a 1 (and the rest of the track is 0).
%What gets marked with a 1 depends on what the reopr_type is.
if strcmpi(reor_type,'all') == 1
   binary_reor(state_matrix ~= 1) = 1; %make anything that isn't forward or pause 1
   binary_reor(state_matrix == 99) = 0; %make the ring zero
   binary_reor(isnan(state_matrix_precise)==1) = NaN; %placing NaNs where the track doesnt exist
   
elseif strcmpi(reor_type,'allnonUpsilon')== 1
   binary_reor(state_matrix ~= 1) = 1; %make anything that isn't forward or pause 1
   binary_reor(abs(state_matrix_precise - 3.1) <0.1) = 0; %Makes upsilons 0
   binary_reor(abs(state_matrix_precise - 7.1) <0.1) = 0; %Makes pure Omega 0
   binary_reor(state_matrix == 99) = 0; %make the ring zero
   binary_reor(isnan(state_matrix_precise)==1) = NaN; %placing NaNs where the track doesnt exist
   
else pir_value = num_state_convert(reor_type);
   binary_reor(abs(state_matrix_precise - pir_value) <0.01) = 1;
   binary_reor(state_matrix == 99) = 0; %make the ring zero
   binary_reor(isnan(state_matrix_precise)==1) = NaN; %placing NaNs where the track doesnt exist
end



%Will now mark with a '1', only the first frame the reorientation
%starts (event_matrix) and will also record that time (time_matrix)

for k = 1: size(state_matrix,2)
    reor_indices = find(binary_reor(:,k) == 1);
    if isempty(reor_indices) == 0
    event_matrix(reor_indices(1),k) = 1;
    time_matrix(reor_indices(1),k) = Time(reor_indices(1),k);
    end
    for m = 2:length(reor_indices)
        if (reor_indices(m) - reor_indices(m-1)) > 1
            event_matrix(reor_indices(m),k) = 1;
            time_matrix(reor_indices(m),k) = Time(reor_indices(m),k);
        end
    end
end





