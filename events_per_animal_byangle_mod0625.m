function vector_events_per_animal = events_per_animal_byangle_mod0625(Tracks,reor_type, t1,t2, FrameRate, anglemin, anglemax)


%This function will give you a vector with the number of reorientations for each animal in a specific time window.
%(each element in the vector will represent an animal and the value of each
%element is the number of reorientations for that animal)
%Will only give animals that don't have NaN's in the entire duration of the t1 to t2 time window.

%Inputs:
%Tracks- an array that contains the output of .linkedTracks, generated
%using join_multiple_sansfood_track_arrays
%reor_type- the reorientation type that you want to flag as an event (eg 'omega',
%'upsilon', etc. - 'all' will take any reorientation as an event). See
%num_state_convert
%t1,t2- beggining and end of time window of interest (in minutes)
%FrameRate- frame rate of acquisition

[event_matrix_reor,~] = mark_reor_events(Tracks,reor_type, FrameRate); %generates event_matrix,. See mark_reor_events
event_matrix_reor = event_matrix_reor';
if isequal(t1,0)==1
    start = 1; %sets start to 1 if t1=0
elseif t1 > 0
    start = t1 * FrameRate * 60; %converts t1 from minutes to frames
end

finish = t2 * FrameRate * 60; %converts t2 from minutes to frames
events = event_matrix_reor(start:finish,:);%selects the events in the region of interest
full = sum(isnan(events))'; % create 'full' which counts the number of NaNs

event_matrix_reor(isnan(event_matrix_reor)) = 0; %get rid of NaNs so we can do logical AND

event_matrix = event_matrix_reor & mark_angle_events_mod032317(Tracks,event_matrix_reor, anglemin,anglemax,3)';


events = event_matrix(start:finish,:);
vector_events_per_animal = sum (events)';
vector_events_per_animal = vector_events_per_animal(full==0); %only keeps animals that have no NaNs in that time window (where full==0)

