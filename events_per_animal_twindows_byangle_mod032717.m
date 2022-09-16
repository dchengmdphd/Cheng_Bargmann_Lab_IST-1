function [events_twindows,time] = events_per_animal_twindows_byangle_mod032717(Tracks,time_window,reor_type,FrameRate, anglemin, anglemax)

%This function will output events_twindows- a matrix where each column is a time window
%and each row is the number of reorientations of each animal in that time
%window
%Will also output time which is just a vector with the corresponding time
%of each time window

%Inputs:
%Tracks- file array
%Time_window- bin window in minutes. Should be a multiple of the total
%time
%reor_type- type of reorientation. To look at all write 'all' (needs to be a
%string)
%FrameRate
%max_time- how long the movie lasts in minutes

%First determining the length of the longest video (in minutes) in the Tracks array
Frame_matrix=track_field_to_matrix_mod032317(Tracks,'frames')'; %matrix of all frames
max_time = length(Frame_matrix)./60./FrameRate; %length of this matrix is equal to the longest video

%Will change the max time of the video to ensure that it is a multiple of
%the the time-window. This way I can have uniform time_windows.

extra_time = rem(max_time,time_window);
max_time = max_time - extra_time;

%Pre-allocating
events_twindows = NaN(size(Frame_matrix,1),max_time/time_window);

%Will first use events_per_animal function to calculate the number of events that each animal does in each time bin.

%angle ranges to test
% angle_ranges = [0 44;45 89;90 134;135 180];

time_bin_num = 1;
for m = 0:time_window:max_time-time_window
    vector_events_per_animal = events_per_animal_byangle_mod032717(Tracks,reor_type, m,m+time_window, FrameRate, anglemin,anglemax);
    events_twindows(1:length(vector_events_per_animal),time_bin_num) = vector_events_per_animal;
    time_bin_num = time_bin_num+1;
end

%Will now trim down the events_twindows
max_events = max(sum(isnan(events_twindows)==0));
events_twindows(max_events+1:end,:) = [];


%Creating correspondign time vector for plotting purposes
time = linspace(time_window,max_time,max_time/time_window); %creates a time vector where each value corresponds to a bin

end
