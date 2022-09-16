function x = track_field_to_matrix_mod032317(Tracks, field, nan_replacement)
% x = track_field_to_matrix(Tracks, field, nan_replacement)
% returns a matrix length(Tracks)-by-max(length(Tracks.field)) filled w/
% the value of the field
% NaN for missing values replaced by nan_replacement

x=[];
if(~isfield(Tracks,field))
    return;
end

x = zeros(length(Tracks), max_struct_array(Tracks,'frames'),'double') + NaN;

for(i=1:length(Tracks))
%     temp = Tracks(i).frames;
%     temp = temp(~isnan(temp));
%     first_frame = temp(1);
%     last_frame = temp(end);
    x(i,:) = Tracks(i).(field);
end

% fill the end of each row with the final value instead of NaN
if(nargin==3)
    if(issstr(nan_replacement))
        if(strcmpi(nan_replacement,'end'))
            for(i=1:length(Tracks))
                x(i,Tracks(i).frames(end):size(x,2)) = Tracks(i).(field)(end);
            end
        end
        return;
    end
    
    if(~issstr(nan_replacement))
        x = matrix_replace(x,'==',NaN,nan_replacement);
    end
end

return;
end
