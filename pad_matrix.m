function output_matrix = pad_matrix( input_matrix, target_height, target_width, stuffitwith )
%PAD_MATRIX This function is essentially a wrapper for the native MATLAB
%function padarray, which just adds "stuffitwith" around an existing matrix
%to fill it out to a desired dimension (must be greater than current
%dimension).
in_height = size(input_matrix,1);
in_width = size(input_matrix,2);
if target_width < in_width || target_height < in_height
    error('Target dimensions must be greater than or equal to input dimensions!');
end
height_to_add = target_height - in_height;
width_to_add = target_width - in_width;
output_matrix = padarray(input_matrix,[height_to_add width_to_add],stuffitwith,'post');
end

