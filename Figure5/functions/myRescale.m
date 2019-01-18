function [ output_re ] = myRescale( stimuli_time, output, output_time )
%Background: 
%1, Camera is working at the limit speed, the read out time for each
%frame may vary a little bit from file to file.
%2, Input and Output are recorded from different files. 
%3, The initiate of the each trial is well aligned (time point of 0) 

%When the time of each frame for input and output are different, use Input
%file's time line, resclae the Output by linear interpolation.

% arrange all the image_times data in a single vector
max_density_time = [];
max_density_time = [max_density_time; stimuli_time; output_time];

%find out the min and max of each image-times vector
m_min = min(min(stimuli_time), min(output_time));
m_max = max(max(stimuli_time), max(output_time));

% find out all the unique elements in the vector
max_density_time = unique(max_density_time);
ind = max_density_time >= m_min & max_density_time<=m_max;
max_density_time =  max_density_time(ind);

% generate the interpolation accounding to the max_density_time.
output_temp = interp1(output_time, output, max_density_time, 'linear'); 

[C,ia,ib] = intersect(max_density_time, stimuli_time);
output_re = output_temp(ia);
end