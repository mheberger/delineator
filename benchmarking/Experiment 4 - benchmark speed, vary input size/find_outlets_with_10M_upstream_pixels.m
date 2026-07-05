%% Find outlets where there are 10 million upstream cells,
% for a benchmarking experiment #4. 

clear
clc
close all

% The number of pixels we're looking for. 10 million
PIX = 10^7;

basins = [11];

n_basins = length(basins);

accum_folder = "C:/Data/GIS/MERITHydro/accum_basins";
mask_folder = "C:/Data/GIS/MERITBasins/masks";

new_accum_folder = "C:/Data/GIS/MERITHydro/accum_basins_rev";


% Let T be my outlet table, with cols basin, lat, lng (coordinates for an outlet)
T = nan(n_basins, 3);

% Let U be the number of upstream cells
% (Store separately, otherwise table gets output as integers not floats)
U = zeros(n_basins, 1);


%% Iterate over the raster files
i = 1;
for basin = 11
    accum_file = sprintf("%s/accum%g.tif", accum_folder, basin);
    [A, R] = readgeoraster(accum_file);

    mask_file = sprintf("%s/mask%g.tif", mask_folder, basin);
    [M, ~] = readgeoraster(mask_file);
    
    % Apply the mask!
    M = logical(M);
    M = ~M;
    A(M) = 0;

    % Write the new, masked geotiff file
    new_accum_file = sprintf("%s/accum%g.tif", new_accum_folder, basin);
    geotiffwrite(new_accum_file, A, R, 'TiffTags', struct('Compression', 'DEFLATE'));

    % Extract some metadata from the raster info
    top = R.LatitudeLimits(2);
    width = R.CellExtentInLatitude;
    half_pix = width / 2;
    left = R.LongitudeLimits(1);

    %A(A==247) = 0;
    % Find the index of a cell where the accum. is 10 million pixels
    % Find elements above the threshold
    indices_above_threshold = find(A > PIX);

    % Get the values above threshold
    values_above_threshold = A(indices_above_threshold);

    [~, idx] = min(abs(values_above_threshold - PIX));

    % Get the index of the element closest to the threshold
    ind = indices_above_threshold(idx);

    % Now figure out the lat, lng coordinates of that cell
    if ~isempty(ind)
        u = A(ind);
        disp(u)
        U(i) = u;
        [y, x] = ind2sub(size(A), ind);
        y = y - 1;
        x = x - 1;
        outlat = top - y * width;
        outlng = left + x * width;
        T(i, :)  = [basin; outlat; outlng];
    else
        T(i,1) = basin;
    end

    disp(T(i,:));
  
    i = i + 1;

end

cd(fileparts(matlab.desktop.editor.getActiveFilename));

writematrix(T, 'expt4_outlets.csv');

fanfare
