% Load and preprocess the image
img = imread('chx rep 2.png');  % Read the image file 
gray_img = rgb2gray(img);  % Convert the image to grayscale
super_gray_img = imgaussfilt(gray_img);  % Apply slight Gaussian smoothing
% Display the image for ROI selection
figure; imshow(super_gray_img);
title('Select ROIs');
% ROI Selection
replicate_names = {'Parental', 'Clone 9', 'Clone 16', 'Clone 177', 'Clone 320'};
rois_per_replicate = 3;  % Number of ROIs per replicate
total_rois = numel(replicate_names) * rois_per_replicate;  % Total ROIs
roi_masks = false(size(gray_img, 1), size(gray_img, 2), total_rois);  % Initialize masks
for i = 1:numel(replicate_names)
   for j = 1:rois_per_replicate
       roi_index = (i-1) * rois_per_replicate + j;
       fprintf('Select ROI %d for %s\n', j, replicate_names{i});  % Display instruction
       roi = drawrectangle;  % Allow the user to draw a rectangle
      
       roi_masks(:,:,roi_index) = createMask(roi);  % Get the binary mask for the current ROI
   end
end
close;  % Close the figure after ROI selection

% Highlighted image for visualizing used pixels
highlighted_image = img;

%% Calculate average darkness and band width for each ROI
mean_darkness = zeros(numel(replicate_names), rois_per_replicate); 
band_area = zeros(numel(replicate_names), rois_per_replicate);    
for i = 1:numel(replicate_names)
   for j = 1:rois_per_replicate
       roi_index = (i-1) * rois_per_replicate + j;
       roi_img = super_gray_img .* uint8(roi_masks(:,:,roi_index)); 
       % Detect edges within the ROI (Canny edge detection)
       edges = edge(roi_img, 'canny', 0.14);
       % Morphological closing to connect edges
       se = strel('disk',3,4);
       edges = imclose(edges, se);
       % Remove ROI rectangle edges from the edge map
       roi_edges = bwperim(roi_masks(:,:,roi_index));
       se = strel('disk', 1); 
       expanded_edges = imdilate(roi_edges, se); 
       edges(expanded_edges) = 0; 
       % Fill the interior bounded by edges
       filled_region = imfill(edges, 'holes');
       % Select the largest connected component
       CC = bwconncomp(filled_region);
       region_sizes = cellfun(@numel, CC.PixelIdxList);
       [~, largest_region_idx] = max(region_sizes); 
       filled_region = false(size(filled_region));
       filled_region(CC.PixelIdxList{largest_region_idx}) = true;
      % Calculate band area
       band_area(i, j) = calculateBandArea(filled_region);
       % Highlight the used pixels on the original image
       for k = 1:3 
           channel = highlighted_image(:,:,k);
           if k == 2 
               channel(filled_region) = 255; 
           else 
               channel(filled_region) = 0;
           end
           highlighted_image(:,:,k) = channel;
       end
       % Calculate the mean darkness of pixels within the filled region
       used_pixels = roi_img(filled_region); 
       % Invert the grayscale values
inverted_pixels = 255 - used_pixels; 
mean_darkness(i,j) = sum(inverted_pixels) / sum(filled_region(:));
   end
end
%%
% Display the highlighted image
figure; imshow(highlighted_image);
title('Highlighted Pixels Used for Analysis');
%%
% Display mean darkness values
disp('Mean Darkness for each ROI:');
for i = 1:numel(replicate_names)
   fprintf('%s: ', replicate_names{i});
   disp(mean_darkness(i,:));
end
%%
% Normalize darkness values to the first band of each replicate AND by band width
concentrations = zeros(numel(replicate_names), rois_per_replicate);
for i = 1:numel(replicate_names)
   concentrations(i,:) = (mean_darkness(i,:) / mean_darkness(i,1)) .* (band_area(i,:) / band_area(i,1));
end
% Fit exponential decay curve for each replicate
time_points = [0 9 18; 0 9 18;0 9 18; 0 9 18; 0 9 18];  % Time points for each replicate (adjust if needed)
half_life = zeros(1, numel(replicate_names));
r_squared_values = zeros(1, numel(replicate_names));
figure; hold on;  % Plot all decay curves on the same figure
for i = 1:numel(replicate_names)
   % Fit an exponential decay curve to the data
   fo = fitoptions('Method', 'NonlinearLeastSquares', ...
               'StartPoint', 0, ...
               'Lower', 0, ...       
               'Upper', Inf);         
   ft = fittype('1*exp(-b*x)', 'options', fo);
    [curve, gof] = fit(time_points(i,:)', concentrations(i,:)', ft);
   % Store R-squared value
   r_squared_values(i) = gof.rsquare;
   % Display goodness-of-fit statistics
   fprintf('Replicate %s:\n', replicate_names{i});
   fprintf('  R-squared: %.4f\n', gof.rsquare);
   fprintf('  SSE: %.4f\n', gof.sse);
   fprintf('  RMSE: %.4f\n', gof.rmse);
   % Calculate half-life for the current replicate
   half_life(i) = log(2) / curve.b;
   % Plot the concentration decay curve
   p = plot(time_points(i,:), curve(time_points(i,:)), 'DisplayName', replicate_names{i});
 end
xlabel('Time (hours)');
ylabel('Normalized Concentration');
title('Concentration Decay after Knockdown');
legend show;
hold off; 
% Analyze half-life across replicates
disp('Half-life for each Replicate:');
for i = 1:numel(replicate_names)
   fprintf('%s: %f hours\n', replicate_names{i}, half_life(i));
end
% Calculate and display average and standard deviation of half-life
avg_half_life = mean(half_life);
std_half_life = std(half_life);
disp(['Average Half-life: ', num2str(avg_half_life), ' +/- ', num2str(std_half_life), ' hours']);

%% Save the plot to a file
saveas(gcf, 'concentration_decay_plot_1.png');

%% Function to calculate band area
function area = calculateBandArea(filled_region)
area = sum(filled_region(:)); % Simply sum the pixels in the filled region
end
