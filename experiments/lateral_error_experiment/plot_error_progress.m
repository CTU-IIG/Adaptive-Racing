% This is to show the effect of ILC
% So I do this on bags purely changing e, not v
stump = "lateral_error_experiment";
errors = load(sprintf("%s_ferror.csv", stump));
len = ceil(max(errors(:, 2)));
fprintf("Estimated trajectory length: %d\n", len);

lap_times = load(sprintf("%s_lap_time.csv", stump));
fprintf("Detected laps: %d\n", length(lap_times));

% I split the data in errors according to the lap time
% I don't care about the actual lap time, but it is an indicator that
% the shifted L was uploaded to the follower
% NOTE: THIS IS ONLY VALID HERE AS IN OTHER EXPERIMENTS THE UPLOAD WAS
%       TIMED!
laps = {};
mses = {};

last_time = 0;

for i = 1:length(lap_times)
    fprintf("Lap %d, lap time: %f\n", i, lap_times(i, 2) / 1e9);
    lap_data = errors(((errors(:, 1) > last_time) & (errors(:, 1) < lap_times(i, 1))), 2:3);
    
    last_time = lap_times(i, 1);

    % Fix first column
    lap_data(:, 1) = mod(round(lap_data(:, 1)), len);
    
    % Sort data
    %lap_data = sortrows(lap_data, 1);
    % I need an order preserving sort
    temp_sort = nan(len, 2);
    temp_sort(:, 1) = 0:(len-1);
    last_index = lap_data(1, 1);
    for j = 1:length(lap_data)
        % This is also done in the ROS node
        % Basically throw away (or in the ROS node case interpolate)
        % everything in between the last and current value.
        if last_index > (len-20) && lap_data(j, 1) < (last_index-20)
            temp_sort((last_index+1):end, 2) = nan;
            temp_sort(1:lap_data(j, 1), 2) = nan;
        else
            temp_sort((last_index+1):(lap_data(j, 1)), 2) = nan;
        end
        temp_sort(lap_data(j, 1)+1, 2) = lap_data(j, 2);
        last_index = lap_data(j, 1)+1;
    end
    
    
    % Select unique
    % https://www.mathworks.com/matlabcentral/answers/486944-how-to-remove-duplicate-rows-in-a-matrix-by-considering-the-repeated-values-in-a-specific-column
    % first occurence
%     [~, uidx] = unique(lap_data(:, 1), 'stable');
%     lap_data = lap_data(uidx, :);
    lap_data = temp_sort(~isnan(temp_sort(:, 2)), :);
    
    
    % ~~ last occurence
%     flipped = flip(lap_data);
%     [~, uidx] = unique(flip(lap_data(:, 1)), 'stable');
%     lap_data = flip(flipped(uidx, :));
    
    laps{i} = lap_data;
    fprintf("mean squared error: %f\n", mse(lap_data(:, 2)));
    mses{i} = mse(lap_data(:, 2));

end


figure
hold on

for i = [1, 2, 3, 4, 5]
    plot(laps{i}(:, 1), laps{i}(:, 2));
    
    csvwrite(sprintf("%s_lap_%d.csv", stump, i), laps{i});
end

csvwrite(sprintf("%s_lap_mse.csv", stump), [mses{:}]');