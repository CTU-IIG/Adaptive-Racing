% Plot /path/visualization gradient
file = "environment_estimation_experiment.visualization.%d.csv";
files = 1:10:897;


range = 1:length(files);


figure

hold on

axis equal

grid on
ax = gca;
ax.GridLineStyle = '--';

color = rgb2hsv(ax.ColorOrder(1, :));
colors = rgb2hsv((jet(length(range))));
% alphas = logspace(-1, 0, length(files));
alphas = linspace(0, 1, length(files));


for i = range
    filename = sprintf(file, files(i));
    
    data = load(filename);
    
%     current_color = color;
%     current_color(2) = (current_color(2)) * (i / length(range)) ;

    current_color = colors(i, :);

    data = [ data; data(1, :) ];
    
    plot(data(:, 1), data(:, 2), 'Color', [hsv2rgb(current_color), alphas(i)]);
%     hsv2rgb(current_color)
   
end



title("/path/visualization gradient");
xlabel("x [m]");
ylabel("y [m]");


return



%% Test txy plot; see plot_error_progress.m
% txy_file = "environment_estimation_experiment_odom.txy.csv";
% lap_time_file = "environment_estimation_experiment_lap_time.csv";
% save = 0;


txy = load(txy_file);
lap_times = importdata(lap_time_file);

fprintf("Detected laps: %d\n", length(lap_times.data));

lap_times = lap_times.data;

num_laps = length(lap_times);

last_time = 0;

figure

hold on

axis equal

grid on
ax = gca;
ax.GridLineStyle = '--';

color = rgb2hsv(ax.ColorOrder(1, :));
colors = rgb2hsv((jet(num_laps)));
% alphas = logspace(-1, 0, length(files));
alphas = linspace(0, 1, num_laps);

last_point = [];

for i = 1:num_laps
    fprintf("Lap %d, lap time: %f\n", i, lap_times(i, 2) / 1e9);
    lap_data = txy(((txy(:, 1) > last_time) & (txy(:, 1) < lap_times(i, 1))), 2:3);
    
    if length(last_point) > 0
        lap_data = [last_point; lap_data];
    end
    
    last_point = lap_data(end, :);
    
    last_time = lap_times(i, 1);
    
    current_color = colors(i, :);
    
    plot(lap_data(:, 1), lap_data(:, 2), 'Color', [hsv2rgb(current_color), alphas(i)]);
    
    if save
        csvwrite(sprintf("environment_estimation_experiment_odom.lap.%d.csv", i), lap_data);
    end
end
