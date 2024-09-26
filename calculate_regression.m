clearvars
clc
% Read the data ready to be regression
data = readtable('PHA_regression.xlsx');


% Fitting the linear model
% with 3 coefficients: Interception Angle (Alpha), Semi-major axis (a),
% Relative Gain rate of BIP to COG (different_delta_v_t_relative)
mdl = fitlm(data, 'ResponseVar', 'Gain_percent', 'PredictorVars', {'Alpha','a','different_delta_v_t_relative'});

% Display the results
disp(mdl);

% To get the coefficients specifically
coefficients = mdl.Coefficients.Estimate;
disp(coefficients);

% Predicted value
y = zeros(length(data.a),1);

for i = 1 : length(data.a)
    y(i) = coefficients(1) + ...
          coefficients(2)*data.Alpha(i) + coefficients(3)*data.a(i) + coefficients(4)*data.different_delta_v_t_relative(i) ;
end



% Plot the graph of actual value and predicted value
x = 1 : length(data.Object);
figure(1)
scatter(x,y,100,'MarkerEdgeColor', 'r', 'LineWidth', 2.5); m1 = 'estimated';
hold on
scatter(x,data.Gain_percent,100,'MarkerEdgeColor', 'b', 'LineWidth', 2.5); m2 = 'actual';
legend(m1 ,m2)
xticks(x);  
xticklabels(data.Object);
ylabel('Gain(%)')
diff = data.Gain_percent - y;




