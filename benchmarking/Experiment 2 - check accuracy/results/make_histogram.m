clear
clc

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["cac_mine", "cac_esri"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
cacscores = readtable("cacs.csv", opts);


%% Clear temporary variables
clear opts

close all
cac = table2array(cacscores);
bins = 0:0.05:1;

histogram(cac(:, 1), 'BinEdges', bins, 'EdgeColor', 'b', 'FaceAlpha', 0.5);
hold on
histogram(cac(:, 2), 'BinEdges', bins, 'EdgeColor', 'r', 'FaceAlpha', 0.5);
legend("delineator.py", "ESRI ArcGIS Online", 'Location', 'northwest');
xlabel("Coefficient of Areal Correspondence");
ylabel("Count of Watersheds");
%title(["USGS gage watersheds", "compared to delineation by USGS"], FontWeight="normal")
fontsize(gcf, 14, "points");

mysavefig(gcf, 'us_watersheds_compare')
