%{
Function help or script explanation.
No need of putting all the '%'.
%}

%% Initialization
close all; clear; clc            

%% Adding to path
% this command add to the path all the folders and subfolders on you're current path
addpath(genpath(fileparts(pwd)))                

%% Figure Setup    
load('MagellanoColorMap.mat'); % make sure you added the colormap file in the current folder or local subfolders
DefaultOrderColor = get(0, 'DefaultAxesColorOrder');
NewOrderColor = [0.9490    0.4745    0.3137
                 0.1020    0.6667    0.74120
                 155/255   155/255   155/255
                 DefaultOrderColor];  
             
set(0, 'DefaultFigureColormap', MagellanoColorMap);
set(0, 'DefaultAxesColorOrder', NewOrderColor);
set(0, 'DefaultLineLineWidth', 2)
set(0, 'DefaultLineMarkerSize', 10)
set(0, 'DefaultFigureUnits', 'normalized');
set(0, 'DefaultFigurePosition', [0 0 1 1]);
set(0, 'DefaultTextFontSize', 18);
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultAxesXGrid', 'on')
set(0, 'DefaultAxesYGrid', 'on')
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultTextInterpreter', 'Latex');

