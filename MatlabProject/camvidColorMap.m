function cmap = camvidColorMap()
% Define the colormap used by CamVid dataset.

cmap = [
    255 255 255%white
    255 0 0   % Red
    0 255 0       % Green
    0 0 0   % Black
    
    ];

% Normalize between [0 1].
cmap = cmap ./ 255;
end