clear; clc;

rootDir = fileparts(mfilename('fullpath'));
cd(rootDir);
addpath(genpath(rootDir));

% name, data file, h, BGDR k, CGDR k
configs = {
    % 'MSRC_v1',     'MSRC_v1_data.mat',     7,  36, 28;
    'Dermatology', 'Dermatology_data.mat', 8,  38, 32;
    % '100leaves',   '100leaves_data.mat',   10, 30, 16;
    % 'mnist4',      'mnist4_data.mat',      10, 14, 10;
    % 'Digit4k',     'Digit4k_data.mat',     10, 40, 40;
    % 'Hdigit',      'Hdigit_data.mat',      12, 40, 32;
    % 'ALOI',        'ALOI_data.mat',        13, 8,  8;
    % 'MNIST',       'MNIST6w_data.mat',     8,  32, 38;
};

initLabel = 'N2HI';
isNormal = 1;
maxIter = 30;

fprintf('Dataset     Method    ACC    NMI    Purity  Precison Recall F-score ARI Time(s)\n');
fprintf('-----------------------------------------------------\n');

for i = 1:size(configs, 1)
    name = configs{i, 1};
    file = configs{i, 2};
    h = configs{i, 3};
    kBGDR = configs{i, 4};
    kCGDR = configs{i, 5};

    data = load(fullfile(rootDir, 'Data', file));
    X = data.X;
    label = data.label(:);

    rng(1);
    [resBGDR,~,~,~,tBGDR] = MvC_BGDR(X,label,kBGDR,h,initLabel,isNormal,maxIter);
    fprintf('%-12s BGDR    %.4f  %.4f  %.4f  %.4f %.4f %.4f %.4f %.2f\n', name, resBGDR(1), resBGDR(2),...
        resBGDR(3), resBGDR(4), resBGDR(5), resBGDR(6), resBGDR(7), tBGDR);

    rng(1);
    [resCGDR,~,~,~,tCGDR] = MvC_CGDR(X,label,kCGDR,h,initLabel,isNormal,maxIter);
    fprintf('%-12s CGDR    %.4f  %.4f  %.4f  %.4f %.4f %.4f %.4f %.2f\n', name, resCGDR(1), resCGDR(2),...
        resCGDR(3), resCGDR(4), resCGDR(5), resCGDR(6), resCGDR(7), tCGDR);
end
