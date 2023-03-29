clear;

%path to GUROBI LP solver
addpath 'C:\gurobi811\win64\matlab'
% addpath '/opt/gurobi810/linux64/matlab'


continents = {'Europe','Europe','North_America','North_America','Asia','Asia','Africa','Africa','Oceania','South_America','South_America','Europe','Europe','Europe','Europe','Europe'};
countries = {'United_Kingdom','Germany','USA','Canada','Japan','India','South Africa','Kenya','Australia','Brazil','Peru','France','Sweden','Spain','Denmark','Italy'};

time_test_begin = datetime('05-01-2020','InputFormat','MM-dd-yyyy');
time_test_end = datetime('11-01-2021','InputFormat','MM-dd-yyyy');
time_points = time_test_begin:caldays(15):time_test_end;

mut_alpha = [501 570 614 681 716 982 1118];
mut_beta = [80 215 417 484 501 614 701];
mut_gamma = [18 20 26 138 190 417 484 501 614 655 1027 1176];
mut_delta = [19 142 158 452 478 614 681 950];
mut_omicron = [67 95 145 212 339 371 373 375 417 440 446 477 478 484 493 496 498 501 505 547 614 655 679 681 764 796 856 954 969 981];
mut_lambda = [75 76 452 490 614 859];
mut_mu = [95 145 346 484 501 614 681 950];
mut_theta = [265 484 501 614 681 1092 1101 1176];
mut_eta = [52 67 484 614 677 888];
mut_kappa = [142 154 452 484 614 681 1071];
mut_zeta = [484 614 1176];
mut_epsilon = [13 152 452 614];
mut_iota = [5 95 614];

nMut = 1273;
voc = false(10,nMut);
voc(1,mut_alpha) = 1;
voc(2,mut_beta) = 1;
voc(3,mut_gamma) = 1;
voc(4,mut_delta) = 1;
voc(5,mut_omicron) = 1;
voc(6,mut_lambda) = 1;
voc(7,mut_mu) = 1;
voc(8,mut_theta) = 1;
voc(9,mut_eta) = 1;
voc(10,mut_kappa) = 1;

% Set the country to be analyzed
cnt = 7;
continent = continents{cnt};
country = countries{cnt};

%% 
outdirRes = ['HELEN data1' filesep 'p-values' filesep country];
p_voc_time = zeros(size(voc,1),length(time_points));
for t = 1:length(time_points)
    for s = 1:size(voc,1)
        outfile = [outdirRes filesep 'res_s' int2str(s) '_t' int2str(t) '.mat'];
        load(outfile)
        p_voc_time(s,t) = p_voc;
    end
end

