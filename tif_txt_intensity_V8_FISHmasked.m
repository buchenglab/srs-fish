%% initialization and file name
%
% clear all;
% close all;
% clc;

%% use as a function

function tif_txt_intensity_V8_FISHmasked(filename_sample_0, grayscale, num, filepath)
close all;
%% specified colorbar color, no tick and tick label
%% FISH threshold & FISH single cell contour
%% set grayscale to 1 if need auto %CD range
%% set num to 2 if need multichannel avg


CH_threshold = 0.02;
CD_tag = 0;
ratio_tag = 1;


if nargin < 2
    grayscale = 0;
    num = 2; % autoaverage
    filepath = pwd;
elseif nargin < 3
    num = 2;
    filepath = pwd;
elseif nargin < 4
    filepath = pwd;
end

filename_sample_CD = [filepath,'/*', filename_sample_0, '*_CD*.txt_aligned.tiff'];
filename_sample_CH = [filepath,'/*', filename_sample_0, '*_CH.txt_aligned.tiff'];
filename_sample_off = [filepath,'/*', filename_sample_0, '*_off.txt_aligned.tiff'];
filename_sample_Cy3 = [filepath,'/*', filename_sample_0, '*_Cy3.txt_masked.tiff'];
filename_sample_Cy5 = [filepath,'/*', filename_sample_0, '*_Cy5.txt_masked.tiff'];

%% resize due to alignments
% size(filename_sample_CD)

% filename_sample = [filepath,'\', filename_sample_0, '.txt']; % Windows
% using backward slash, mac using forward slash
l = 3;  % 3rd largest point to avoid extreme data

color = 1;  % jet plot or gray scale
range = 1;  % fix color bar or adjust according to max & min of data
%%%%%%%%% change range here
if CD_tag == 1
%     min_data = 0.1;
%     max_data = 1.5;
elseif ratio_tag == 1
    min_data = 0;
    max_data = 0.25;
else
%     min_data = 0;
%     max_data = 30;
end

if grayscale == 1
    color = 0;
    range = 0;
end
%% import

% sample_CD = read_hSRS(filename_sample_CD,CHANNEL_NUM);
% sample_CH = read_hSRS(filename_sample_CH,CHANNEL_NUM);
% sample_off = read+hSRS(filename_sample_off,CHANNEL_NUM);

t1 = dir(filename_sample_CD);
t1 = Tiff(t1.name);
sample_CD = read(t1);
t2 = dir(filename_sample_CH);
t2 = Tiff(t2.name);
sample_CH = read(t2);
t3 = dir(filename_sample_off);
t3 = Tiff(t3.name);
sample_off = read(t3);
t4 = dir(filename_sample_Cy3);
t4 = Tiff(t4.name);
sample_Cy3 = read(t4);
t5 = dir(filename_sample_Cy5);
t5 = Tiff(t5.name);
sample_Cy5 = read(t5);

CompFact = 1; % off-resonance subtraction compensation factor
sample_off = CompFact*sample_off;
[row_s, column_s, lambda_s] = size(sample_CD);

% if num == 2
%     sample_CD(:,:,num) = sum(sample_CD,3)/CHANNEL_NUM;
%     sample_CH(:,:,num) = sum(sample_CH,3)/CHANNEL_NUM;    
% end

%% ratio and plotting mask creation: thresholding
CD_ratio = (sample_CD-sample_off)./(sample_CD+sample_CH-2*sample_off+1e-4);
mask1 = sample_CH > CH_threshold;
% mask2 = sample_Cy5(:,:,num) > Cy5_threshold;
mask3 = sample_Cy3|sample_Cy5;
% imshow(CD_ratio);
% CD_ratio_masked = CD_ratio.*(mask1|mask2);
CD_ratio_masked = CD_ratio.*mask1.*mask3;

%% plot

stepsize = 0.520e0*0.006/0.02; % um

X = [1:row_s]*stepsize;
Y = [1:column_s]*stepsize;

% [X,Y] = meshgrid(Y,X);
% mesh(X,Y,sample_data(:,:,num));

set(0, 'DefaultTextFontSize', 20);
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesLIneWidth', 2);
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultAxesFontName', 'Arial');

f = figure;
RatioMap = imagesc(CD_ratio_masked);
axis image;
h = gca;
f.Units = 'pixels';


%% color and contrast
cyan = [0, 0.75, 0.75];
lightBlue = [0, 0.5, 1];
lightOrange = [0.9372, 0.7480, 0.2549];
deepGreen = [0, 0.5686, 0];
deepBlue = [0 0.1569 0.5608];
purple = [1, 0, 1];
darkRed = [0.6350, 0.0780, 0.1840];
red = [1 0 0];
yellow = [1 1 0];
deepOrange = [0.5373 0.4196 0];

CD_top = [1 0.9137 0.8353];
CD = lightOrange;

CD_top = [0.9804 0.8078 0.2196];
CD = deepOrange;

CH_top = [0.8353 1 0.8353];
CH = deepGreen;

ratio_top = CD_top;
ratio_color = CH;


if color == 1
    %     colormap(gca,colMapGen(Cy3_570_top, [0 0 0], 2000, 'midCol',Cy3_570));
    %     colormap(gca,colMapGen(Cy5_670_top, [0 0 0], 2000, 'midCol',Cy5_670));
    if CD_tag == 1
        colormap(h,colMapGen(CD_top, [0 0 0], 2000, 'midCol',CD));
    else
        if ratio_tag ~= 1
            colormap(h,colMapGen(CH_top, [0 0 0], 2000, 'midCol',CH));
        else
%             colormap(gca,colMapGen([0,0.4,0],deepBlue, 2000, 'midCol', cyan));
            colormap(magma(2000));
        end
    end
end



if range == 0
    A1 = maxk(sample_CH(:,:,num), l, 2);
    A2 = A1(:,end);
    max_data = max(A2);
    % max_data = maxk(maxk(sample_data(:,:,num), k),k);
    min_data = min(min(sample_CH(:,:,num)));
end
caxis([min_data, max_data]);


%% hide edges
% h = gca;
% % ylabel('\mum');
h.XAxis.Visible = 'off';
h.YAxis.Visible = 'off';
h.Position=[0 0 1 1];   % default for gca is normalized position

%%% bare figure version
f.Position = [0 0 row_s, column_s];

% saveas(RatioMap, [filepath,'/RatioPlot/', filename_sample_0,'Ratio.tif'],'tiffn');
% exportgraphics(gca, [filepath,'/RatioPlot/', filename_sample_0,'Ratio.tif'], 'Resolution', 200);


% save like this may lead to size/pixel number different in the saved
% format
% saveas(RatioMap, [filepath,'\RatioPlot\',filename_sample_0,'Ratio.tif']);
saveas(RatioMap, [filepath,'/',filename_sample_0,'Ratio.tif']);
% %%% colorbar without box version
% colorbar('Box', 'off', 'XTickLabel', {},'XTick',[], 'TickLength', 0, 'XColor', [1 1 1]);
% f.Position = [0 0 row_s+50 column_s+50]; % for colorbar


pause(1.5);
%% contour
f_contour_Cy3 = figure;
interp_ratio = 1;
Cy3_contour_ROI = interp2(sample_Cy3,interp_ratio)>0.5;
Cy3_contour_line = edge(Cy3_contour_ROI);
Cy5_contour_ROI = interp2(sample_Cy5,interp_ratio)>0.5;
Cy5_contour_line = edge(Cy5_contour_ROI);
% visualization
imagesc(Cy3_contour_line, 'AlphaData',0.7);
colormap(gray);
h_contour = gca;
axis square;
f_contour_Cy3.Units = 'pixels';
h_contour.XAxis.Visible = 'off';
h_contour.YAxis.Visible = 'off';
h_contour.Position=[0 0 1 1];   % default for gca is normalized position
[row_cy3, column_cy3] = size(Cy3_contour_ROI);
f_contour_Cy3.Position = [0 0 row_cy3, column_cy3];
pause(1.5);

f_contour_Cy5 = figure;;
imagesc(Cy5_contour_line, 'AlphaData',0.7);
colormap(gray);
h_contour = gca;
axis square;
f_contour_Cy5.Units = 'pixels';
h_contour.XAxis.Visible = 'off';
h_contour.YAxis.Visible = 'off';
h_contour.Position=[0 0 1 1];   % default for gca is normalized position
[row_cy5, column_cy5] = size(Cy5_contour_ROI);
f_contour_Cy5.Position = [0 0 row_cy5, column_cy5];



% if you want to save as specific pixel size figure
% imwrite(imresize(im2uint16(CH_contour_line),[200 200]), [filepath,'/RatioPlot/',filename_sample_0,'C.tif']);

% did not imresize and did not specify, the tif is in compressed mode
imwrite(im2uint16(Cy3_contour_line), [filepath,'/',filename_sample_0,'Cy3.tif']);
imwrite(im2uint16(Cy5_contour_line), [filepath,'/',filename_sample_0,'Cy5.tif']);

%% save data
% filename_output=strcat(filename_sample_0, '_pattern_reomoved.tif' );
% pattern_removed_sample_data = double(pattern_removed_sample_data);
% for K=1:lambda_s
%     imwrite(pattern_removed_sample_data(:, :, K), filename_output, 'WriteMode', 'append', 'Compression', 'none');
% end


%% functions
% read_hSRS('FA3 1_FISH.txt.tif')
% when use as a function, this should be commented
function hSRS_data = read_hSRS(filename,CHANNEL_NUM)
        [fPath, fName, fExt] = fileparts(filename);
        
        switch lower(fExt)
            case '.txt'
                %read from txt
                hSRS_data = dlmread(filename);
                %                 CHANNEL_NUM = 10;
                hSRS_data=permute(reshape(hSRS_data,[size(hSRS_data,1),size(hSRS_data,2)/CHANNEL_NUM,CHANNEL_NUM]),[1,2,3]);
                %     hSRS_data=hSRS_data(:,1+5:end-5,:);
                %             hSRS_data_resize=hSRS_data(:,1+5:end-5,:);
            case '.tif'
                %read from tif
                info = imfinfo(filename);
                num_images = numel(info);
                for k = 1:num_images
                    hSRS_data(:,:,k) = imread(filename, k, 'Info', info);
                    hSRS_data=single(hSRS_data);
                end
            otherwise  % Under all circumstances SWITCH gets an OTHERWISE!
                error('Unexpected file extension: %s', fExt);
        end
    end

function map = magma(N)
% Perceptually uniform sequential colormap from MatPlotLib.
%
% Copyright (c) 2017-2020 Stephen Cobeldick
%
%%% Syntax:
%  map = magma
%  map = magma(N)
%
% Colormap designed by Nathaniel J. Smith and Stefan van der Walt.
%
% For MatPlotLib 2.0 improved colormaps were created in the perceptually
% uniform colorspace CAM02-UCS. The new colormaps are introduced here:
% <http://matplotlib.org/2.0.0rc2/users/dflt_style_changes.html>
% The RGB data is from here: <https://bids.github.io/colormap/>
%
% Note VIRIDIS replaces the awful JET as the MatPlotLib default colormap.
%
%% Examples %%
%
%%% Plot the scheme's RGB values:
% rgbplot(magma(256))
%
%%% New colors for the COLORMAP example:
% load spine
% image(X)
% colormap(magma)
%
%%% New colors for the SURF example:
% [X,Y,Z] = peaks(30);
% surfc(X,Y,Z)
% colormap(magma)
% axis([-3,3,-3,3,-10,5])
%
%% Input and Output Arguments %%
%
%%% Inputs (*=default):
% N = NumericScalar, N>=0, an integer to define the colormap length.
%   = *[], use the length of the current figure's colormap (see COLORMAP).
%
%%% Outputs:
% map = NumericMatrix, size Nx3, a colormap of RGB values between 0 and 1.
%
% See also CIVIDIS INFERNO PLASMA VIRIDIS TWILIGHT TAB10 SET LINES COLORMAP PARULA

if nargin<1 || isnumeric(N)&&isequal(N,[])
	N = size(get(gcf,'colormap'),1);
else
	assert(isscalar(N)&&isfinite(N)&&isreal(N),...
		'SC:magma:N:NotRealFiniteScalarNumeric',...
		'First argument <N> must be a real finite numeric scalar.')
end
%
raw = [0.001462,0.000466,0.013866; 0.002258,0.001295,0.018331; 0.003279,0.002305,0.023708; 0.004512,0.003490,0.029965; 0.005950,0.004843,0.037130; 0.007588,0.006356,0.044973; 0.009426,0.008022,0.052844; 0.011465,0.009828,0.060750; 0.013708,0.011771,0.068667; 0.016156,0.013840,0.076603; 0.018815,0.016026,0.084584; 0.021692,0.018320,0.092610; 0.024792,0.020715,0.100676; 0.028123,0.023201,0.108787; 0.031696,0.025765,0.116965; 0.035520,0.028397,0.125209; 0.039608,0.031090,0.133515; 0.043830,0.033830,0.141886; 0.048062,0.036607,0.150327; 0.052320,0.039407,0.158841; 0.056615,0.042160,0.167446; 0.060949,0.044794,0.176129; 0.065330,0.047318,0.184892; 0.069764,0.049726,0.193735; 0.074257,0.052017,0.202660; 0.078815,0.054184,0.211667; 0.083446,0.056225,0.220755; 0.088155,0.058133,0.229922; 0.092949,0.059904,0.239164; 0.097833,0.061531,0.248477; 0.102815,0.063010,0.257854; 0.107899,0.064335,0.267289; 0.113094,0.065492,0.276784; 0.118405,0.066479,0.286321; 0.123833,0.067295,0.295879; 0.129380,0.067935,0.305443; 0.135053,0.068391,0.315000; 0.140858,0.068654,0.324538; 0.146785,0.068738,0.334011; 0.152839,0.068637,0.343404; 0.159018,0.068354,0.352688; 0.165308,0.067911,0.361816; 0.171713,0.067305,0.370771; 0.178212,0.066576,0.379497; 0.184801,0.065732,0.387973; 0.191460,0.064818,0.396152; 0.198177,0.063862,0.404009; 0.204935,0.062907,0.411514; 0.211718,0.061992,0.418647; 0.218512,0.061158,0.425392; 0.225302,0.060445,0.431742; 0.232077,0.059889,0.437695; 0.238826,0.059517,0.443256; 0.245543,0.059352,0.448436; 0.252220,0.059415,0.453248; 0.258857,0.059706,0.457710; 0.265447,0.060237,0.461840; 0.271994,0.060994,0.465660; 0.278493,0.061978,0.469190; 0.284951,0.063168,0.472451; 0.291366,0.064553,0.475462; 0.297740,0.066117,0.478243; 0.304081,0.067835,0.480812; 0.310382,0.069702,0.483186; 0.316654,0.071690,0.485380; 0.322899,0.073782,0.487408; 0.329114,0.075972,0.489287; 0.335308,0.078236,0.491024; 0.341482,0.080564,0.492631; 0.347636,0.082946,0.494121; 0.353773,0.085373,0.495501; 0.359898,0.087831,0.496778; 0.366012,0.090314,0.497960; 0.372116,0.092816,0.499053; 0.378211,0.095332,0.500067; 0.384299,0.097855,0.501002; 0.390384,0.100379,0.501864; 0.396467,0.102902,0.502658; 0.402548,0.105420,0.503386; 0.408629,0.107930,0.504052; 0.414709,0.110431,0.504662; 0.420791,0.112920,0.505215; 0.426877,0.115395,0.505714; 0.432967,0.117855,0.506160; 0.439062,0.120298,0.506555; 0.445163,0.122724,0.506901; 0.451271,0.125132,0.507198; 0.457386,0.127522,0.507448; 0.463508,0.129893,0.507652; 0.469640,0.132245,0.507809; 0.475780,0.134577,0.507921; 0.481929,0.136891,0.507989; 0.488088,0.139186,0.508011; 0.494258,0.141462,0.507988; 0.500438,0.143719,0.507920; 0.506629,0.145958,0.507806; 0.512831,0.148179,0.507648; 0.519045,0.150383,0.507443; 0.525270,0.152569,0.507192; 0.531507,0.154739,0.506895; 0.537755,0.156894,0.506551; 0.544015,0.159033,0.506159; 0.550287,0.161158,0.505719; 0.556571,0.163269,0.505230; 0.562866,0.165368,0.504692; 0.569172,0.167454,0.504105; 0.575490,0.169530,0.503466; 0.581819,0.171596,0.502777; 0.588158,0.173652,0.502035; 0.594508,0.175701,0.501241; 0.600868,0.177743,0.500394; 0.607238,0.179779,0.499492; 0.613617,0.181811,0.498536; 0.620005,0.183840,0.497524; 0.626401,0.185867,0.496456; 0.632805,0.187893,0.495332; 0.639216,0.189921,0.494150; 0.645633,0.191952,0.492910; 0.652056,0.193986,0.491611; 0.658483,0.196027,0.490253; 0.664915,0.198075,0.488836; 0.671349,0.200133,0.487358; 0.677786,0.202203,0.485819; 0.684224,0.204286,0.484219; 0.690661,0.206384,0.482558; 0.697098,0.208501,0.480835; 0.703532,0.210638,0.479049; 0.709962,0.212797,0.477201; 0.716387,0.214982,0.475290; 0.722805,0.217194,0.473316; 0.729216,0.219437,0.471279; 0.735616,0.221713,0.469180; 0.742004,0.224025,0.467018; 0.748378,0.226377,0.464794; 0.754737,0.228772,0.462509; 0.761077,0.231214,0.460162; 0.767398,0.233705,0.457755; 0.773695,0.236249,0.455289; 0.779968,0.238851,0.452765; 0.786212,0.241514,0.450184; 0.792427,0.244242,0.447543; 0.798608,0.247040,0.444848; 0.804752,0.249911,0.442102; 0.810855,0.252861,0.439305; 0.816914,0.255895,0.436461; 0.822926,0.259016,0.433573; 0.828886,0.262229,0.430644; 0.834791,0.265540,0.427671; 0.840636,0.268953,0.424666; 0.846416,0.272473,0.421631; 0.852126,0.276106,0.418573; 0.857763,0.279857,0.415496; 0.863320,0.283729,0.412403; 0.868793,0.287728,0.409303; 0.874176,0.291859,0.406205; 0.879464,0.296125,0.403118; 0.884651,0.300530,0.400047; 0.889731,0.305079,0.397002; 0.894700,0.309773,0.393995; 0.899552,0.314616,0.391037; 0.904281,0.319610,0.388137; 0.908884,0.324755,0.385308; 0.913354,0.330052,0.382563; 0.917689,0.335500,0.379915; 0.921884,0.341098,0.377376; 0.925937,0.346844,0.374959; 0.929845,0.352734,0.372677; 0.933606,0.358764,0.370541; 0.937221,0.364929,0.368567; 0.940687,0.371224,0.366762; 0.944006,0.377643,0.365136; 0.947180,0.384178,0.363701; 0.950210,0.390820,0.362468; 0.953099,0.397563,0.361438; 0.955849,0.404400,0.360619; 0.958464,0.411324,0.360014; 0.960949,0.418323,0.359630; 0.963310,0.425390,0.359469; 0.965549,0.432519,0.359529; 0.967671,0.439703,0.359810; 0.969680,0.446936,0.360311; 0.971582,0.454210,0.361030; 0.973381,0.461520,0.361965; 0.975082,0.468861,0.363111; 0.976690,0.476226,0.364466; 0.978210,0.483612,0.366025; 0.979645,0.491014,0.367783; 0.981000,0.498428,0.369734; 0.982279,0.505851,0.371874; 0.983485,0.513280,0.374198; 0.984622,0.520713,0.376698; 0.985693,0.528148,0.379371; 0.986700,0.535582,0.382210; 0.987646,0.543015,0.385210; 0.988533,0.550446,0.388365; 0.989363,0.557873,0.391671; 0.990138,0.565296,0.395122; 0.990871,0.572706,0.398714; 0.991558,0.580107,0.402441; 0.992196,0.587502,0.406299; 0.992785,0.594891,0.410283; 0.993326,0.602275,0.414390; 0.993834,0.609644,0.418613; 0.994309,0.616999,0.422950; 0.994738,0.624350,0.427397; 0.995122,0.631696,0.431951; 0.995480,0.639027,0.436607; 0.995810,0.646344,0.441361; 0.996096,0.653659,0.446213; 0.996341,0.660969,0.451160; 0.996580,0.668256,0.456192; 0.996775,0.675541,0.461314; 0.996925,0.682828,0.466526; 0.997077,0.690088,0.471811; 0.997186,0.697349,0.477182; 0.997254,0.704611,0.482635; 0.997325,0.711848,0.488154; 0.997351,0.719089,0.493755; 0.997351,0.726324,0.499428; 0.997341,0.733545,0.505167; 0.997285,0.740772,0.510983; 0.997228,0.747981,0.516859; 0.997138,0.755190,0.522806; 0.997019,0.762398,0.528821; 0.996898,0.769591,0.534892; 0.996727,0.776795,0.541039; 0.996571,0.783977,0.547233; 0.996369,0.791167,0.553499; 0.996162,0.798348,0.559820; 0.995932,0.805527,0.566202; 0.995680,0.812706,0.572645; 0.995424,0.819875,0.579140; 0.995131,0.827052,0.585701; 0.994851,0.834213,0.592307; 0.994524,0.841387,0.598983; 0.994222,0.848540,0.605696; 0.993866,0.855711,0.612482; 0.993545,0.862859,0.619299; 0.993170,0.870024,0.626189; 0.992831,0.877168,0.633109; 0.992440,0.884330,0.640099; 0.992089,0.891470,0.647116; 0.991688,0.898627,0.654202; 0.991332,0.905763,0.661309; 0.990930,0.912915,0.668481; 0.990570,0.920049,0.675675; 0.990175,0.927196,0.682926; 0.989815,0.934329,0.690198; 0.989434,0.941470,0.697519; 0.989077,0.948604,0.704863; 0.988717,0.955742,0.712242; 0.988367,0.962878,0.719649; 0.988033,0.970012,0.727077; 0.987691,0.977154,0.734536; 0.987387,0.984288,0.742002; 0.987053,0.991438,0.749504];
%
num_ = size(raw,1);
% With small extrapolation when N>num:
vec = linspace(0,num_+1,N+2);
map = interp1(1:num_,raw,vec(2:N+1),'linear','extrap');
% Interpolation only for all values of N:
%map = interp1(1:num,raw,linspace(1,num,N),'spline')
% Range limits:
map = max(0,min(1,map));
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%magma
end