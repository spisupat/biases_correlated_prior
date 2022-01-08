close all
% clear

% Try loading in summary statistics from current colder
try
    load('/data/Voc_DataMatrix.mat','ex_voc','ex_env','df')
    load('/data/RaterScore.mat'); % load inclusion/exclusion decisions 
    load('/data/SyntheticSounds.mat')
catch
    warning('Vocalization data or Rater Decisions not found in current folder');
end

% clearvars -except datmatrix R df Fs
lowerFreq=14000;  % Lower bound of peak frequency
upperFreq=34000; % Upper bound of peak frequency


% Index frequency bounds
centerFreq = df.PeakFreq;
fbounds = centerFreq>lowerFreq & centerFreq<upperFreq & RaterScores; % Index for vocalizations that are within bounds
dfv = df(fbounds,:); % summary statistics table with non 22kHz vocalizations removed

dur = 1000.*(dfv.thoff-dfv.thon);   % Calculation of duration
onset_slope = (((dfv.Ahon-dfv.Aon)./dfv.AUC)./abs(dfv.thon-dfv.ton)); % Calculation of onset slope
offset_slope = (((dfv.Ahoff-dfv.Aoff)./dfv.AUC)./abs(dfv.toff-dfv.thoff)); % Calculation of offset slope

idx = ~isinf(onset_slope) & ~isinf(offset_slope);
onset_slope = onset_slope(idx);
offset_slope = offset_slope(idx);
dur = dur(idx);

AUC_slow_condition = trapz([1:Synthetic.Fs]./Synthetic.Fs,Synthetic.Slow)/2;
AUC_fast_condition = trapz([1:Synthetic.Fs]./Synthetic.Fs,Synthetic.Fast)/2;
slowDerivative =(diff(Synthetic.Slow)./AUC_slow_condition).*Synthetic.Fs; % Derivative of synthetic sounds (Fast slope condition)
fastDerivative = (diff(Synthetic.Fast)./AUC_fast_condition).*Synthetic.Fs; % Derivative of synthetic sounds (Slow slope condition)
 
x = [1:Synthetic.Fs]./Synthetic.Fs;

% Onset/offset slopes of synthetic sounds (both slope conditions: fast and slow)
slow = max(slowDerivative); % Slope rate at inflection point of synthetic sounds
fast = max(fastDerivative);    % Slope rate at inflection point of synthetic sounds


given = [slow, fast]; % Conditional slope values from synthetic sounds used from the discrimination task
color = {'r','b'};  % Color scheme (slope conditions: red=slow; blue=fast)
slopes = [onset_slope, offset_slope]; % Matrix for all three slope variables

model_type = {'onset_slope', 'offset_slope', 'joint_slope'}'; % Slope variables used for different prior types
parameter_label = {'mu','sigma2','mu|slow','sigma2|slow','mu|fast','sigma2|fast'}; % Labels of parameters for bivariate and conditional gaussians for empirical prior (mu = means, sigma2 = variance)
prior = array2table(NaN(length(model_type),length(parameter_label)),'variablenames',parameter_label,'rownames',model_type); % construct table

prior.('mu') = {mean([dur,onset_slope],'omitnan'),mean([dur,offset_slope],'omitnan'),mean([dur,onset_slope,offset_slope],'omitnan')}'; % Means of multivariate gaussians
prior.('sigma2') = {cov([dur,onset_slope],'omitrows'),cov([dur,offset_slope],'omitrows'), cov([dur,onset_slope,offset_slope],'omitrows')}'; % Covariances of multivariate gaussians
for mdl_idx =1:length(model_type)
    [prior.('mu|slow')(mdl_idx), prior.('sigma2|slow')(mdl_idx)] = condGauss(prior.('mu'){mdl_idx},prior.('sigma2'){mdl_idx},slow);  % Conditional Gaussian parameters for slow slope condition 
    [prior.('mu|fast')(mdl_idx), prior.('sigma2|fast')(mdl_idx)] = condGauss(prior.('mu'){mdl_idx},prior.('sigma2'){mdl_idx},fast); %Conditional Gaussian parameters for fast slope condition 
end
 

%%%%%%%%%%%%%%%%%%%%%%%
%  Large figure panel %
%%%%%%%%%%%%%%%%%%%%%%%

figure('units','normalized','outerposition',[0 0 1 1])
CondLab = {'Onset','Offset'};
x = [0:1500]';

% Overlay regression line w/ confidence intervals on scatter plot

    x_domain =  [0:5:1500];
    y_domain =  [0:2.5:225];
for i = 1:2
    subplot(4,4,8+i)
    [yhat,ci] = predict(fitlm(dur,slopes(:,i)),x,'Alpha',0.05);
    h = scatter(dur,slopes(:,i),40,'.','filled');
    hold on
    h.MarkerFaceColor = [.6,.6,.6];
    h.MarkerEdgeColor = [.6,.6,.6];
    plot(x,yhat,'k-','linewidth',2)
    plot(x,ci,'k:','linewidth',1)
    ylim([0,225])
    xlim([0,1500])
    set(gca,'TickDir','out');
    colorbar
    set(colorbar,'visible','off')  % colorbar adjustment to align subplots
    xlabel('Plateau Duration (ms)')
    ylabel('Slope (A/s)')
    title(CondLab{i});
    axis square

    % Plot bivariate probability density functions as a smoothed surface
    ax=subplot(4,4,12+i);
    mu = prior.('mu'){i}; sigma = prior.('sigma2'){i};
    [X, Y] = meshgrid(x_domain,y_domain);
    variable_space = [X(:), Y(:)];
    Px = mvnpdf(variable_space,mu,sigma);
    Px = reshape(Px,length(y_domain),length(x_domain));
    Px = Px./max(max(Px));
    s = surf(x_domain,y_domain,zeros(size(Px)),Px);
    s.EdgeColor = 'none';
    hold on
    set(gca,'YDir','normal')
    view(2)
    axis square
    colormap(ax,brewermap([],'greys'))
    set(gca,'TickDir','out');
    ylim([0,225])
    ylabel('Slope (A/ms)')
    xlabel('Duration (ms)')
    title(CondLab{i});
end

% Plot conditional gaussians for each type of prior
for i = [1:3]
subplot(4,4,[11,12,15,16])
x_domain =  [-2500:2:2500];
y = normpdf(x_domain,prior.('mu|slow')(i)-175,sqrt(prior.('sigma2|slow')(i))); hold on
plot(x_domain,(y./max(y))./1.3+i,'r-','linewidth',2)
plot([prior.('mu|slow')(i)-175, prior.('mu|slow')(i)-175],[i, max((y./max(y))./1.3+i)],'r--','linewidth',1.5)

y = normpdf(x_domain,prior.('mu|fast')(i)-175,sqrt(prior.('sigma2|fast')(i)));
plot(x_domain,(y./max(y))./1.3+i,'b-','linewidth',2)
plot([prior.('mu|fast')(i)-175, prior.('mu|fast')(i)-175],[i, max((y./max(y))./1.3+i)],'b--','linewidth',1.5)

plot([0,3000],[i,i],'-','color',[.6,.6,.6])
title('Priors')
end
yticks([1:3])
ylim([0.75,4])
xlim([-2500,2500])
yticklabels({'Onset','Offset','Both(Joint)'})
xlabel('Duration Relative to Categorical Boundary')

%Plot spectrogram of vocalization
ax = subplot(4,5,[6,7]);
xf = 1:length(ex_voc);
Fs = 250000;
Time = xf/Fs;
spectrogram(ex_voc,128*15,120*15,128*15,Fs,'yaxis') % Spectrogram plot display
hold on
set(colorbar,'visible','on')  % Display color bar
colormap(ax,brewermap([],'Greys'))
xticks([0:3])
xlim([0 3])
yticks([])
ylim([0 80])
xlabel('Time (s)','Fontsize',12)
ylabel('Frequency (kHz)','fontsize',12)
caxis([-110 -90])        % set colorbar axis for contrast
box off
yticks([0:10:80])

%Plot raw signal of example vocalizations
subplot(4,5,[1,2])
plot(Time,ex_voc,'k-')
colorbar
set(colorbar,'visible','off')  % colorbar aligns spectrogram with signal
ylabel('Amplitude (a.u.)','fontsize',12)
xlabel('Time (s)','Fontsize',12)
xticks([0:1:3])
yticks([-1:1])
box off

% plot envelope of first example vocalization
a = 834.5;
Idx = 500;
subplot(4,4,[3,4])
max_env = max(ex_env(Time<1.75 & Time>.2));
plot(Time(Time<1.75 & Time>.2)-.2,ex_env(Time<1.75 & Time>.2)'./max_env,'-','linewidth',1.5,'color',[.2 .6 .6])
hold on
xlabel('Time (s)','Fontsize',12)

%plot onset, offset and half-max points onto envelope
tonex = df.ton(Idx)-a;
Aonex = df.Aon(Idx)./max_env;
thonex = df.thon(Idx)-a;
Ahonex = df.Ahon(Idx)./max_env;

toffex = df.toff(Idx)-a;
Aoffex = df.Aoff(Idx)./max_env;
thoffex = df.thoff(Idx)-a;
Ahoffex = df.Ahoff(Idx)./max_env;
for call=1:length(Idx)
plot([tonex(call),thonex(call),thoffex(call),toffex(call)],[Aonex(call),Ahonex(call),Ahoffex(call),Aoffex(call)],'ko-','linewidth',2)
end
xlim([0,1.5])
box off

% plot envelope of synthetic sounds
subplot(4,4,[7,8])
Xs = [1:length(Synthetic.Short)]./Synthetic.Fs;
plot(Xs(Xs<.5),Synthetic.Short(Xs<.5),'-','linewidth',2,'color',[.3 .7 .7]); hold on
plot(Xs(Xs<.5),Synthetic.Middle(Xs<.5),'-','linewidth',2,'color',[.2 .6 .6])
plot(Xs(Xs<.5),Synthetic.Long(Xs<.5),'-','linewidth',2,'color',[.1 .5 .5])
xlim([0,1.5])
ylim([0,1.1])
xlabel('Time (s)','Fontsize',12)
box off

clearvars -except prior
h=gcf;
set(h,'PaperSize',[20 10]); %set the paper size to what you want  
% print(h,'/results/figure_panel.pdf') % then print it

saveas(h,'/results/figure_panel.pdf')

function [condMu,condSigma] = condGauss(Mu,Sigma,a)

% conditional gaussian calculation (handles 2-3 variables)
if length(Mu) == 2
    condMu = Mu(1)+Sigma(1,2)/Sigma(2,2)*(a-Mu(2));
    condSigma = Sigma(1,1)-(Sigma(1,2)/Sigma(2,2))*Sigma(2,1);
elseif length(Mu) == 3
    
    Mu1 = Mu(1:2)+(Sigma(3,1:2)/Sigma(3,3))*(a-Mu(3));
    Sigma1 = Sigma(1:2,1:2)-(Sigma(1:2,3)/Sigma(3,3))*Sigma(3,1:2);
    
    condMu = Mu1(1)+(Sigma1(1,2)/Sigma1(2,2))*(a-Mu1(2));
    condSigma = Sigma1(1,1)-(Sigma1(1,2)/Sigma1(2,2))*Sigma1(2,1);
    
end

end

%%% function used for colormapping

function [map,num,typ,scheme] = brewermap(N,scheme)
% (c) 2014-2020 Stephen Cobeldick
persistent scm
%
raw = bmColors();
%
err = 'First input must be a real positive scalar numeric or [] or NaN.';
if nargin==0&&nargout==0
	hdr = {   'Type'; 'Scheme'; 'Nodes'};
	tsn = [{raw.typ};{raw.str};{raw.num}];
	fprintf('%-12s %-9s %s\n',hdr{:});
	fprintf('%-12s %-9s %u\n',tsn{:});
	return
elseif nargin==0 || isnumeric(N)&&isequal(N,[])
	% Default is the same as MATLAB colormaps:
	N = size(get(gcf,'colormap'),1);
	if nargin<2
		assert(~isempty(scm),'SC:colorbrewer:SchemeNotPreset',...
			'Scheme must be preset before this call: BREWERMAP(SCHEME)')
		scheme = scm;
	end
elseif nargin==1&&ischar(N)&&ndims(N)==2&&size(N,1)==1
	if strcmpi(N,'list')
		map = {raw.str};
		num = [raw.num];
		typ = {raw.typ};
		return
	end
	scheme = N; % preset
else
	assert(isnumeric(N)&&isscalar(N),...
		'SC:brewermap:NotScalarNumeric',err)
	assert(isnan(N)||isreal(N)&&isfinite(N)&&fix(N)==N&&N>=0,...
		'SC:brewermap:NotRealPositiveNotNaN',err)
end
%
assert(ischar(scheme)&&ndims(scheme)==2&&size(scheme,1)==1,...
	'SC:brewermap:NotCharacterVector',...
	'Second input must be a character vector (the scheme name).')
isr = strncmp(scheme,'*',1);
ids = strcmpi(scheme(1+isr:end),{raw.str});
assert(any(ids),'SC:brewermap:UnknownScheme','Unknown scheme name: %s',scheme)
%
num = raw(ids).num;
typ = raw(ids).typ;
%
if ischar(N)
	map = scm;
	scm = N;
	return
elseif N==0
	map = ones(0,3);
	return
elseif isnan(N)
	N = num;
end
%
% Downsample:
[idx,itp] = bmIndex(N,num,typ);
map= raw(ids).rgb(idx,:)/255;
% Interpolate:
if itp
	M = [... sRGB to XYZ
		0.4124564,0.3575761,0.1804375;...
		0.2126729,0.7151522,0.0721750;...
		0.0193339,0.1191920,0.9503041];
	wpt = [0.95047,1,1.08883]; % D65
	%
	map = bmRGB2Lab(map,M,wpt); % optional
	%
	% Extrapolate a small amount beyond end nodes:
	%ido = linspace(0,num+1,N+2);
	%ido = ido(2:end-1);
	% Interpolation completely within end nodes:
	ido = linspace(1,num,N);
	%
	switch typ
		case 'Diverging'
			mid = ceil(num/2);
			ida =   1:mid;
			idz = mid:num;
			map = [...
				interp1(ida,map(ida,:),ido(ido<=mid),'pchip');...
				interp1(idz,map(idz,:),ido(ido>mid),'pchip')];
		case 'Sequential'
			map = interp1(1:num,map,ido,'pchip');
		otherwise
			error('SC:brewermap:NoInterp','Cannot interpolate this type.')
	end
	%
	map = bmLab2RGB(map,M,wpt); % optional
end
% Limit output range:
map = max(0,min(1,map));
% Reverse row order:
if isr
	map = map(end:-1:1,:);
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%brewermap
function lab = bmRGB2Lab(rgb,M,wpt)
% Convert a matrix of sRGB values to Lab.
%applycform(rgb,makecform('srgb2lab','AdaptedWhitePoint',wpt))
% RGB2XYZ:
xyz = bmGammaInv(rgb) * M.';
% Remember to include my license when copying my implementation.
% XYZ2Lab:
xyz = bsxfun(@rdivide,xyz,wpt);
idx = xyz>(6/29)^3;
F = idx.*(xyz.^(1/3)) + ~idx.*(xyz*(29/6)^2/3+4/29);
lab(:,2:3) = bsxfun(@times,[500,200],F(:,1:2)-F(:,2:3));
lab(:,1) = 116*F(:,2) - 16;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%bmRGB2Lab
function rgb = bmGammaInv(rgb)
% Inverse gamma correction of sRGB data.
idx = rgb <= 0.04045;
rgb(idx) = rgb(idx) / 12.92;
rgb(~idx) = real(((rgb(~idx) + 0.055) / 1.055).^2.4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%bmGammaInv
function rgb = bmLab2RGB(lab,M,wpt)
% Convert a matrix of Lab values to sRGB.
%applycform(lab,makecform('lab2srgb','AdaptedWhitePoint',wpt))
% Lab2XYZ
tmp = bsxfun(@rdivide,lab(:,[2,1,3]),[500,Inf,-200]);
tmp = bsxfun(@plus,tmp,(lab(:,1)+16)/116);
idx = tmp>(6/29);
tmp = idx.*(tmp.^3) + ~idx.*(3*(6/29)^2*(tmp-4/29));
xyz = bsxfun(@times,tmp,wpt);
% Remember to include my license when copying my implementation.
% XYZ2RGB
rgb = max(0,min(1, bmGammaCor(xyz / M.')));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%cbLab2RGB
function rgb = bmGammaCor(rgb)
% Gamma correction of sRGB data.
idx = rgb <= 0.0031308;
rgb(idx) = 12.92 * rgb(idx);
rgb(~idx) = real(1.055 * rgb(~idx).^(1/2.4) - 0.055);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%bmGammaCor
function [idx,itp] = bmIndex(N,num,typ)
% Ensure exactly the same colors as the online ColorBrewer colorschemes.
%
itp = N>num;
switch typ
	case 'Qualitative'
		itp = false;
		idx = 1+mod(0:N-1,num);
	case 'Diverging'
		switch N
			case 1 % extrapolated
				idx = 8;
			case 2 % extrapolated
				idx = [4,12];
			case 3
				idx = [5,8,11];
			case 4
				idx = [3,6,10,13];
			case 5
				idx = [3,6,8,10,13];
			case 6
				idx = [2,5,7,9,11,14];
			case 7
				idx = [2,5,7,8,9,11,14];
			case 8
				idx = [2,4,6,7,9,10,12,14];
			case 9
				idx = [2,4,6,7,8,9,10,12,14];
			case 10
				idx = [1,2,4,6,7,9,10,12,14,15];
			otherwise
				idx = [1,2,4,6,7,8,9,10,12,14,15];
		end
	case 'Sequential'
		switch N
			case 1 % extrapolated
				idx = 6;
			case 2 % extrapolated
				idx = [4,8];
			case 3
				idx = [3,6,9];
			case 4
				idx = [2,5,7,10];
			case 5
				idx = [2,5,7,9,11];
			case 6
				idx = [2,4,6,7,9,11];
			case 7
				idx = [2,4,6,7,8,10,12];
			case 8
				idx = [1,3,4,6,7,8,10,12];
			otherwise
				idx = [1,3,4,6,7,8,10,11,13];
		end
	otherwise
		error('SC:brewermap:UnknownType','Unknown type string.')
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%bmIndex
function raw = bmColors()
% Return a structure of all colorschemes: name, scheme type, RGB values, number of nodes.
% Order: first sort by <typ>, then case-insensitive sort by <str>:

raw(23).str = 'Greys';
raw(23).typ = 'Sequential';
raw(23).rgb = [255,255,255;247,247,247;240,240,240;217,217,217;204,204,204;189,189,189;150,150,150;115,115,115;99,99,99;82,82,82;37,37,37;37,37,37;0,0,0];
% number of nodes:
for k = 1:numel(raw)
	switch raw(k).typ
		case 'Diverging'
			raw(k).num = 11;
		case 'Qualitative'
			raw(k).num = size(raw(k).rgb,1);
		case 'Sequential'
			raw(k).num = 9;
		otherwise
			error('SC:brewermap:UnknownType','Unknown type string.')
	end
end
%
end

