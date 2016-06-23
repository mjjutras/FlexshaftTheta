%%
BRnam = 'JN141024002';
Panddir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data\2014\JN_14_10_24\JN_14_10_24_12_36';

%%
datdir = 'C:\Data\VR_Blackrock\';

NS2 = openNSx(fullfile(datdir,[BRnam '.ns2']),'read');

chanLabels = {NS2.ElectrodesInfo.Label}';

LFPind = strncmp('A_chan',chanLabels,6) | ...
    strncmp('B_chan',chanLabels,6) | ...
    strncmp('C_chan',chanLabels,6);
eyeXind = strncmp('ainp1',chanLabels,5);
eyeYind = strncmp('ainp2',chanLabels,5);

%% pEpisode
shoulderMS = 500;
samplerate = 1000;
shoulder = round(shoulderMS*samplerate/1000); % shoulder in samples
freqs = (2^(1/8)).^(8:42);
width = 7;

chans = find(LFPind);
UV = cell(size(chans));
for chnlop = 1:length(chans)
    % get energy vector
    B = single(multienergyvec(double(NS2.Data(chans(chnlop),:)),freqs,samplerate,width));
    % calc the mean fit
    Blog = log10(double(B));
    Pm = mean(Blog,2);
    % get the fit
    % IMPORTANT: chi_squarefit assumes that frequencies are
    % logarithmically spaced!
    [all,R2] = chi_squarefit(freqs,Pm);
    all = all';
    % set the threshold
    thresh = all(:,951);
    % loop through each frequency and save the unionvector
    unionvector = single(zeros(size(freqs,2),size(B,2)));
    for f = 1:size(freqs,2)
      % get the pepisode
      unionvector(f,:) = single(episodeid(B(f,:),thresh(f),3*samplerate/freqs(f),shoulder,0));
    end
    UV{chnlop} = unionvector;
end

datdir = 'R:\Buffalo Lab\Mike\VirtualNavigationProject\MATFiles\pEpisode\Flexshaft_JN2014implant';
save(fullfile(datdir,'JN141024002_NS2_foraging_pEpi_160622.mat'),'UV','-v7.3')

%% plot pEpisode

% Array A
figure
% 1
subplot(3,4,9); imagesc(squeeze(nanmean(UVt{1}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 1')
% 2
subplot(3,4,5); imagesc(squeeze(nanmean(UVt{2}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 2')
% 3
subplot(3,4,1); imagesc(squeeze(nanmean(UVt{3}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 3')
% 4
subplot(3,4,10); imagesc(squeeze(nanmean(UVt{4}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 4')
% 5
subplot(3,4,6); imagesc(squeeze(nanmean(UVt{5}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 5')
% 6
subplot(3,4,2); imagesc(squeeze(nanmean(UVt{6}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 6')
% 7
subplot(3,4,11); imagesc(squeeze(nanmean(UVt{7}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 7')
% 8
subplot(3,4,7); imagesc(squeeze(nanmean(UVt{8}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 8')
% 9
subplot(3,4,3); imagesc(squeeze(nanmean(UVt{9}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 9')
% 10
subplot(3,4,12); imagesc(squeeze(nanmean(UVt{10}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 10')
% 11
subplot(3,4,8); imagesc(squeeze(nanmean(UVt{11}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 11')
% 12
subplot(3,4,4); imagesc(squeeze(nanmean(UVt{12}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 12')

% run following line after maximizing figure
suplabel('Array A','t')


% Array B
figure
% 1
subplot(3,4,9); imagesc(squeeze(nanmean(UVt{13}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 1')
% 2
subplot(3,4,5); imagesc(squeeze(nanmean(UVt{14}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 2')
% 3
subplot(3,4,1); imagesc(squeeze(nanmean(UVt{15}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 3')
% 4
subplot(3,4,10); imagesc(squeeze(nanmean(UVt{16}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 4')
% 5
subplot(3,4,6); imagesc(squeeze(nanmean(UVt{17}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 5')
% 6
subplot(3,4,2); imagesc(squeeze(nanmean(UVt{18}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 6')
% 7
subplot(3,4,11); imagesc(squeeze(nanmean(UVt{19}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 7')
% 8
subplot(3,4,7); imagesc(squeeze(nanmean(UVt{20}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 8')
% 9
subplot(3,4,3); imagesc(squeeze(nanmean(UVt{21}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 9')
% 10
subplot(3,4,12); imagesc(squeeze(nanmean(UVt{22}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 10')
% 11
subplot(3,4,8); imagesc(squeeze(nanmean(UVt{23}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 11')
% 12
subplot(3,4,4); imagesc(squeeze(nanmean(UVt{24}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 12')

% run following line after maximizing figure
suplabel('Array B','t')


% Array C
figure
% 1
subplot(3,4,9); imagesc(squeeze(nanmean(UVt{25}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 1')
% 2
subplot(3,4,5); imagesc(squeeze(nanmean(UVt{26}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 2')
% 3
subplot(3,4,1); imagesc(squeeze(nanmean(UVt{27}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 3')
% 4
subplot(3,4,10); imagesc(squeeze(nanmean(UVt{28}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 4')
% 5
subplot(3,4,6); imagesc(squeeze(nanmean(UVt{29}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 5')
% 6
subplot(3,4,2); imagesc(squeeze(nanmean(UVt{30}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 6')
% 7
subplot(3,4,11); imagesc(squeeze(nanmean(UVt{31}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 7')
% 8
subplot(3,4,7); imagesc(squeeze(nanmean(UVt{32}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 8')
% 9
subplot(3,4,3); imagesc(squeeze(nanmean(UVt{33}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 9')
% 10
subplot(3,4,12); imagesc(squeeze(nanmean(UVt{34}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 10')
% 11
subplot(3,4,8); imagesc(squeeze(nanmean(UVt{35}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 11')
% 12
subplot(3,4,4); imagesc(squeeze(nanmean(UVt{36}(:,:,1:5000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 12')

% run following line after maximizing figure
suplabel('Array C','t')

%%
fltord = 40;
lowpasfrq = 40;
lim = 70000;
artpadding = 0.01;
Fs = 1/fs; % sampling rate in Hz

clear xind yind
for k=1:length(NS2.ElectrodesInfo)
    if strmatch('ainp1',NS2.ElectrodesInfo(k).Label)
        xind = k;
    elseif strmatch('ainp2',NS2.ElectrodesInfo(k).Label)
        yind = k;
    end
end

% data_eye.artifact = cell(size(data_eye.time));
e_x = double(NS2.Data(xind,:));
e_y = double(NS2.Data(yind,:));
% e_x = double(NS2.Data{1}(xind,:));
% e_y = double(NS2.Data{1}(yind,:));

%low pass filter the eye position data_eye
nyqfrq = Fs ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);
e_x_lowpassfilter=filtfilt(flt,1, e_x);
e_y_lowpassfilter=filtfilt(flt,1, e_y);

%differentiate and multiply with sampling rate to get velocity as deg/sec
x_vF = diff(e_x_lowpassfilter) .* Fs;
y_vF = diff(e_y_lowpassfilter) .* Fs;

% % differentiate and multiply with sampling rate without filtering
% % this gives the eye velocity in the horizontal and vertical domains
x_v = diff(e_x) .* Fs;
y_v = diff(e_y) .* Fs;

% combine x- and y-velocity to get eye velocity in degrees/second
vel = abs(complex(x_v,y_v));
velF = abs(complex(x_vF,y_vF));

%detect saccade begins and saccade ends
sacbeg = find(diff(velF > lim) > 0);
sacend = find(diff(velF > lim) < 0);
if velF(end)>lim
    sacbeg = sacbeg(1:end-1); % changed this line from artifact_xysaccade120420.m
end
if velF(1)>lim
    sacend = sacend(2:end);
end

if size(sacbeg,1)
    sacbeg=sacbeg';
    sacend=sacend';
end
    
artifact = round([sacbeg(:) - artpadding*Fs sacend(:) + artpadding*Fs]);

% % find saccade start/end in task time
% sacdum = [NS2_timestamp(artifact(:,1))' NS2_timestamp(artifact(:,2))'];

% merge artifacts when inter-artifact interval is less than 10 ms
sacarr = [];
iai=(artifact(2:end,1)-artifact(1:end-1,2))>0.01; sacdumcol1=artifact(2:end,1); sacdumcol2=artifact(1:end-1,2);
sacarr(:,1) = [artifact(1,1); sacdumcol1(iai,1)];
sacarr(:,2) = [sacdumcol2(iai,1); artifact(end,2)];

% figure;plot(velF);hold on;for k=1:100;line([sacarr(k,1) sacarr(k,1)],ylim,'Color','g');line([sacarr(k,2) sacarr(k,2)],ylim,'Color','r');end
