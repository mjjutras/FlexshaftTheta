% This code gets trial-based data from VR foraging task
% Data streams are aligned using eye movement data, calculates pEpisode
% across trials, keeps the neural data in a Fieldtrip-friendly structure

%%
BRnam = 'JN141024002';
lognam = 'JN_14_10_24_12_36_trldat.mat';

datdir = 'C:\Data\VR_Blackrock\';
trldatdir = 'R:\Buffalo Lab\Mike\VirtualNavigationProject\MATFiles\trldat';

%% load Blackrock data
NS2 = openNSx(fullfile(datdir,[BRnam '.ns2']),'read');

chanLabels = {NS2.ElectrodesInfo.Label}';

LFPind = strncmp('A_chan',chanLabels,6) | ...
    strncmp('B_chan',chanLabels,6) | ...
    strncmp('C_chan',chanLabels,6);
eyeXind = strncmp('ainp1',chanLabels,5);
eyeYind = strncmp('ainp2',chanLabels,5);

%% pEpisode (only run once, saves results)
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

UVdir = 'R:\Buffalo Lab\Mike\VirtualNavigationProject\MATFiles\pEpisode\Flexshaft_JN2014implant';
save(fullfile(UVdir,'JN141024002_NS2_foraging_pEpi_160622.mat'),'UV','-v7.3')

%% align data streams using eyetracking data
MATdir = 'C:\Data\MAT';
load(fullfile(MATdir,'JN141024002_NS2_foraging_pEpi_160622.mat'))

clear trial sampleinfo
c = 1;
UVt = cell(size(UV));
for trllop = 1:length(trldat.time)
    
    fprintf('Processing trial %d of %d', trllop, length(trldat.time));
    fprintf('\n')

    if trllop==1
        [Xs1, lag] = xcorr(double(NS2.Data(eyeXind,:)), trldat.eyedat{trllop}(1,:));
        Xs2 = xcorr(double(NS2.Data(eyeXind,:)), trldat.eyedat{trllop}(1,300:end));
        Xs3 = xcorr(double(NS2.Data(eyeYind,:)), trldat.eyedat{trllop}(2,:));
        Xs4 = xcorr(double(NS2.Data(eyeYind,:)), trldat.eyedat{trllop}(2,300:end));
        [~,Is1] = max(Xs1.*(Xs1-Xs2).*Xs3.*(Xs3-Xs4));
        
        [~,locs] = findpeaks(Xs1.*Xs3,1000);
        [~,Is2] = min(abs(lag(round(locs*1000))-lag(Is1)));
        
        trl_start = lag(round(locs(Is2)*1000))+1;
        trl_end = trl_start+size(trldat.eyedat{trllop},2)-1;
    else
        lagind1 = trl_start;
        lagind2 = trl_end+(size(trldat.eyedat{trllop}(1,:),2)*2);
        [Xs1, lag] = xcorr(double(NS2.Data(eyeXind,lagind1:min([size(NS2.Data,2) lagind2]))), trldat.eyedat{trllop}(1,:));
        Xs2 = xcorr(double(NS2.Data(eyeYind,lagind1:min([size(NS2.Data,2) lagind2]))), trldat.eyedat{trllop}(2,:));
        [~,Is1] = max(Xs1.*Xs2);
        
        trl_start = lag(Is1)+lagind1-1;
        trl_end = trl_start+size(trldat.eyedat{trllop},2)-1;
    end
    
    % Save as a new file structure containing all the necessary info.
    trial{c} = double(NS2.Data(1:end-2,trl_start:trl_end));
    
    % append the eye data from the Python log file
    trial{c}(end+1:end+2,:) = trldat.eyedat{trllop};
    
    % include both Blackrock and Python analog data for now, until we
    % decide which version has cleaner A2D conversion
    
    % corresponds to samples in NS2 file
    sampleinfo(c,:) = [trl_start trl_end];
    
    for chnlop = 1:length(UV)
        dum = UV{chnlop}(:,trl_start:trl_end);
        if c==1
            UVt{chnlop}(c,:,:) = dum;
        else
            if size(UVt{chnlop},3)>size(dum,2)
                UVt{chnlop}(c,:,:) = cat(2,dum,nan(size(dum,1),size(UVt{chnlop},3)-size(dum,2)));
            elseif size(UVt{chnlop},3)<size(dum,2)
                UVt{chnlop} = cat(3,UVt{chnlop},nan(size(UVt{chnlop},1),size(dum,1),size(dum,2)-size(UVt{chnlop},3)));
                UVt{chnlop}(c,:,:) = dum;
            end
        end
    end
    
    c = c+1;
    
end

% plot eye data trial by trial to ensure accuracy
eyexerr = nan(size(trial)); % error between normalized eyetraces per trial
for trllop = 1:length(trial)

    xp = trial{trllop}(end-1,:);
    xb = trial{trllop}(end-3,:);
    
    xpn = xp/(max(xp)-min(xp));
    xbn = xb/(max(xb)-min(xb));
    
    figure
    plot([xpn; xbn]')
%     pause;close

    eyexerr(trllop) = mean(xpn-xbn);

end

% create the data structure for Fieldtrip
clear data
data.trial = trial;
data.sampleinfo = sampleinfo;
for trllop = 1:size(sampleinfo,1)
    data.time{trllop} = 0:0.001:(diff(sampleinfo(trllop,:),1,2)/1000);
end
data.fsample = 1000;
data.label = {'A01'; 'A02'; 'A03'; 'A04'; 'A05'; 'A06'; 'A07'; ...
    'A08'; 'A09'; 'A10'; 'A11'; 'A12'; 'B01'; 'B02'; 'B03'; 'B04'; ...
    'B05'; 'B06'; 'B07'; 'B08'; 'B09'; 'B10'; 'B11'; 'B12'; 'C01'; ...
    'C02'; 'C03'; 'C04'; 'C05'; 'C06'; 'C07'; 'C08'; 'C09'; 'C10'; ...
    'C11'; 'C12'; 'eyeXBlk'; 'eyeYBlk'; 'eyeXPyt'; 'eyeYPyt'};
clear trial sampleinfo

NSdir = 'R:\Buffalo Lab\Mike\VirtualNavigationProject\MATFiles\NSdat\';
save(fullfile(NSdir,'JN141024002_NS2_foraging_NSdat.mat'),'data')
UVdir = 'R:\Buffalo Lab\Mike\VirtualNavigationProject\MATFiles\pEpisode\Flexshaft_JN2014implant';
save(fullfile(UVdir,'JN141024002_NS2_foraging_pEpiTrial_160628.mat'),'UVt','-v7.3')

%% plot pEpisode

% Array A
figure
% 1
subplot(3,4,9); imagesc(squeeze(nanmean(UVt{1}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 1')
% 2
subplot(3,4,5); imagesc(squeeze(nanmean(UVt{2}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 2')
% 3
subplot(3,4,1); imagesc(squeeze(nanmean(UVt{3}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 3')
% 4
subplot(3,4,10); imagesc(squeeze(nanmean(UVt{4}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 4')
% 5
subplot(3,4,6); imagesc(squeeze(nanmean(UVt{5}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 5')
% 6
subplot(3,4,2); imagesc(squeeze(nanmean(UVt{6}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 6')
% 7
subplot(3,4,11); imagesc(squeeze(nanmean(UVt{7}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 7')
% 8
subplot(3,4,7); imagesc(squeeze(nanmean(UVt{8}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 8')
% 9
subplot(3,4,3); imagesc(squeeze(nanmean(UVt{9}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 9')
% 10
subplot(3,4,12); imagesc(squeeze(nanmean(UVt{10}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 10')
% 11
subplot(3,4,8); imagesc(squeeze(nanmean(UVt{11}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 11')
% 12
subplot(3,4,4); imagesc(squeeze(nanmean(UVt{12}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 12')

% run following line after maximizing figure
suplabel('Array A','t')


% Array B
figure
% 1
subplot(3,4,9); imagesc(squeeze(nanmean(UVt{13}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 1')
% 2
subplot(3,4,5); imagesc(squeeze(nanmean(UVt{14}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 2')
% 3
subplot(3,4,1); imagesc(squeeze(nanmean(UVt{15}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 3')
% 4
subplot(3,4,10); imagesc(squeeze(nanmean(UVt{16}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 4')
% 5
subplot(3,4,6); imagesc(squeeze(nanmean(UVt{17}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 5')
% 6
subplot(3,4,2); imagesc(squeeze(nanmean(UVt{18}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 6')
% 7
subplot(3,4,11); imagesc(squeeze(nanmean(UVt{19}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 7')
% 8
subplot(3,4,7); imagesc(squeeze(nanmean(UVt{20}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 8')
% 9
subplot(3,4,3); imagesc(squeeze(nanmean(UVt{21}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 9')
% 10
subplot(3,4,12); imagesc(squeeze(nanmean(UVt{22}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 10')
% 11
subplot(3,4,8); imagesc(squeeze(nanmean(UVt{23}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 11')
% 12
subplot(3,4,4); imagesc(squeeze(nanmean(UVt{24}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 12')

% run following line after maximizing figure
suplabel('Array B','t')


% Array C
figure
% 1
subplot(3,4,9); imagesc(squeeze(nanmean(UVt{25}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 1')
% 2
subplot(3,4,5); imagesc(squeeze(nanmean(UVt{26}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 2')
% 3
subplot(3,4,1); imagesc(squeeze(nanmean(UVt{27}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 3')
% 4
subplot(3,4,10); imagesc(squeeze(nanmean(UVt{28}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 4')
% 5
subplot(3,4,6); imagesc(squeeze(nanmean(UVt{29}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 5')
% 6
subplot(3,4,2); imagesc(squeeze(nanmean(UVt{30}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 6')
% 7
subplot(3,4,11); imagesc(squeeze(nanmean(UVt{31}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 7')
% 8
subplot(3,4,7); imagesc(squeeze(nanmean(UVt{32}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 8')
% 9
subplot(3,4,3); imagesc(squeeze(nanmean(UVt{33}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 9')
% 10
subplot(3,4,12); imagesc(squeeze(nanmean(UVt{34}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 10')
% 11
subplot(3,4,8); imagesc(squeeze(nanmean(UVt{35}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 11')
% 12
subplot(3,4,4); imagesc(squeeze(nanmean(UVt{36}(:,:,1:50000),1)))
axis xy; set(gca,'YTickLabel',num2str(freqs(5:5:35)')); colormap hot
title('Ch. 12')

% run following line after maximizing figure
suplabel('Array C','t')

%%
fltord = 40;
lowpasfrq = 40;
lim = 70000;
artpadding = 0.01;
Fs = 1000; % sampling rate in Hz

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
