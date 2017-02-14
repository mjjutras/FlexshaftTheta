%% uses code from Seth's ImportVREyeData function

%%
% 12:20 manual calibration 001
% 12:23 starting w/ photos + calibration (manual first)
% 12:35 Starting in square. 002

% images from this set are located in these folders
% R:\Buffalo Lab\eblab\Cortex Programs\SCM Relational Memory Project\Image Sets\SCM01
% or
% R:\Buffalo Lab\eblab\Cortex Programs\List Relational Memory Project\Image Sets

img_dir = 'R:\Buffalo Lab\eblab\Cortex Programs\SCM Relational Memory Project\Image Sets\SCM01';

%%
% % load the Blackrock NS6 file
% NS6 = openNSx('C:\Data\VR_Blackrock\JN141024001.ns6');
% load the Blackrock NS2 file
NS2 = openNSx('C:\Data\Blackrock_VR\JN141024001.ns2');

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
save(fullfile(datdir,'JN141024001_NS2_picviewing_pEpi_160616.mat'),'UV','-v7.3')

%% ID saccades
fs = 1/NS2.MetaTags.SamplingFreq;
fltord = 40;
lowpasfrq = 40;
lim = 30000;
artpadding = 0.01;
Fs = 1/fs; % sampling rate in Hz

%low pass filter the eye position data_eye
nyqfrq = Fs ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);
e_x_lowpassfilter = filtfilt(flt,1, double(NS2.Data(eyeXind,:)));
e_y_lowpassfilter = filtfilt(flt,1, double(NS2.Data(eyeYind,:)));

%differentiate and multiply with sampling rate to get velocity as deg/sec
x_vF = diff(e_x_lowpassfilter) .* Fs;
y_vF = diff(e_y_lowpassfilter) .* Fs;

% combine x- and y-velocity to get eye velocity in degrees/second
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

%%
% load the eye calibration file
fildir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data\2014\JN_14_10_24';
eyecal = 'eye_cal2_14_10_24_12_23';
timcal = 'time_cal2_14_10_24_12_23';

[eye1, eye2, eye3] = textread(fullfile(fildir,eyecal), ...
    '%f%f%f','headerlines',1,'delimiter',',','whitespace','\n');

try %sometimes eye data file gets cutoff at end?
    eye = [eye1 eye2 eye3];
catch
    eye = [eye1(1:end-1) eye2 eye3];
end

t0 = eye(1,1); %time 0
eye(:,1) = eye(:,1)-t0; %make time relative to 0

% event codes:
% 1: square_position, 2: square moved/on, 3: reward,
% 4: break fixation, 5: clrchng, 6: square off,
% 7: image on, 8: image off, 9: image_name

time = [];
square_pos = [];
events = [];
image_name = {};
img_count = 1;

fid = fopen(fullfile(fildir,timcal));
linecount = 1;
tline = fgets(fid); %skip line 1
tline = fgets(fid); %skip line 2
while ischar(tline)
    tline = fgets(fid); %get next line
    if tline ~= -1 %end of file is noted by -1
        C = textscan(tline,'%s'); %parse line by spaces and comas
        time(linecount,1) = str2double(C{1}{1});
        if strcmpi(C{1}{2},'Square')
            if strcmpi(C{1}{3},'Position,') %gives position of square in pixels
                events(linecount,1) = 1;
                square_pos(linecount,1) = str2double(C{1}{4});
                square_pos(linecount,2) = str2double(C{1}{5});
            elseif strcmpi(C{1}{3},'moved') || strcmpi(C{1}{3},'on') %square moved to different location and displayed
                events(linecount,1) = 2;
            elseif strcmpi(C{1}{3},'dims') %color change occured
                events(linecount,1) = 5;
            elseif strcmpi(C{1}{3},'off') %square diseappears
                events(linecount,1) = 6;
            else
                disp('unknown events parameter in column 3 of event array') %so we known we didn't miss anything
            end
        elseif strcmpi(C{1}{2},'reward') %monkey received reward
            events(linecount,1) = 3;
        elseif strcmpi(C{1}{2},'no') %trial aborted due to break fixation, etc.
            events(linecount,1) = 4;
        elseif ~isempty(strfind(lower(C{1}{2}),'photo')) %image trials
            if length(C{1}) == 3
                if  strcmpi(C{1}{3},'on') %image on
                    events(linecount,1) = 7;
                elseif strcmpi(C{1}{3},'off') %image off
                    events(linecount,1) = 8;
                end
            else
                if~isempty(strfind(C{1}{2},'/')) %image name
                    events(linecount,1) = 9;
                    slash = strfind(C{1}{2},'/');
                    image_name{img_count,1} = C{1}{2}(slash+1:end);
                    img_count = img_count+1;
                else
                    disp('Unknown parameter for column 2 for picture trials')  %so we known we didn't miss anything
                end
            end
        end
        linecount = linecount+1; %to go to the next line
    end
end
fclose(fid)

time = time-t0; %use t0 to ensure time is the same but should be

%---Get Set number for images...for later use---%
if strcmpi(image_name{1}(1),'S')
    setnum = str2double(image_name{1}(2:3)); %get the set number, assumes set number is the same for all images
else %VR pilot sets
    setnum = str2double(image_name{1}(4)); %get the set number, assumes set number is the same for all images
end

%% calibrate (don't really need to do this right now)

% %---Get the calibration Function---%
% %doing for individual sessions of calibration. Currently assumes stable
% %calibration within short session. Does not track drift
% 
% % get unique calibration points
% xsquare_poses = unique(square_pos(:,1)); xsquare_poses(isnan(xsquare_poses)) = [];
% ysquare_poses = unique(square_pos(:,2)); ysquare_poses(isnan(ysquare_poses)) = [];
% 
% rewards = find(events == 3); %find where a reward occured
% caldat_x = cell(length(xsquare_poses),length(ysquare_poses));
% caldat_y = cell(length(xsquare_poses),length(ysquare_poses));
% square_dims = find(events == 5);
% square_pos_events = find(events == 1);
% for rew = 1:length(rewards)
%     rewardind = rewards(rew); %index of event array in which reward was delivered
%     
%     square_change = square_dims(square_dims < rewardind); %find the last time the color change
%     square_change = square_change(end);
%     clrchngtim = time(square_change); %look at time in which square changes color (dims)
%     start_data_collection = clrchngtim - 0.5; %500 eyes
%     
%     %find indices associated with times between 500 ms before color change
%     %and color change
%     eyeind = find(eye(:,1) > start_data_collection & eye(:,1) <= clrchngtim);
%     
%     %determine x and y eye position for those 500 ms
%     xsquare_pos = eye(eyeind,2);
%     ysquare_pos = eye(eyeind,3);
%     
%     %position of square
%     sq_pos = square_pos_events(square_pos_events < rewardind);
%     sq_pos = sq_pos(end);
%     square_x = square_pos(sq_pos,1);
%     square_y = square_pos(sq_pos,2);
%     
%     xind = find(xsquare_poses == square_x);
%     yind = find(ysquare_poses == square_y);
%     
%     caldat_x{xind,yind} = [caldat_x{xind,yind} mean(xsquare_pos)];
%     caldat_y{xind,yind} = [caldat_y{xind,yind} mean(ysquare_pos)];
% end
% 
% % control are unique combinations x and y position of square in pixel coordinates
% control_x = [];
% control_y = [];
% for xi = 1:length(xsquare_poses);
%     for yi = 1:length(ysquare_poses);
%         control_x = [control_x; 1.5*xsquare_poses(xi)];
%         control_y = [control_y; ysquare_poses(yi)];
%     end
% end
% 
% % input are estimated position of the color change squares in mV (the
% % output of eyescan).
% input_x = [];
% input_y = [];
% for xi = 1:length(xsquare_poses);
%     for yi = 1:length(ysquare_poses);
%         input_x = [input_x; nanmean(caldat_x{xi,yi})];
%         input_y = [input_y; nanmean(caldat_y{xi,yi})];
%     end
% end
% 
% % tform automatically transforms estimated eye positions in mV into pixel
% % coordinates in a "magical" way.
% tform = cp2tform([control_x control_y], [input_x input_y],'polynomial',4);
% tform.forward_fcn = tform.inverse_fcn;
% 
% %---Create Plot of Calibration Quality---%
% %Visualize position of estimated square position (blue *) compared to
% %actual position of color change squares (red +)
% [xp,yp] = tformfwd(tform,input_x,input_y);%transform inputs from mV to pixels to plot
% figure
% hold on
% for xi = 1:length(xsquare_poses);
%     for yi = 1:length(ysquare_poses);
%         plot(1.5*xsquare_poses(xi),ysquare_poses(yi),'r+')
%     end
% end
% for xi = 1:length(xp);
%     plot(xp(xi),yp(xi),'b*')
% end
% hold off
% % save_and_close_fig(set_dir,['Calibration for ' data_file])
% 
% %---Recalibrate Collected Eye Data---%
% eyedat = {};
% img_on = find(events == 7);
% img_off = find(events == 8);
% for img = 1:length(image_name)
%     eye_ind = find(eye(:,1) > time(img_on(img)) & eye(:,1) <= time(img_off(img)));
%     
%     [xp,yp] = tformfwd(tform,eye(eye_ind,2)',eye(eye_ind,3)');%transform inputs from mV to pixels to plot
%     eyedat{img} = [xp;yp];
% end

%% get image numbers and presentation order

%---Get image numbers and the order they appeared in---%
imgnum = NaN(1,36);
if strcmpi(image_name{1}(1),'S')
    for img = 1:length(image_name);
        imgnum(img) = str2double(image_name{img}(5:6));
    end
else %VR pilot sets
    for img = 1:length(image_name);
        imgnum(img) = str2double(image_name{img}(6:7));
    end
end

imgnum(imgnum == 0) = NaN;%don't know why zero shows up

pairings = NaN(2,36);
for img = min(imgnum):max(imgnum) %assumes images don't cross over sets
    ind = find(imgnum == img);
    if length(ind) == 2 %otherwise don't care
        pairings(1,img) = ind(1); %novel
        pairings(2,img) = ind(2); %repeat
    end
end

%% plot eye data (need to run calibration to do this)

% %---Plot Eye Data on Corresponding Figure---%
% for img = 1:size(pairings,2);
%     if ~isnan(pairings(1,img))
%         im1 = imread(fullfile(img_dir,image_name{pairings(1,img)}));
%         im2 = imread(fullfile(img_dir,image_name{pairings(2,img)}));
%         
%         figure
%         subplot(1,2,1)
%         imshow(im1)
%         hold on
%         plot(eyedat{pairings(1,img)}(1,:)+400,...
%             600-(eyedat{pairings(1,img)}(2,:)+300),'y')
%         hold off
%         title('Novel Image')
%         
%         subplot(1,2,2)
%         imshow(im2)
%         hold on
%         plot(eyedat{pairings(2,img)}(1,:)+400,...
%             600-(eyedat{pairings(2,img)}(2,:)+300),'y')
%         hold off
%         title('Repeat Image')
%         
%         suplabel(['Set ' num2str(setnum) ' Image ' num2str(pairings(1,img)) ...
%             ' Spacing ' num2str(pairings(2,img)-pairings(1,img))],'t',[.08 .08 .84 .74]);
%         
% %         save_and_close_fig(set_dir,['Set_' num2str(setnum) '_img_' ...
% %             num2str(pairings(1,img)) '_' data_file]);
%     end
% end
% 
% % it looks like this might remove data that fall outside of the image
% % border; commenting out for the purposes of this analysis
% % for eye = 1:length(eyedat)
% %     x = eyedat{eye}(1,:);
% %     y = eyedat{eye}(2,:);
% %     x = x+400;
% %     y = y+300;
% %     
% %     y(x < -50) = [];
% %     x(x < -50) = [];
% %     x(y < -50) = [];
% %     y(y < -50) = [];
% %     
% %     y(x > 850) = [];
% %     x(x > 850) = [];
% %     x(y > 650) = [];
% %     y(y > 650) = [];
% %     
% %     eyedat{eye} = [x;y];
% % end
% 
% % fixationstats = ClusterFixation_Short(eyedat);
% % save([eyedata_dir monk  data_file '-fixation.mat'],'fixationstats',...
% %     'image_name','pairings','imgnum','setnum')

%% get eye data and interpolate for alignment with neural signals

% eye time (first column) is only resolved down to 10 ms, let's fix this
eyetimeReal = eye(1,1):eye(end,1)/size(eye,1):eye(end,1);
%  make the lengths match by dropping extra elements
eyetimeReal = eyetimeReal(1:size(eye,1));

eyedatPyt = {};
eyetimePyt = {};
img_on = find(events == 7);
img_off = find(events == 8);
for img = 1:length(image_name)
    eye_ind = find(eyetimeReal > time(img_on(img)) & eyetimeReal <= time(img_off(img)));
    eyedatPyt{img} = eye(eye_ind,2:3)';
    eyetimePyt{img} = eyetimeReal(eye_ind);
end

% interpolate eye data
eyedatInterp = cell(1,length(eyedatPyt));
for trllop = 1:length(eyedatPyt)
    eyetime_interp = (round(eyetimePyt{trllop}(1)*1000):round(eyetimePyt{trllop}(end)*1000))/1000;
    eye_interp = nan(2,length(eyetime_interp));
    for k=1:length(eyetimePyt{trllop})
        eye_interp(:,ft_nearest(eyetime_interp,eyetimePyt{trllop}(k))) = eyedatPyt{trllop}(:,k);
    end
    for ii = find(~isnan(eye_interp(1,:)))
        I = eye_interp(:,ii);
        subind = ii+1:find(~isnan(eye_interp(1,ii+1:end)),1,'first')+ii-1;
        eye_interp(:,subind) = repmat(I,1,length(subind));
    end
    eyedatInterp{trllop} = eye_interp;
end

%% align data streams using eyetracking data

% cross-correlate entire eye signal (after downsampling NS2 signals)
% to get a rough idea of where the lag is
% sampling rate in Python log should be 240
% use this to check: round(1/median(diff(eyetimePyt{1})))
eyeXdownsamp = resample(double(NS2.Data(eyeXind,:)),240,1000);
[Xs1, lag] = xcorr(eyeXdownsamp, eye(:,2));
[~,Is1] = max(Xs1);
eyedatLag = lag(Is1)*(1000/240); % lag in ms
firstTrialLag = round(eyedatLag+(eyetimePyt{1}(1)*1000)); % approx. lag of first trial

clear trial sampleinfo
c = 1;
UVt = cell(size(UV));
for img = 1:size(pairings,2)
    if ~isnan(pairings(1,img))
        % eyedat{pairings(1,img)} contains eyedata for first viewing
        
        if pairings(1,img)==1
            [Xs1, lag] = xcorr(double(NS2.Data(eyeXind,firstTrialLag-1000:firstTrialLag+1000+length(eyedatInterp{pairings(1,img)}(1,:)))), eyedatInterp{pairings(1,img)}(1,:));
            %                 Xs2 = xcorr(double(NS2.Data(eyeXind,firstTrialLag-1000:firstTrialLag+1000+length(eyedatInterp{pairings(1,img)}(1,:)))), eyedatInterp{pairings(1,img)}(1,300:end));
            Xs3 = xcorr(double(NS2.Data(eyeYind,firstTrialLag-1000:firstTrialLag+1000+length(eyedatInterp{pairings(1,img)}(1,:)))), eyedatInterp{pairings(1,img)}(2,:));
            %                 Xs4 = xcorr(double(NS2.Data(eyeYind,firstTrialLag-1000:firstTrialLag+1000+length(eyedatInterp{pairings(1,img)}(1,:)))), eyedatInterp{pairings(1,img)}(2,300:end));
            %                 [~,Is1] = max(Xs1.*(Xs1-Xs2).*Xs3.*(Xs3-Xs4));
            [~,Is1] = max(Xs1.*Xs3);
            
            [~,locs] = findpeaks(Xs1.*Xs3,1000);
            [~,Is2] = min(abs(lag(round(locs*1000))-lag(Is1)));
            
            trl_start = lag(round(locs(Is2)*1000))+firstTrialLag-999;
            trl_end = trl_start+length(eyedatInterp{pairings(1,img)})-1;
        else
            lagind1 = trl_end+round(((eyetimePyt{pairings(1,img)}(1)-1-eyetimePyt{pairings(1,img)-1}(end)))*1000);
            lagind2 = trl_end+round(((eyetimePyt{pairings(1,img)}(end)+2-eyetimePyt{pairings(1,img)-1}(end)))*1000);
            [Xs1, lag] = xcorr(double(NS2.Data(eyeXind,lagind1:min([size(NS2.Data,2) lagind2]))), eyedatInterp{pairings(1,img)}(1,:));
            Xs2 = xcorr(double(NS2.Data(eyeYind,lagind1:min([size(NS2.Data,2) lagind2]))), eyedatInterp{pairings(1,img)}(2,:));
            [~,Is1] = max(Xs1.*Xs2);
            
            trl_start = lag(Is1)+lagind1-1;
            trl_end = trl_start+length(eyedatInterp{pairings(1,img)})-1;
        end
        
        % Save as a new file structure containing all the necessary info.
        trial{c} = double(NS2.Data(1:end-2,trl_start:trl_end));
        
        % append the eye data from the Python log file
        trial{c}(end+1:end+2,:) = eyedatInterp{pairings(1,img)};
        
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
end

% plot eye data trial by trial to ensure accuracy
for trllop = 1:length(trial)

    xp = trial{trllop}(end-1,:);
    xb = trial{trllop}(end-3,:);
    
    xpn = xp/(max(xp)-min(xp));
    xbn = xb/(max(xb)-min(xb));
    
    figure
    plot([xpn; xbn]')
%     pause;close
    
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
sactrl = [];
presac = 200;
postsac = 250;
for trllop = 1:length(data.trial)
    % start and stop times for saccades in this trial
    trlsacdef = sacarr(sacarr(:,1)>data.sampleinfo(trllop,1) & sacarr(:,2)<data.sampleinfo(trllop,2),:);
    for saclop = 1:size(trlsacdef,1)
        % make sure pre- and post-saccade period doesn't coincide with
        % start/end of trial
        if (trlsacdef(saclop,1)-presac)>=data.sampleinfo(trllop,1) ...
                && (trlsacdef(saclop,1)+postsac)<=data.sampleinfo(trllop,2)
%             % make sure pre-saccade period doesn't overlap with the
%             % previous saccade
%             if saclop == 1 || ...
%                     trlsacdef(saclop,1)-trlsacdef(saclop-1,2)>presac
                % make sure post-saccade period doesn't overlap with the next
                % saccade
                if saclop == size(trlsacdef,1) || ...
                        trlsacdef(saclop+1,1)-trlsacdef(saclop,1)>postsac
                    sactrl = [sactrl; ...
                        trlsacdef(saclop,1)-presac trlsacdef(saclop,1)+postsac -presac];
                end
%             end
        end
    end
end

cfg = [];
cfg.trl = sactrl;
sacdata = ft_redefinetrial(cfg, data);

cfg = [];
cfg.demean = 'yes';
cfg.detrend = 'yes';
sacdata = ft_preprocessing(cfg,sacdata);

% do the time-locked analysis
cfg = [];
cfg.channel       = sacdata.label;
cfg.covariance    = 'no';
cfg.keeptrials    = 'yes';
cfg.vartrllength  =  2;
timelock = ft_timelockanalysis(cfg,sacdata);

chans = find(strncmp('A_chan',chanLabels,6) | ...
    strncmp('B_chan',chanLabels,6) | ...
    strncmp('C_chan',chanLabels,6));

% calculate the phase using the hilbert transform
for chnlop = 1:length(chans)
    angraw = nan(size(timelock.trial,1),size(timelock.trial,3));
    for trllop = 1:size(timelock.trial,1)
        h = hilbert(squeeze(timelock.trial(trllop,chnlop,:)));
        angraw(trllop,:) = mod(angle(h)+2*pi,2*pi);
    end
    
    % calculate rayleigh statistic
    phasebins = [0:pi/12:2*pi];
    
    % for all trials
    anglp = histc(angraw,phasebins);
    anglp = anglp ./ repmat(sum(anglp),size(anglp,1),1);
    
    unitv = exp(1i*angraw);
    sumv = squeeze(mean(unitv,1));
    
    rayleigh=abs(sumv);
    critval=sqrt(-log(0.01)/size(angraw,1));

end

