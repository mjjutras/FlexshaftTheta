% get LFPs aligned to image onset, across trials including all sessions,
% and do CSD

fileMat = {'2014\JN_14_10_24' '14_10_24_12_23' 'JN141024001';
    '2014\JN_14_10_27' '14_10_27_12_41' 'JN141027002';
    '2014\JN_14_10_27' '14_10_27_13_45' 'JN141027004';
    '2014\JN_14_10_27' '14_10_27_14_43' 'JN141027006';
    '2014\JN_14_10_29' '14_10_29_12_06' 'JN141029001';
    '2014\JN_14_10_29' '14_10_29_13_14' 'JN141029003';
    '2014\JN_14_10_29' '14_10_29_14_23' 'JN141029005';
    '2014\JN_14_10_30' '14_10_30_14_17' 'JN141030001';
    '2014\JN_14_10_30' '14_10_30_15_35' 'JN141030003';
    '2014\JN_14_10_30' '14_10_30_16_35' 'JN141030005';
    '2014\JN_14_11_03' '14_11_03_13_49' 'JN141103001';
    '2014\JN_14_11_03' '14_11_03_14_47' 'JN141103003';
    '2014\JN_14_11_03' '14_11_03_15_50' 'JN141103005'};

datDir = 'C:\Data\Blackrock_VR\'; % hard drive directory containing Blackrock recording datafiles
pEpiDir = 'R:\Mike\VirtualNavigationProject\MATFiles\pEpisode\Flexshaft_JN2014implant';
NSdir = 'R:\Mike\VirtualNavigationProject\MATFiles\NSdat';
pandaDir = 'R:\VR Task Data UW\Giuseppe\panda data\';

datmat = nan(1,36,1000);
sesmat = nan(size(fileMat,1),36,1000);
d = 1;

for fillop = 1:size(fileMat,1)
    
    logDir = fileMat{fillop,1}; % Panda log subdirectory (within pandaDir)
    calNam = fileMat{fillop,2}; % calibration folder name
    BlackrockNam = fileMat{fillop,3}; % Blackrock recording filename
    
    NS2 = openNSx(fullfile(datDir,[BlackrockNam '.ns2']));
    
    chanLabels = {NS2.ElectrodesInfo.Label}';
    
    LFPind = strncmp('A_chan',chanLabels,6) | ...
        strncmp('B_chan',chanLabels,6) | ...
        strncmp('C_chan',chanLabels,6);
    eyeXind = strncmp('ainp1',chanLabels,5);
    eyeYind = strncmp('ainp2',chanLabels,5);
    
    % align data streams
    if exist(fullfile(NSdir,[BlackrockNam '_' calNam '_NS2_picviewing_NSdat.mat']),'file')
        load(fullfile(NSdir,[BlackrockNam '_' calNam '_NS2_picviewing_NSdat.mat']))
        disp(['Loaded ' BlackrockNam '_' calNam '_NS2_picviewing_NSdat.mat'])
    else
        
        disp('Aligning data streams using eye data')
        
        % ID saccades
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
        
        % load the eye calibration file
        fildir = fullfile(pandaDir,logDir);
        eyecal = ['eye_cal2_' calNam];
        timcal = ['time_cal2_' calNam];
        
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
        
        % get image numbers and presentation order
        
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
        
        %     pairings = NaN(2,36);
        clear pairings
        for img = min(imgnum):max(imgnum) %assumes images don't cross over sets
            ind = find(imgnum == img);
            if length(ind) == 2 %otherwise don't care
                pairings(1,img) = ind(1); %novel
                pairings(2,img) = ind(2); %repeat
            end
        end
        pairings(pairings==0) = nan;
        
        % get eye data and interpolate for alignment with neural signals
        
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
        
        % cross-correlate entire eye signal (after downsampling NS2 signals)
        % to get a rough idea of where the lag is
        % sampling rate in Python log should be 240
        % use this to check: round(1/median(diff(eyetimePyt{1})))
        eyeXdownsamp = resample(double(NS2.Data(eyeXind,:)),240,1000);
        [Xs1, lag] = xcorr(eyeXdownsamp, eye(:,2));
        [~,Is1] = max(Xs1);
        eyedatLag = lag(Is1)*(1000/240); % lag in ms
        firstTrialLag = round(eyedatLag+(eyetimePyt{1}(1)*1000)); % approx. lag of first trial
        
        % align data streams using eyetracking data
        clear trial sampleinfo
        c = 1;
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
                
                c = c+1;
                
            end
        end
        
%         % plot eye data trial by trial to ensure accuracy
%         alignFigDir = 'R:\Buffalo Lab\Mike\VirtualNavigationProject\Figures\Flexshaft_firstGizImplant2014\VR_PicFreeviewing\trialalign_eyetrace';
%         if ~isdir(fullfile(alignFigDir,calNam))
%             mkdir(alignFigDir,calNam)
%         end
%         for trllop = 1:length(trial)
%             
%             xp = trial{trllop}(end-1,:);
%             xb = trial{trllop}(end-3,:);
%             
%             xpn = xp/(max(xp)-min(xp));
%             xbn = xb/(max(xb)-min(xb));
%             
%             figure
%             plot([xpn; xbn]')
%             saveas(gcf,fullfile(alignFigDir,calNam,[num2str(trllop) '.png']),'png')
%             close
%             
%         end
        
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
        
        save(fullfile(NSdir,[BlackrockNam '_' calNam '_NS2_picviewing_NSdat.mat']),'data')
        disp(['Created ' BlackrockNam '_' calNam '_NS2_picviewing_NSdat.mat'])
        
    end
    
    for trllop = 1:length(data.trial)
        datmat(d,:,:) = data.trial{trllop}(1:36,1:1000);
        d=d+1;
    end
    
    cfg = [];
    cfg.vartrllength = 2;
    tl = ft_timelockanalysis(cfg,data);
    sesmat(fillop,:,:) = tl.avg(1:36,1:1000);

end

% % plot ERP for all detrended channels:
datmatDT = datmat-(repmat(nanmean(datmat,3),1,1,size(datmat,3)));
figure;plot(squeeze(nanmean(datmatDT,1))')


%% do CSD
% code taken from CalcCSD_MJJ taken in turn from CalcCSD from Nathan Killian

keepsplines = 1;

bad = [2 5 8 9 10 12];
% bad = [];

CSDmethod = 'iCSD_splines';
% CSDmethod = 'standard';

avgForBad = 1; % should definitely do this instead of leaving them out!

dttrl = 'yes'; % detrend the trials, should probably do this too
% dttrl = 'no';


goodset = 1:12; % array A
% goodset = 13:24; % array B
% goodset = 25:36; % array C
good = setdiff(goodset,bad);
labels = {'A01'; 'A02'; 'A03'; 'A04'; 'A05'; 'A06'; 'A07'; ...
    'A08'; 'A09'; 'A10'; 'A11'; 'A12'; 'B01'; 'B02'; 'B03'; 'B04'; ...
    'B05'; 'B06'; 'B07'; 'B08'; 'B09'; 'B10'; 'B11'; 'B12'; 'C01'; ...
    'C02'; 'C03'; 'C04'; 'C05'; 'C06'; 'C07'; 'C08'; 'C09'; 'C10'; ...
    'C11'; 'C12'};

data= [];
c=1;
for trllop = 1:size(datmat,1)
    if isempty(find(isnan(squeeze(datmat(trllop,1,:))), 1))
        data.trial{c} = squeeze(datmat(trllop,good,:));
        data.time{c} = 0:0.001:0.999;
        c=c+1;
    end
end

data.label = labels(good);

LFPdatatemp = [];

shift = 100;
locstmp = ((0:200:2200)+shift)/1000; % mm

% for missing channels, a dummy variable of all zeros needs to be created
% so that a replacement average can be taken
% this also assumes that these 'bad' channnels are in the spreadsheet
possibleLFPs = labels(goodset);
for k = 1:length(possibleLFPs)
    if ~any(strcmp(possibleLFPs{k},data.label))
        data.label = [data.label; possibleLFPs{k}];
        [data.label,sortind] = sort(data.label);
        for kk = 1:length(data.trial)
            data.trial{kk} = [data.trial{kk}; zeros(1,size(data.trial{kk},2))];
            data.trial{kk} = data.trial{kk}(sortind,:);
        end
    end
end

if strmatch(dttrl,'yes') % remove temporal average from each channel
    for k = 1:length(data.trial)
        data.trial{k} = data.trial{k}-repmat(nanmean(data.trial{k},2),1,size(data.trial{k},2));
    end
end

% Assumes that voltage is in Volts to start with, around 0.4 uA/mm3 is the max CSD you'd expect to
% see after hundreds of trials on the average, plexon is in microvolts, need to do *1e-6
%                 cond = 0.3;               % cortical conductivity, default: 0.3; S/m
%                 cond_top = 0.3;           % conductivity on top of cortex, default: cond
%                 diam = 0.5*1e-3;          % activity diameter, default: 500e-6; mm->m
%                 gauss_sigma = 0.15*1e-3;  % standard deviation of the gaussian for gaussian filtering; mm->m
%                 filter_range = 5*gauss_sigma;
% the CSD can be thought of as:
% CSD (of positive charge) for an external/arbitrary location   =   V''/R
% CSD (of negative charge wrt the outside of the membrane) or of positive charge
%     wrt the internal wrt membrane         = - V''/R
%       this is the 'CSD', so the sign of the LFP is the same as the sign of the CSD
%       meaning you should negate the CSD to put the sinks as being more positive
% basically, a locally negative LFP should mean local excitation because the internal
% part of the cell is becoming positive
%  V'' has -2Vo, so a neg. LFP leads to a pos. CSD (ext. source/internal sink)
% -V'' means neg. LFP -> neg. CSD
% standard method as coded here is defined correctly, so that it needs to be negated for sinks
% (excitation) to be positive
%
% so here more positive values mean excitation
% SET PARAMETERS FOR THE CONDUCTANCE AND THE SMOOTHING:
%---------------------------------------------------------
cond = 0.3; % S/m S = 1 A/V, this is 0.3 mhos/meter = 3.33... ohms/meter
cond_top = 0.3;
diam = 0.034*1e-3; % mm->m
% IMPORTANT NOTE:
% standard CSD assumes infinite disc sources so that you can reduce to 1-dimension
% as disc diameter of the presumed sources approaches infinity, this becomes the standard method
% as the diameter becomes very small (tiny sources), the CSD is then proportional to the LFP!

gauss_sigma = 0.1*1e-3; % mm->m
filter_range = 5*gauss_sigma;
% num_zs = 200; % number of out-parameters for make_cubic_splines
num_zs = 50; % number of out-parameters for make_cubic_splines

if avgForBad
    goodkeep = good;
    good = goodset;
else
    goodkeep = good;
end

el_pos = locstmp(good)*1e-3;%mm->m
[el_pos si] = sort(el_pos);
oldchannelorder = good;
newchannelorder = good(si);
newgoodorder = newchannelorder;
newgoodorder(ismember(newgoodorder,bad)) = nan;

Fcs = F_cubic_spline(el_pos,diam,cond,cond_top);

Nrpt = length(data.trial);
for trli = 1:Nrpt %should do for each trial if want to use the trials later, e.g. for SFC
%     dtmp = data.trial{trli}(rows(si),:);
    dtmp = data.trial{trli}(si,:);
    
    if avgForBad,
        % average potentials for bad channels to create surrogate data there
        % to improve the estimate
        for bi = 1:length(bad)
            bind = find(newchannelorder==bad(bi));
            avgind1 = nearest(newgoodorder,bad(bi)-1);
            avgind2 = nearest(newgoodorder,bad(bi)+1);
            if avgind1~=avgind2
                dtmp(bind,:) = nanmean(dtmp([avgind1 avgind2],:),1);
            else
                dtmp(bind,:) = dtmp(avgind1,:);
            end
            data.trial{trli}(oldchannelorder == bad(bi),:) = dtmp(bind,:);
        end
        
    end
    
    % do the CSD calculation
    [zs,CSD_cs] = make_cubic_splines(el_pos,dtmp,Fcs,num_zs);
    if gauss_sigma~=0 %filter iCSD
        [zs,CSD_cs] = gaussian_filtering(zs,CSD_cs,gauss_sigma,filter_range);
    end;
    unit_scale  = 1e-3 *1e-6; % 1e-6 is to correct for the data being in microvolts, A/m^3 -> muA/mm^3, schroeder uses mV/mm^2, CSD = mV/mm^2 * g, g = 0.3 1/ohm-m
    CSD_cs      = CSD_cs*unit_scale;
    
    if keepsplines
        timetmp = data.time{trli};
        LFPdatatemp.trial{trli} = [];
        LFPdatatemp.trial{trli} = CSD_cs;
        LFPdatatemp = var2field(LFPdatatemp,el_pos,zs,gauss_sigma,filter_range,cond,cond_top,diam,0);
        if trli == 1
            labels2 = cell(length(zs),1);
            for k = 1:length(zs)
                labels2{k,1} = num2str(rsig(zs(k)*1e6,0));%microns
            end
            LFPdatatemp.label = labels2;
        end
    else
        % choose nearest points
        %                     if trli == 1
        %                         newrows = zeros(size(el_pos))';
        %                         for pi = 1:length(el_pos)
        %                             newrows(pi) = nearest(zs,el_pos(pi));
        %                         end
        %                     end
        %do interpolation -better than nearest
        timetmp =   data.time{trli};% this allows for variable trial lengths!!
        [X Y]   =   meshgrid(timetmp,zs);%current interpolated data points
        [XI YI] =   meshgrid(timetmp,el_pos);%points to extract
        CSD_new =   interp2(X,Y,CSD_cs,XI,YI);
        LFPdatatemp.trial{trli} = [];
        %     LFPdatatemp.trial{trli} = CSD_cs(newrows,:);%for nearest
        LFPdatatemp.trial{trli} = CSD_new;
        if trli == 1,LFPdatatemp.label = labels(rows(si));end
    end
    
end

h=[];
for k=1:length(LFPdatatemp.trial)
    h(k,:,:)=LFPdatatemp.trial{k};
end
figure;imagesc(0:0.001:0.999,LFPdatatemp.zs,squeeze(nanmean(h,1))); axis xy
set(gca,'YTick',el_pos)
colormap jet

