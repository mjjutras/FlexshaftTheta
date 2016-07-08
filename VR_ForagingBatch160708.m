fileMat = {'2014\JN_14_10_24\JN_14_10_24_12_36' 'JN141024002';
    '2014\JN_14_10_27\JN_14_10_27_12_56' 'JN141027003';
    '2014\JN_14_10_27\JN_14_10_27_13_59' 'JN141027005';
    '2014\JN_14_10_27\JN_14_10_27_14_58' 'JN141027007';
    '2014\JN_14_10_29\JN_14_10_29_12_22' 'JN141029002';
    '2014\JN_14_10_29\JN_14_10_29_13_27' 'JN141029004';
    '2014\JN_14_10_30\JN_14_10_30_14_43' 'JN141030002';
    '2014\JN_14_10_30\JN_14_10_30_15_48' 'JN141030004';
    '2014\JN_14_11_03\JN_14_11_03_14_03' 'JN141103002';
    '2014\JN_14_11_03\JN_14_11_03_15_02' 'JN141103004';
    '2014\JN_14_11_03\JN_14_11_03_16_02' 'JN141103006'};

homeDir = 'C:\Users\michael.jutras\Documents\MATLAB';
datDir = 'C:\Data\VR_Blackrock\';
trldatDir = 'R:\Buffalo Lab\Mike\VirtualNavigationProject\MATFiles\trldat';
pandaDir = 'R:\Buffalo Lab\VR Task Data UW\Giuseppe\panda data\';
NSDir = 'R:\Buffalo Lab\Mike\VirtualNavigationProject\MATFiles\NSdat';
pEpiDir = 'R:\Buffalo Lab\Mike\VirtualNavigationProject\MATFiles\pEpisode\Flexshaft_JN2014implant';

UVtAll = {};
for fillop = 1:size(fileMat,1)

    BRnam = fileMat{fillop,2};
    logDir = fileMat{fillop,1};
    [~,logNam] = fileparts(logDir);
    trldatNam = [logNam '_trldat.mat'];
    
    if exist(fullfile(trldatDir,trldatNam),'file')
        load(fullfile(trldatDir,trldatNam))
        disp(['Loaded ' trldatNam])
    else
        
        disp('Running PARtoTRLDAT')
        cd('VirtualNav\PARtoTRLDAT')
        
        sessionDir = fullfile(pandaDir,logDir);
        trldat = PARtoTRLDAT(sessionDir,0);
        movefile(fullfile(sessionDir,trldatNam),trldatDir)
        cd(homeDir)
        
    end
    
    % load Blackrock data
    NS2 = openNSx(fullfile(datDir,[BRnam '.ns2']),'read');
    
    chanLabels = {NS2.ElectrodesInfo.Label}';
    
    LFPind = strncmp('A_chan',chanLabels,6) | ...
        strncmp('B_chan',chanLabels,6) | ...
        strncmp('C_chan',chanLabels,6);
    eyeXind = strncmp('ainp1',chanLabels,5);
    eyeYind = strncmp('ainp2',chanLabels,5);
    
    % pEpisode (only run once, saves results)
    if exist(fullfile(pEpiDir,[BRnam '_NS2_foraging_pEpi_160708.mat']),'file')
        load(fullfile(pEpiDir,[BRnam '_NS2_foraging_pEpi_160708.mat']))
        disp(['Loaded ' BRnam '_NS2_foraging_pEpi_160708.mat'])
    else
     
        disp('Calculating pEpisode')
       
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
        
        save(fullfile(pEpiDir,[BRnam '_NS2_foraging_pEpi_160708.mat']),'UV','-v7.3')
        disp(['Created ' BRnam '_NS2_foraging_pEpi_160708.mat'])
        
    end
    
    % align data streams using eyetracking data
    if exist(fullfile(NSDir,[BRnam '_' logNam '_NS2_foraging_NSdat.mat']),'file')
        load(fullfile(NSDir,[BRnam '_' logNam '_NS2_foraging_NSdat.mat']))
        disp(['Loaded ' BRnam '_' logNam '_NS2_foraging_NSdat.mat'])
    else
        
        disp('Aligning data streams using eye data')
        
        clear trial sampleinfo
        c = 1;
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
            
            c = c+1;
            
        end
        
        % plot eye data trial by trial to ensure accuracy
        alignFigDir = 'R:\Buffalo Lab\Mike\VirtualNavigationProject\Figures\Flexshaft_firstGizImplant2014\VR_Foraging\trialalign_eyetrace';
        if ~isdir(fullfile(alignFigDir,logNam))
            mkdir(alignFigDir,logNam)
        end
        eyexerr = nan(size(trial)); % error between normalized eyetraces per trial
        for trllop = 1:length(trial)
            
            xp = trial{trllop}(end-1,:);
            xb = trial{trllop}(end-3,:);
            
            xpn = xp/(max(xp)-min(xp));
            xbn = xb/(max(xb)-min(xb));
            
            figure
            plot([xpn; xbn]')
            saveas(gcf,fullfile(alignFigDir,logNam,[num2str(trllop) '.png']),'png')
            close
            
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
        
        save(fullfile(NSDir,[BRnam '_' logNam '_NS2_foraging_NSdat.mat']),'data')
        
        disp(['Created ' BRnam '_' logNam '_NS2_foraging_NSdat.mat'])
        
    end
    
    % create UVt: pEpisode across trial
    if exist(fullfile(pEpiDir,[BRnam '_' logNam '_NS2_foraging_pEpiTrial_160708.mat']),'file')
        load(fullfile(pEpiDir,[BRnam '_' logNam '_NS2_foraging_pEpiTrial_160708.mat']))
        disp(['Loaded ' BRnam '_' logNam '_NS2_foraging_pEpiTrial_160708.mat'])
    else
        
        UVt = cell(size(UV));
        for trllop = 1:size(data.sampleinfo,1)
            for chnlop = 1:length(UV)
                dum = UV{chnlop}(:,data.sampleinfo(trllop,1):data.sampleinfo(trllop,2));
                if trllop==1
                    UVt{chnlop}(trllop,:,:) = dum;
                else
                    if size(UVt{chnlop},3)>size(dum,2)
                        UVt{chnlop}(trllop,:,:) = cat(2,dum,nan(size(dum,1),size(UVt{chnlop},3)-size(dum,2)));
                    elseif size(UVt{chnlop},3)<size(dum,2)
                        UVt{chnlop} = cat(3,UVt{chnlop},nan(size(UVt{chnlop},1),size(dum,1),size(dum,2)-size(UVt{chnlop},3)));
                        UVt{chnlop}(trllop,:,:) = dum;
                    end
                end
            end
        end
        
        save(fullfile(pEpiDir,[BRnam '_' logNam '_NS2_foraging_pEpiTrial_160708.mat']),'UVt','-v7.3')
        disp(['Created ' BRnam '_' logNam '_NS2_foraging_pEpiTrial_160708.mat'])
        
    end
    
    for chnlop = 1:length(UVt)
        if fillop==1
            UVtAll{chnlop} = UVt{chnlop}(:,:,1:50000);
        else
            UVtAll{chnlop} = cat(1,UVtAll{chnlop},UVt{chnlop}(:,:,1:50000));
        end
    end
    
end
