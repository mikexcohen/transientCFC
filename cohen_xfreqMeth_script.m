% Code for analyzing transient cross-frequency coupling EEG data. 
% Note that you must load in the data yourself. The script assumes you have
% one channel of data in a vector (called 'data'). Note that this script
% relies on functions from two toolboxes: (1) 'eegfilt' from the EEGLAB
% toolbox (free to download); (2) Matlab's filter toolbox. It's possible to
% do your own filtering, but make sure you use a good one with zero
% phase-shift. 
%
% Feel free to email Mike Cohen (mikexcohen@gmail.com) with any questions. 
% If you use this script or these methods, please cite the following
% article: Cohen, MX (2008) Assessing transient cross-frequency coupling in
% EEG data. J Neurosci Methods. March 15; 168(2):494-9. PubMed ID: 18061683



%% -- initializations -- %%
% list frequencies in Hz
freqs=20:4:150;

% time decimation (in samples, not time units)
decfactor=10;

% sampling rate in Hz
srate=1000;

% size of time window for analysis. I used 400 ms, but see the paper for
% a discussion of how to choose this time window
analwindow=400;

% initialize phase coherence and other matrices
pc=zeros(length(freqs),length(data)/decfactor);
fluxfreq=pc;
fluxfreqPos=1; fs=9999;
clear i


%% analyze data... loop through frequencies and then through time. 
% Note that it might be a good idea to break up the time series into small
% chunks of a few (I used 4) seconds before running this loop. 
for fi=1:length(freqs)

    % get upper freq power time series
    upperfreq=abs(hilbert(eegfilt(data,srate,freqs(fi)-2.5,freqs(fi)+2.5))).^2;

    p=1;
    % loop through time
    for tx=1:decfactor:length(data)

        % get a little bit of data... 
        tempTS=upperfreq(tx:tx+analwindow);
        last=fs(fluxfreqPos); % remember the last fluxfrequency

        % first, find max frequency component of upper powerTS
        powTSfft=abs(fft(tempTS)); powTSfft=powTSfft(1:round(length(powTSfft)/2));
        fs=linspace(0,srate/2,length(powTSfft));
        nyquist=srate/(analwindow/2); % nyquist of small snippet
        [junk,fluxfreqPos]=max(powTSfft(find(fs>nyquist & fs<40)));
        fs=fs(find(fs>nyquist & fs<40));
        fluxfreq(fi,p)=squeeze(fluxfreq(fi,p)) + fs(fluxfreqPos);

        % get phase times series of flux... timesaver if same as last
        if fs(fluxfreqPos)~=last
            fluxphaseall=angle(hilbert(eegfilt(data,srate,fs(fluxfreqPos)-1.5,fs(fluxfreqPos)+1.5)));
        end
        fluxphase=fluxphaseall(tx:tx+analwindow);


        % calculate PhaseCoherence
        powerTSphase=angle(hilbert(zscore(tempTS)))';
        pc(fi,p)=squeeze(pc(fi,p)) + mean( exp(i*(powerTSphase-fluxphase)) );


        % bootstrap for significance at each point
        f=fft(tempTS);
        for b=1:200 % 200 iterations
            
            % there are several methods for bootstrapping. Three
            % alternatives are presented here, but others could be
            % acceptable as well. Make sure two of the three methods are 
            % commented out, otherwise you'll run all of them! 
            
            % Method 1: create surrogate dataset by shuffling phases while
            % preserving power structure. 
            A=abs(f);
            zphs=cos(angle(f))+i*sin(angle(f));
            phasec=angle(hilbert(real(ifft(A.*zphs(randperm(length(zphs)))))));

            
            % Method 2: break up phase time series once and reorder
            stime=randperm(length(powerTSphase)-50);
            stime=stime(1)+25;
            phasec=[ powerTSphase(stime:end) powerTSphase(1:stime-1) ];
            
            
            % Method 3: Randomize phase time series (fastest method, and
            % perhaps the best null hypothesis)
            phasec=powerTSphase(randperm(length(powerTSphase)));

            
            % compute boot-strapped PC
            bpc(b)=abs(mean(exp(i*(phasec-fluxphase))));
        end
        
        % p-value is the percent of bootstrapped values larger than the
        % observed phase coherence
        pvals(fi,p)=length(find(abs(pc(fi,p))<bpc))/b;
        
        % --  end bootstrapping -- %
        

        % convert real time to decimated time;
        downsamptime(p)=EEG.times(tx);
        p=p+1;

    end % end time loop

    % friendly progress updater
    if mod(fi,3)==0, fprintf('%i ',freqs(fi)); end
end % end frequency loop


