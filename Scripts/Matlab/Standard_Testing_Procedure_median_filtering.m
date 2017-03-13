%% Data Gathering

% Initialize
clear
close all
clear figure
clc
startfolder = pwd;

%Load CSV file
% for the mac daddy
cd('/Volumes/data/Projects/BassConnections_2015-2016/Mobile Headset/EOG_Headset/Test Measurements');
% for PCing
%cd('Z:\Projects\BassConnections_2015-2016\Mobile Headset\EOG_Headset\Test Measurements');
trial_number = 1;
eye_link = 1;
%% iterated
for j = 1:trial_number
    [fit(j).data,path] = uigetfile('*.csv', 'Select a CSV file to load');
    cd(path);
    fit(j).signal = csvread([fit(j).data],3,0);
end

for k = 1:eye_link
    [fit(k).data_2,path] = uigetfile('*.csv', 'Select a CSV file to load');
    cd(path);
    fit(k).eyelink = csvread([fit(k).data_2],3,0);
    fit(k).t1 = fit(k).eyelink(:,1);
    fit(k).angle = fit(k).eyelink(:,2);
end
%% Getting Data

for i = 1:trial_number
    %Variable creation
    fit(i).time = fit(i).signal(:,1);
    fit(i).voltage = fit(i).signal(:,2) - mean(fit(i).signal(1:1000,2));
    
    % Raw data plot
    figure(1);
    plot(fit(i).time, fit(i).voltage, 'k-');
    hold on
    % Filtering

    % Low Pass Filtering
    Fs = 10000;
    [fit(i).Pxx,fit(i).f] = periodogram(fit(i).voltage,rectwin(length(fit(i).voltage)),length(fit(i).voltage),Fs);

    % Design the big filter   
    
    Fp = 50;     % Passband frequency in Hz
    Fst = 55;    % Stopband frequency in Hz
    Fp1 = 9;     % Passband frequency in Hz for no movement regions
    Fst1 = 10;    % Stopband frequency in Hz for no movement regions
    Ap = 1;      % Passband ripple in dB
    Ast = 40;    % Stopband attenuation in dB

    LPF = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',Fst1/(Fs/2),'DesignMethod','butter'); 
    fit(i).v = filtfilt(LPF,fit(i).voltage);
    
    [fit(i).PXX,fit(i).F] = periodogram(fit(i).v,rectwin(length(fit(i).v)),length(fit(i).v),Fs);       
    

    % Rap and Filtered Data Plots
    % Rap and Filtered Data Plots
    figure(2)
    subplot(211)
    plot(fit(i).time,fit(i).v)
    hold on
    grid on
    xlabel('t [s]')
    ylabel('Signal [V]')
    title('Filtered Data')
    xlim([0 50])
    subplot(212)
    semilogx(fit(i).F,10*log(fit(i).PXX))
    hold on
    ax=gca;
    ax.XMinorTick = 'on';
    grid on
    xlim([0.1 5000])
    xlabel('Log Frequency (Hz)')
    ylabel('Poper/Frequency (dB/Hz)')
    title('PSD of Filtered Data')
    
    figure(3)
    subplot(211)
    plot(fit(i).time,fit(i).voltage)
    hold on
    grid on
    xlabel('t [s]')
    ylabel('Signal [V]')
    title('Unfiltered Data')
    xlim([0 50])
    subplot(212)
    semilogx(fit(i).f,10*log(fit(i).Pxx))
    hold on
    ax=gca;
    ax.XMinorTick = 'on';
    grid on
    xlim([0.1 5000])
    xlabel('Log Frequency (Hz)')
    ylabel('Poper/Frequency (dB/Hz)')
    title('PSD of Unfiltered Data')
    
    figure(4);
    plot(fit(i).time,fit(i).v)
    grid on 
    xlabel('t [s]')
    ylabel('Signal [V]')
    legendInfo{i} = ['Trial ' num2str(i)];
    legend(legendInfo, 'location', 'best');
    %legend(num2str(i), 'location', 'southpest');
    hold on
end   

%% Getting differences for median filtering (not finished)
for i = 1:trial_number
    fit(i).down_time = downsample(fit(i).time,1);
    fit(i).down_v = downsample(fit(i).v,1);
    fit(i).volt_diff = diff(fit(i).down_v);
    figure(5);
    plot(fit(i).down_time(2:end),fit(i).volt_diff);
    hold on;
    
    averaging_time = 0.5;
    median_window = averaging_time.*(Fs);
    
    fit(i).detect_change = find(fit(i).time>1); 
    eye_movement = 0; %when this is zero, median filtering
    r = 1; % provides matrix indice for start_change
    t = 1; % provides matrix indice for end_change
    
    for j = 1:length(fit(i).volt_diff)
        if (abs(fit(i).volt_diff(j))>(3e-6)) && eye_movement == 0 && (j > fit(i).detect_change(1))
            fit(i).start_change(r,1) = j;
            eye_movement = 1;
            r = r+1;
        elseif (abs(fit(i).volt_diff(j))<(3e-6)) && eye_movement == 1 && (j > fit(i).detect_change(1))
            fit(i).end_change(t,1) = j;
            eye_movement = 0;
            t = t+1;
        end
    end
    start_point = 1;
    end_point = 1;
    for u = 1:length(fit(i).start_change)
        while start_point == 1
            if sign(fit(i).volt_diff(fit(i).start_change(1)))-sign(fit(i).volt_diff(fit(i).start_change(1)-1)) == 0
                fit(i).start_change(u) = fit(i).start_change(u) - 1;
            elseif sign(fit(i).volt_diff(fit(i).start_change(1)))-sign(fit(i).volt_diff(fit(i).start_change(1)-1)) ~= 0
                start_point = 0;
            end 
        end
    end
    for u = 1:length(fit(i).end_change)
        while end_point == 1
            if sign(fit(i).volt_diff(fit(i).end_change(1)))-sign(fit(i).volt_diff(fit(i).end_change(1)-1)) == 0
                fit(i).end_change(u) = fit(i).end_change(u) + 1;
            elseif sign(fit(i).volt_diff(fit(i).end_change(1)))-sign(fit(i).volt_diff(fit(i).end_change(1)-1)) ~= 0
                end_point = 0;
            end 
        end
    end
    fit(i).median_values = medfilt1(fit(i).v,median_window);
    fit(i).new_voltage = fit(i).v - fit(i).median_values;
    for p = 1:length(fit(i).end_change)
        fit(i).change_indices(p).x = [fit(i).start_change(p):1:fit(i).end_change(p)]';
    end
    for p = 1:(length(fit(i).end_change)-1)
        for q = 1:length(fit(i).change_indices(p).x)
            fit(i).new_voltage(fit(i).change_indices(p).x(q)) = fit(i).new_voltage(fit(i).change_indices(p).x(q)-1) + fit(i).volt_diff(fit(i).change_indices(p).x(q)+1);
        end
        for q = (fit(i).change_indices(p).x(end)+1):(fit(i).change_indices(p+1).x(1)-1)
            fit(i).new_voltage(q) = fit(i).new_voltage(fit(i).change_indices(p).x(end)) + fit(i).new_voltage(q);
        end
    end
    figure(6);
    plot(fit(i).time,4000*fit(i).new_voltage, fit(i).t1, fit(i).angle);
    hold on;
end