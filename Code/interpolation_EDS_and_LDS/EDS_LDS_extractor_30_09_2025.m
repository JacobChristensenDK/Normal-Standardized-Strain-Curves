function [EDS_TO_GLS,EDS,LDS] = EDS_LDS_extractor_30_09_2025(S,Time,final_data,AVC_strain,fpass)

%% Code dexirption
% 1. Strain curves are smoothed with a low-pass filter

% 2. i_EDS is located first by differentiating the curves two times (acceleration).
% i_EDS is defined as the first min peak of the diastolic acceleration curve.

% 3. i_LDS is defined as the last peak of the acceleration curve during diastole

%4. for both i_EDS og i_LDS if there are no peaks they are defined as the last sample of the curves

%%
Time = table2array(Time);
FR = round(1./mean(diff(Time,1,2),2,'omitnan'));

[N, M] = size(S);
end_time = nan(N,1);
DS1 = nan(N,M);
DDS1 = nan(N,M);
Time_dd = nan(N,M);
local_min = nan(N,M);

idx_AVC = round(FR.*AVC_strain./1000); %converts to seconds and find the index
idx_missing_AVC = find(isnan(S(:,1)) < isnan(idx_AVC)); % a few AVCs are missing and therefore they are replaced by the AVC meassured by the TDI method
idx_TDI_AVC = round(FR.*final_data.AVC);
idx_AVC(idx_missing_AVC) = idx_TDI_AVC(idx_missing_AVC);

DDS1_dia = nan(N,M);
i_EDS = nan(N,1);
i_LDS = nan(N,1);
EDS = nan(N,1);
LDS = nan(N,1);

%% s
tic
S_LP = nan(N,M);
for n = 1:N
    if ~isnan(S(n,1)) && ~isnan(Time(n,1))
                end_time(n) = length(Time(n,~isnan(Time(n,:))));
                S_LP(n,1:end_time(n)) = lowpass(S(n,1:end_time(n)),fpass,FR(n)); 
    end
end
toc
%%

for n = 1:N
    if ~isnan(S_LP(n,1)) && ~isnan(Time(n,1))
        DS1(n,:) = gradient(S_LP(n,:));
        DDS1(n,:) = gradient(DS1(n,:));
        dia_len = end_time(n)-idx_AVC(n);
        DDS1_dia(n,1:dia_len) = DDS1(n,idx_AVC(n)+1:end_time(n));
        Time_dd(n,1:end_time(n)-2) = Time(n,1:end_time(n)-2); % this and Timd_d_dia and Time_dd_dia is made because gradient() makes the vectors one shorter.
    end
end

%%
[GLS_LP,i_GLS_LP] = min(S_LP,[],2);
[GLS,i_GLS] = min(S,[],2);

%%

[local_min] = islocalmin(DDS1_dia,2);
for n = 1:N
    if ~isnan(S_LP(n,1)) %& ~isnan(Time(n,1))
    temp2 = min(find(find(islocalmin(DDS1_dia(n,:),2)) + idx_AVC(n) > i_GLS_LP(n)));
    temp3 = find(local_min(n,:));
    temp4 = temp3(temp2);
    if isnan(temp4) == 0
        i_EDS(n) = temp4 + idx_AVC(n);
    else
        i_EDS(n) = end_time(n);
    end
    end
end
%%
[local_max] = islocalmax(DDS1_dia,2);

for n = 1:N
    if ~isnan(S_LP(n,1)) && ~isnan(Time(n,1))
    temp2 = max(find(local_max(n,:))); % finds the last index of max of DDS1_dia
    if isnan(temp2) == 0 & (temp2 + idx_AVC(n)) > i_EDS(n)
        i_LDS(n) = temp2 + idx_AVC(n);
    else
        i_LDS(n) = end_time(n);
    end
    end
end

%%

for n = 1:N
    if ~isnan(GLS(n))
       EDS(n,1) = abs(GLS(n))-abs(S(n,i_EDS(n)));
       LDS(n,1) = abs(S(n,i_LDS(n)));
    end
end

EDS_TO_GLS = EDS./GLS_LP;

%% random sample to visualize
% Visualisation
% figure()
% 
% for n = 200:200
%     if ~isnan(S(n,1))
%         plot(Time(n,:),S(n,:),'b')
%         hold on
%         plot(Time(n,:),S_LP(n,:),'r')
%         hold on
%         %plot(Time(n,idx_AVC(n)),S(n,idx_AVC(n)),'*k')
%         %hold on
%         plot(Time(n,i_GLS(n)),GLS(n),'*k')
%         hold on
%         plot(Time(n,i_EDS(n)),EDS(n),'*b')
%         hold on
%         plot(Time(n,i_LDS(n)),LDS(n),'*r')
%         hold on
%         plot(Time_dd(n,:),DDS1(n,:)*50+10,'k')
% 
%     end
% end 


end
