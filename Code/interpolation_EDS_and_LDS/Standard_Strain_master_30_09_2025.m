clc
clear all
close all
cd 'code location'
%% Load patient data and metadata
df = readtable('patient_data_location.csv');

%code needed in the cophenhagen city heart study
%idx_over4466 = df.QBUS5NR > 4466;
%df(idx_over4466,:) = [];

df.QBUS5NR = string(df.QBUS5NR); %replace QBUS5NR with patient number
idx_empty=all(cellfun(@isempty,df{:,1}),2);
df(idx_empty,:) = [];
df.Properties.RowNames = strcat('x',df.QBUS5NR);

%% load strain data
% strain_DATA must be a struct organized as in the following way: strain_DATA -->
% patient ID --> 3 fields called x2CH (2 chamber), x4CH (four chamber) and
% APLAX (3 chamber). --> up to 9 dobules with Time, 1-6 segment names,
% left_marker and right_marker.
Traces = 'strain_trace_location.json';
fid = fopen(Traces, 'r');
str = fread(fid, '*char').';
fclose(fid);
strain_DATA = jsondecode(str);

%% Define segment names (can also be writen as colors for older echo pac versions)
%2CH:
Segment_type{1} = 'apInf'; 
Segment_type{2} = 'midInf'; 
Segment_type{3} = 'basInf'; 
Segment_type{4} = 'apAnt'; 
Segment_type{5} = 'midAnt';
Segment_type{6} = 'basAnt';
 
%4CH:
Segment_type{7} = 'basSept';
Segment_type{8} = 'midSept';
Segment_type{9} = 'apSept'; 
Segment_type{10} = 'apLat'; 
Segment_type{11} = 'midLat';
Segment_type{12} = 'basLat';

%APLAX:
Segment_type{13} = 'basPost'; 
Segment_type{14} = 'midPost'; 
Segment_type{15} = 'apPost';
Segment_type{16} = 'apAntSept'; 
Segment_type{17} = 'midAntSept';
Segment_type{18} = 'basAntSept';

%% load event times (AVO)
data.event_times = df(:,{'AVO','AVC','MVO'});

is_event_time_unsorted = true(size(df,1),1);
for n = 1:size(df,1)
    is_event_time_unsorted(n,1) = issorted(table2array(data.event_times(n,:))) == 0;
end

is_event_times_nan = isnan(data.event_times.AVO);

%% Exclude subjects with wrong event times and missing strain cuves
names_1 = df.Properties.RowNames;
names_2 = names_1(~is_event_times_nan & ~is_event_time_unsorted);
names_strain = fieldnames(strain_DATA);
names_3 = intersect(names_2,names_strain);
final_names = table('Size',[length(names_3),1],'VariableTypes',{'cell'});
final_names.Properties.RowNames = names_3;
final_df = innerjoin(final_names,df,'Keys','Row');

%% Interpolation
[pol_apInf,pol_midInf,pol_basInf,pol_apAnt,pol_midAnt,pol_basAnt,x2CH_Average_strain,Names_2CH] = Standard_Strain_project_Interpolation_2CH_30_09_2025(strain_DATA,Segment_type,final_df);
[pol_basSept,pol_midSept,pol_apSept,pol_apLat,pol_midLat,pol_basLat,x4CH_Average_strain,Names_4CH] = Standard_Strain_project_Interpolation_4CH_30_09_2025(strain_DATA,Segment_type,final_df);
[pol_basPost,pol_midPost,pol_apPost,pol_apAntSept,pol_midAntSept,pol_basAntSept,APLAX_CH_Average_strain,Names_APLAX_CH] = Standard_Strain_project_Interpolation_APLAX_30_09_2025(strain_DATA,Segment_type,final_df);

%% combine file with all interpolated segments and df (including outcomes)
All_segments_struct = struct('pol_apAnt',pol_apAnt,'pol_apAntSept',pol_apAntSept,'pol_apInf',pol_apInf,'pol_apLat',pol_apLat,'pol_apPost',pol_apPost,'pol_apSept',pol_apSept,'pol_basAnt',pol_basAnt,'pol_basAntSept',pol_basAntSept,'pol_basInf',pol_basInf,'pol_basLat',pol_basLat,'pol_basPost',pol_basPost,'pol_basSept',pol_basSept,'pol_midAnt',pol_midAnt,'pol_midAntSept',pol_midAntSept,'pol_midInf',pol_midInf,'pol_midLat',pol_midLat,'pol_midPost',pol_midPost,'pol_midSept',pol_midSept);
pol_segment_type = fieldnames(All_segments_struct);

all_segments = All_segments_struct.(pol_segment_type{1});
for n = 2:size(pol_segment_type,1)
all_segments = outerjoin(all_segments,All_segments_struct.(pol_segment_type{n}),'Keys','Row');
end

all_segments_Outcomes = innerjoin(all_segments,df,'Keys','Row');
%% Save all segments of interpolated strain curves as csv file. 
%writetable(all_segments_Outcomes,'location.csv');

%% EDS preparation

Segment_type = [{'basInf'} {'midInf'} {'apInf'} {'apAnt'} {'midAnt'} {'basAnt'} {'basSept'} {'midSept'} {'apSept'} {'apLat'} {'midLat'} {'basLat'} {'basPost'} {'midPost'} {'apPost'} {'apAntSept'} {'midAntSept'} {'basAntSept'}];

% reorganize data
[basInf, midInf, apInf,apAnt,midAnt,basAnt,Time_2CH] = Standard_strain_reorg_2CH_30_09_2025(strain_DATA,Segment_type(1:6),names_3);
[basSept,midSept,apSept,apLat,midLat,basLat,Time_4CH] = Standard_strain_reorg_4CH_30_09_2025(strain_DATA,Segment_type(7:12),names_3);
[basPost,midPost,apPost,apAntSept,midAntSept,basAntSept,Time_APLAX] = Standard_strain_reorg_APLAX_30_09_2025(strain_DATA,Segment_type(13:18),names_3);

all_strain_segments = struct('basInf',basInf,'midInf',midInf,'apInf',apInf,'apAnt',apAnt,'midAnt',midAnt,'basAnt',basAnt,'basSept',basSept,'midSept',midSept,'apSept',apSept,'apLat',apLat,'midLat',midLat,'basLat',basLat,'basPost',basPost,'midPost',midPost,'apPost',apPost,'apAntSept',apAntSept,'midAntSept',midAntSept,'basAntSept',basAntSept);
strain_times  = struct('Time_2CH',Time_2CH,'Time_4CH',Time_4CH,'Time_APLAX',Time_APLAX);
segments = fieldnames(all_strain_segments);
strain_time_names = fieldnames(strain_times);
AVC_strain_names = {'g_2CH_AVC__ms_','g_4CH_AVC__ms_','g_APL_AVC__ms_'};
clear(segments{:});
clear(strain_time_names{:});
final_data = innerjoin(final_names,data.event_times,'Keys','Row');
final_data.Var1 = [];

fpass = 4; % set HZ low pass filter
final_AVC_strain(:,1) = final_df.(AVC_strain_names{1});
final_AVC_strain(:,2) = final_df.(AVC_strain_names{2});
final_AVC_strain(:,3) = mean([final_df.(AVC_strain_names{1}),final_df.(AVC_strain_names{2})],2,'omitnan'); %the aplax chamber does not have preciese AVC-times and are replaced by the mean of 2CH and 4CH avc.

for i = 1:length(segments)
    T = all_strain_segments.(segments{i});
    temp_table = innerjoin(final_names,T,'Keys','Row');
    temp_table.Var1 = [];
    final_all_strain.(segments{i}) = temp_table;
end

for i = 1:length(strain_time_names)
    T = strain_times.(strain_time_names{i});
    temp_table = innerjoin(final_names,T,'Keys','Row');
    temp_table.Var1 = [];
    final_strain_times.(strain_time_names{i}) = temp_table;
end

%% EDS/LDS extraction
% loops through all chambers and segments. Might take some time
for m = 1:3 %set to number of chambers to include
    Time = final_strain_times.(strain_time_names{m});
    AVC_strain = final_AVC_strain(:,m);
        for n = 1:6
            i = (m-1)*6 + n;
            S = table2array(final_all_strain.(segments{i}));
            [EDS_TO_GLS(:,i),EDS(:,i),LDS(:,i)] = EDS_LDS_extractor_30_09_2025(S,Time,final_data,AVC_strain,fpass);
        end
end

%% save EDS and LDS
EDS_mean = mean(EDS(:,1:18),2,'omitnan');
LDS_mean = mean(LDS(:,1:18),2,'omitnan');

EDS_LDS_mean = [final_names table(EDS_mean,LDS_mean,'VariableNames',{'EDS','LDS'})];
EDS_LDS_mean.Var1 = final_names.Properties.RowNames;
EDS_LDS_mean = renamevars(EDS_LDS_mean, 'Var1', 'ID');

EDS_final = [final_names array2table(EDS,'VariableNames',segments)];
LDS_final = [final_names array2table(LDS,'VariableNames',segments)];

%save as CSV
%writetable(EDS_final,'location.csv');
%writetable(LDS_final,'location.csv');