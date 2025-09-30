function [pol_apInf,pol_midInf,pol_basInf,pol_apAnt,pol_midAnt,pol_basAnt,x2CH_Av_strain_ML,Names_2CH,pol_basAnt_dia,pol_midAnt_dia,pol_apAnt_dia,pol_basInf_dia,pol_midInf_dia,pol_apInf_dia] = Standard_Strain_project_Interpolation_2CH_30_09_2025(strain_DATA,Segment_type,final_df)
%% create name-vector
Names = final_df.Properties.RowNames;

%% Removing names without 2CH

%creates af vector where 1 is patients with x2CH and 0 is without
for n = 1:length(Names)
x2CH(n) = isfield(strain_DATA.(Names{n,1}),'x2CH');
end

Names_2CH = Names(x2CH == 1);
Names_2CH_table = table('Size',[length(Names_2CH),1],'VariableTypes',{'cell'});
Names_2CH_table.Properties.RowNames = Names_2CH;
A = innerjoin(Names_2CH_table,final_df,'Keys','Row');

%% Times: T0, T_MVO, T1

% Finds T0 (left marker), T_MVO (time betweeen systole and diastole) and T1
% (right marker). These are used to devide the cuves into systole and
% diastole. 
for n = 1:length(Names_2CH)
    t0(n) = strain_DATA.(Names_2CH{n}).x2CH.left_marker;
end

for n = 1:length(Names_2CH)
    if ischar(strain_DATA.(Names_2CH{n}).x2CH.right_marker) == 1;
        strain_DATA.(Names_2CH{n}).x2CH.right_marker = str2double(strain_DATA.(Names_2CH{n}).x2CH.right_marker); 
    end
end

for n = 1:length(Names_2CH)
    t1(n) = strain_DATA.(Names_2CH{n}).x2CH.right_marker;
end

t_MVO = A.MVO';
%% finding wrong time interals
t_false = find(t0 < t_MVO & t_MVO < t1);
t0 = t0(t_false);
t_MVO = t_MVO(t_false);
t1 = t1(t_false);
A = A(t_false,:);
Names_2CH = Names_2CH(t_false);

N = length(Names_2CH);


%% Finds index of T0, MVO and T1 

for n = 1:length(Names_2CH)
    [~,I_t0(n)] = min(abs(strain_DATA.(Names_2CH{n}).x2CH.Time-t0(n)));
end

for n = 1:length(Names_2CH)
    [~,I_MVOtemp(n)] = min(abs(strain_DATA.(Names_2CH{n}).x2CH.Time-t_MVO(n)));
    I_MVO(n) = I_MVOtemp(n) + I_t0(n);
end 

for n = 1:length(Names_2CH)
    [d,I_t1(n)] = min(abs(strain_DATA.(Names_2CH{n}).x2CH.Time-t1(n)));
end 


%% Plots a random example (nr 16 but can be replaced)
cmap = copper(1);
for n = 16:16
    if isfield(strain_DATA.(Names_2CH{n}).x2CH,Segment_type{1}) == 1
    figure()
    plot(strain_DATA.(Names_2CH{n}).x2CH.Time(1:I_t1(n)-I_t0(n)+1),strain_DATA.(Names_2CH{n}).x2CH.(Segment_type{1})(I_t0(n):I_t1(n)),'color', cmap(1,:),'LineWidth',4)
    hold on
    
    AVO = table2array(final_df(Names_2CH{n},'AVO'));
    xline(AVO,':',{'AVO'},'color', cmap(1,:),'LineWidth',4)
    hold on   

    AVC = table2array(final_df(Names_2CH{n},'AVC'));
    xline(AVC,'--',{'AVC'},'color', cmap(1,:),'LineWidth',4)
    hold on 
    
    MVO = table2array(final_df(Names_2CH{n},'MVO'));
    xline(MVO,'-',{'MVO'},'color', cmap(1,:),'LineWidth',4)
    hold on  
    
    MVC = table2array(final_df(Names_2CH{n},'MVC'));
    xline(MVC,'-.',{'MVC'},'color', cmap(1,:),'LineWidth',4)
    hold on  
    pause(0.8)
    end
    %close all
end

%% cycle length and sys/dia length (defined by T0, MVO and T1)
% the length are defined as the mean of the copenhagen cityheart study
% which was 30 for the systole and 31 for the diastole. 

for n = 1:length(Names_2CH)
     len_cyc(n) = I_t1(n)-I_t0(n)+1;
     len_sys(n) = I_MVO(n)-I_t0(n)+1;
     len_dia(n) = I_t1(n)-I_MVO(n)+1;
end

mean_len_sys = 30; 
mean_len_dia = 31; 
mean_len_cyc = 61; 
%% Segments prep
Seg_apAnt_sys = NaN(300,length(Names_2CH));
Seg_basInf_sys = NaN(300,length(Names_2CH));
Seg_basAnt_sys = NaN(300,length(Names_2CH));
Seg_apInf_sys = NaN(300,length(Names_2CH));
Seg_midAnt_sys = NaN(300,length(Names_2CH));
Seg_midInf_sys = NaN(300,length(Names_2CH));

Seg_apAnt_dia = NaN(300,length(Names_2CH));
Seg_basInf_dia = NaN(300,length(Names_2CH));
Seg_basAnt_dia = NaN(300,length(Names_2CH));
Seg_apInf_dia = NaN(300,length(Names_2CH));
Seg_midAnt_dia = NaN(300,length(Names_2CH));
Seg_midInf_dia = NaN(300,length(Names_2CH));

for n = 1:N
    if isfield(strain_DATA.(Names_2CH{n}).x2CH,Segment_type{1}) == 1
        Seg_apInf_sys(1:len_sys(n),n) = strain_DATA.(Names_2CH{n}).x2CH.(Segment_type{1})(I_t0(n):I_MVO(n));
        Seg_apInf_dia(1:len_dia(n),n) = strain_DATA.(Names_2CH{n}).x2CH.(Segment_type{1})(I_MVO(n):I_t1(n));
    end
    
    if isfield(strain_DATA.(Names_2CH{n}).x2CH,Segment_type{2}) == 1
        Seg_midInf_sys(1:len_sys(n),n) = strain_DATA.(Names_2CH{n}).x2CH.(Segment_type{2})(I_t0(n):I_MVO(n));
        Seg_midInf_dia(1:len_dia(n),n) = strain_DATA.(Names_2CH{n}).x2CH.(Segment_type{2})(I_MVO(n):I_t1(n));
    end
    
    if isfield(strain_DATA.(Names_2CH{n}).x2CH,Segment_type{3}) == 1
        Seg_basInf_sys(1:len_sys(n),n) = strain_DATA.(Names_2CH{n}).x2CH.(Segment_type{3})(I_t0(n):I_MVO(n));
        Seg_basInf_dia(1:len_dia(n),n) = strain_DATA.(Names_2CH{n}).x2CH.(Segment_type{3})(I_MVO(n):I_t1(n));
    end
    
    if isfield(strain_DATA.(Names_2CH{n}).x2CH,Segment_type{4}) == 1
        Seg_apAnt_sys(1:len_sys(n),n) = strain_DATA.(Names_2CH{n}).x2CH.(Segment_type{4})(I_t0(n):I_MVO(n));
        Seg_apAnt_dia(1:len_dia(n),n) = strain_DATA.(Names_2CH{n}).x2CH.(Segment_type{4})(I_MVO(n):I_t1(n));
    end
    
    if isfield(strain_DATA.(Names_2CH{n}).x2CH,Segment_type{5}) == 1
        Seg_midAnt_sys(1:len_sys(n),n) = strain_DATA.(Names_2CH{n}).x2CH.(Segment_type{5})(I_t0(n):I_MVO(n));
        Seg_midAnt_dia(1:len_dia(n),n) = strain_DATA.(Names_2CH{n}).x2CH.(Segment_type{5})(I_MVO(n):I_t1(n));
    end
    
     if isfield(strain_DATA.(Names_2CH{n}).x2CH,Segment_type{6}) == 1
        Seg_basAnt_sys(1:len_sys(n),n) = strain_DATA.(Names_2CH{n}).x2CH.(Segment_type{6})(I_t0(n):I_MVO(n));
        Seg_basAnt_dia(1:len_dia(n),n) = strain_DATA.(Names_2CH{n}).x2CH.(Segment_type{6})(I_MVO(n):I_t1(n));
    end
end

%% Interpolate all segments to 61 length

pol_apAnt_sys = NaN(mean_len_sys,length(Names_2CH));
pol_basInf_sys = NaN(mean_len_sys,length(Names_2CH));
pol_basAnt_sys = NaN(mean_len_sys,length(Names_2CH));
pol_apInf_sys = NaN(mean_len_sys,length(Names_2CH));
pol_midAnt_sys = NaN(mean_len_sys,length(Names_2CH));
pol_midInf_sys = NaN(mean_len_sys,length(Names_2CH));

pol_apAnt_dia = NaN(mean_len_dia,length(Names_2CH));
pol_basInf_dia = NaN(mean_len_dia,length(Names_2CH));
pol_basAnt_dia = NaN(mean_len_dia,length(Names_2CH));
pol_apInf_dia = NaN(mean_len_dia,length(Names_2CH));
pol_midAnt_dia = NaN(mean_len_dia,length(Names_2CH));
pol_midInf_dia = NaN(mean_len_dia,length(Names_2CH));

for n = 1:N
x = linspace(0,1,len_sys(n));
v = Seg_apAnt_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_apAnt_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_apAnt_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_apAnt_dia(:,n) = interp1(x,v,xq);

% 2
x = linspace(0,1,len_sys(n));
v = Seg_basInf_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_basInf_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_basInf_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_basInf_dia(:,n) = interp1(x,v,xq);

%3
x = linspace(0,1,len_sys(n));
v = Seg_basAnt_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_basAnt_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_basAnt_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_basAnt_dia(:,n) = interp1(x,v,xq);

%4
x = linspace(0,1,len_sys(n));
v = Seg_apInf_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_apInf_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_apInf_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_apInf_dia(:,n) = interp1(x,v,xq);

%5
x = linspace(0,1,len_sys(n));
v = Seg_midAnt_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_midAnt_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_midAnt_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_midAnt_dia(:,n) = interp1(x,v,xq);

%6
x = linspace(0,1,len_sys(n));
v = Seg_midInf_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_midInf_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_midInf_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_midInf_dia(:,n) = interp1(x,v,xq);

end 

pol_apAnt = [pol_apAnt_sys(1:mean_len_sys-1,:);pol_apAnt_dia];
pol_basInf = [pol_basInf_sys(1:mean_len_sys-1,:);pol_basInf_dia];
pol_basAnt = [pol_basAnt_sys(1:mean_len_sys-1,:);pol_basAnt_dia];
pol_apInf = [pol_apInf_sys(1:mean_len_sys-1,:);pol_apInf_dia];
pol_midAnt = [pol_midAnt_sys(1:mean_len_sys-1,:);pol_midAnt_dia];
pol_midInf = [pol_midInf_sys(1:mean_len_sys-1,:);pol_midInf_dia];

%% mean (average) of all segments
for i = 1:length(Names_2CH)
for n = 1:(mean_len_cyc-1)
Av_strain(n,i) = mean([pol_apAnt(n,i),pol_basInf(n,i),pol_basAnt(n,i),pol_apInf(n,i),pol_midAnt(n,i),pol_midInf(n,i)],'omitnan');
end 
end

x2CH_Av_strain_ML = Av_strain';

%% plot average strain
figure ()

for n = 1:N
plot(linspace(0,1,mean_len_cyc-1),Av_strain(1:mean_len_cyc-1,n),'Color',[0, 0, 0, 0.05])
hold on
end
title('2CH mean of all segments')
hold off

%% organize data
x2CH_Av_strain_ML = array2table(x2CH_Av_strain_ML);
x2CH_Av_strain_ML.Properties.RowNames = Names_2CH;

pol_apAnt = pol_apAnt';
pol_basInf = pol_basInf';
pol_basAnt = pol_basAnt';
pol_apInf = pol_apInf';
pol_midAnt = pol_midAnt';
pol_midInf = pol_midInf';

pol_apAnt = array2table(pol_apAnt);
pol_basInf = array2table(pol_basInf);
pol_basAnt = array2table(pol_basAnt);
pol_apInf = array2table(pol_apInf);
pol_midAnt = array2table(pol_midAnt);
pol_midInf = array2table(pol_midInf);

pol_apAnt.Properties.RowNames = Names_2CH;
pol_basInf.Properties.RowNames = Names_2CH;
pol_basAnt.Properties.RowNames = Names_2CH;
pol_apInf.Properties.RowNames = Names_2CH;
pol_midAnt.Properties.RowNames = Names_2CH;
pol_midInf.Properties.RowNames = Names_2CH;

pol_basAnt_dia = pol_basAnt_dia';
pol_midAnt_dia = pol_midAnt_dia';
pol_apAnt_dia = pol_apAnt_dia';
pol_basInf_dia = pol_basInf_dia';
pol_midInf_dia = pol_midInf_dia';
pol_apInf_dia = pol_apInf_dia';

pol_basAnt_dia = array2table(pol_basAnt_dia);
pol_midAnt_dia = array2table(pol_midAnt_dia);
pol_apAnt_dia = array2table(pol_apAnt_dia);
pol_basInf_dia = array2table(pol_basInf_dia);
pol_midInf_dia = array2table(pol_midInf_dia);
pol_apInf_dia = array2table(pol_apInf_dia);

pol_basAnt_dia.Properties.RowNames = Names_2CH;
pol_midAnt_dia.Properties.RowNames = Names_2CH;
pol_apAnt_dia.Properties.RowNames = Names_2CH;
pol_basInf_dia.Properties.RowNames = Names_2CH;
pol_midInf_dia.Properties.RowNames = Names_2CH;
pol_apInf_dia.Properties.RowNames = Names_2CH;

end
