function [pol_basSept,pol_midSept,pol_apSept,pol_apLat,pol_midLat,pol_basLat,FourCH_Av_strain_ML,Names_FourCH,pol_basLat_dia,pol_midLat_dia,pol_apLat_dia,pol_apSept_dia,pol_midSept_dia,pol_basSept_dia] = Standard_Strain_project_Interpolation_4CH_30_09_2025(strain_DATA,Segment_type,final_df)
%% info
% for details please refer to the similar function for 2CH: Standard_Strain_project_Interpolation_2CH_30_09_2025
%% create name-vector
Names = final_df.Properties.RowNames;

%% Removing names without 2CH

%creates af vector where 1 is patients with x4CH and 0 is without
for n = 1:length(Names)
x4CH(n) = isfield(strain_DATA.(Names{n,1}),'x4CH');
end

Names_FourCH = Names(x4CH == 1);
Names_FourCH_table = table('Size',[length(Names_FourCH),1],'VariableTypes',{'cell'});
Names_FourCH_table.Properties.RowNames = Names_FourCH;
A = innerjoin(Names_FourCH_table,final_df,'Keys','Row');

%% Times: T0, T_MVO, T1
for n = 1:length(Names_FourCH)
    t0(n) = strain_DATA.(Names_FourCH{n}).x4CH.left_marker;
end

for n = 1:length(Names_FourCH)
    if ischar(strain_DATA.(Names_FourCH{n}).x4CH.right_marker) == 1;
        strain_DATA.(Names_FourCH{n}).x4CH.right_marker = str2double(strain_DATA.(Names_FourCH{n}).x4CH.right_marker); 
    end
end

for n = 1:length(Names_FourCH)
    t1(n) = strain_DATA.(Names_FourCH{n}).x4CH.right_marker;
end

t_MVO = A.MVO';
%% excluding subjects with wrong time interals
t_false = find(t0 < t_MVO & t_MVO < t1);
t0 = t0(t_false);
t_MVO = t_MVO(t_false);
t1 = t1(t_false);
A = A(t_false,:);
Names_FourCH = Names_FourCH(t_false);

N = length(Names_FourCH);

for n = 1:length(Names_FourCH)
    [~,I_t0(n)] = min(abs(strain_DATA.(Names_FourCH{n}).x4CH.Time-t0(n)));
end

for n = 1:length(Names_FourCH)
    [~,I_MVOtemp(n)] = min(abs(strain_DATA.(Names_FourCH{n}).x4CH.Time-t_MVO(n)));
    I_MVO(n) = I_MVOtemp(n) + I_t0(n);
end 

for n = 1:length(Names_FourCH)
    [d,I_t1(n)] = min(abs(strain_DATA.(Names_FourCH{n}).x4CH.Time-t1(n)));
end 
%% plot times intervals for random example
figure
cmap = copper(1);
for n = 16:16
    if isfield(strain_DATA.(Names_FourCH{n}).x4CH,Segment_type{7}) == 1
    figure()
    plot(strain_DATA.(Names_FourCH{n}).x4CH.Time(1:I_t1(n)-I_t0(n)+1),strain_DATA.(Names_FourCH{n}).x4CH.(Segment_type{7})(I_t0(n):I_t1(n)),'color', cmap(1,:),'LineWidth',4)
    hold on
    
    AVO = table2array(final_df(Names_FourCH{n},'AVO'));
    xline(AVO,':',{'AVO'},'color', cmap(1,:),'LineWidth',4)
    hold on   

    AVC = table2array(final_df(Names_FourCH{n},'AVC'));
    xline(AVC,'--',{'AVC'},'color', cmap(1,:),'LineWidth',4)
    hold on 
    
    MVO = table2array(final_df(Names_FourCH{n},'MVO'));
    xline(MVO,'-',{'MVO'},'color', cmap(1,:),'LineWidth',4)
    hold on  
    
    MVC = table2array(final_df(Names_FourCH{n},'MVC'));
    xline(MVC,'-.',{'MVC'},'color', cmap(1,:),'LineWidth',4)
    hold on  
    pause(0.8)
    end
    %close all
end

%% finding systole/diastole length

for n = 1:length(Names_FourCH)
     len_cyc(n) = I_t1(n)-I_t0(n)+1;
     len_sys(n) = I_MVO(n)-I_t0(n)+1;
     len_dia(n) = I_t1(n)-I_MVO(n)+1;
end

mean_len_sys = 30; 
mean_len_dia = 31; 
mean_len_cyc = 61; 
%% Segments 
Seg_apLat_sys = NaN(300,length(Names_FourCH));
Seg_apSept_sys = NaN(300,length(Names_FourCH));
Seg_basLat_sys = NaN(300,length(Names_FourCH));
Seg_basSept_sys = NaN(300,length(Names_FourCH));
Seg_midLat_sys = NaN(300,length(Names_FourCH));
Seg_midSept_sys = NaN(300,length(Names_FourCH));

Seg_apLat_dia = NaN(300,length(Names_FourCH));
Seg_apSept_dia = NaN(300,length(Names_FourCH));
Seg_basLat_dia = NaN(300,length(Names_FourCH));
Seg_basSept_dia = NaN(300,length(Names_FourCH));
Seg_midLat_dia = NaN(300,length(Names_FourCH));
Seg_midSept_dia = NaN(300,length(Names_FourCH));

for n = 1:N
    if isfield(strain_DATA.(Names_FourCH{n}).x4CH,Segment_type{7}) == 1
        Seg_basSept_sys(1:len_sys(n),n) = strain_DATA.(Names_FourCH{n}).x4CH.(Segment_type{7})(I_t0(n):I_MVO(n));
        Seg_basSept_dia(1:len_dia(n),n) = strain_DATA.(Names_FourCH{n}).x4CH.(Segment_type{7})(I_MVO(n):I_t1(n));
    end
    
    if isfield(strain_DATA.(Names_FourCH{n}).x4CH,Segment_type{8}) == 1
        Seg_midSept_sys(1:len_sys(n),n) = strain_DATA.(Names_FourCH{n}).x4CH.(Segment_type{8})(I_t0(n):I_MVO(n));
        Seg_midSept_dia(1:len_dia(n),n) = strain_DATA.(Names_FourCH{n}).x4CH.(Segment_type{8})(I_MVO(n):I_t1(n));
    end
    
    if isfield(strain_DATA.(Names_FourCH{n}).x4CH,Segment_type{9}) == 1
        Seg_apSept_sys(1:len_sys(n),n) = strain_DATA.(Names_FourCH{n}).x4CH.(Segment_type{9})(I_t0(n):I_MVO(n));
        Seg_apSept_dia(1:len_dia(n),n) = strain_DATA.(Names_FourCH{n}).x4CH.(Segment_type{9})(I_MVO(n):I_t1(n));
    end
    
    if isfield(strain_DATA.(Names_FourCH{n}).x4CH,Segment_type{10}) == 1
        Seg_apLat_sys(1:len_sys(n),n) = strain_DATA.(Names_FourCH{n}).x4CH.(Segment_type{10})(I_t0(n):I_MVO(n));
        Seg_apLat_dia(1:len_dia(n),n) = strain_DATA.(Names_FourCH{n}).x4CH.(Segment_type{10})(I_MVO(n):I_t1(n));
    end
    
    if isfield(strain_DATA.(Names_FourCH{n}).x4CH,Segment_type{11}) == 1
        Seg_midLat_sys(1:len_sys(n),n) = strain_DATA.(Names_FourCH{n}).x4CH.(Segment_type{11})(I_t0(n):I_MVO(n));
        Seg_midLat_dia(1:len_dia(n),n) = strain_DATA.(Names_FourCH{n}).x4CH.(Segment_type{11})(I_MVO(n):I_t1(n));
    end
    
     if isfield(strain_DATA.(Names_FourCH{n}).x4CH,Segment_type{12}) == 1
        Seg_basLat_sys(1:len_sys(n),n) = strain_DATA.(Names_FourCH{n}).x4CH.(Segment_type{12})(I_t0(n):I_MVO(n));
        Seg_basLat_dia(1:len_dia(n),n) = strain_DATA.(Names_FourCH{n}).x4CH.(Segment_type{12})(I_MVO(n):I_t1(n));
    end
end

%% Interpolate all segments to 61

pol_apLat_sys = NaN(mean_len_sys,length(Names_FourCH));
pol_apSept_sys = NaN(mean_len_sys,length(Names_FourCH));
pol_basLat_sys = NaN(mean_len_sys,length(Names_FourCH));
pol_basSept_sys = NaN(mean_len_sys,length(Names_FourCH));
pol_midLat_sys = NaN(mean_len_sys,length(Names_FourCH));
pol_midSept_sys = NaN(mean_len_sys,length(Names_FourCH));

pol_apLat_dia = NaN(mean_len_dia,length(Names_FourCH));
pol_apSept_dia = NaN(mean_len_dia,length(Names_FourCH));
pol_basLat_dia = NaN(mean_len_dia,length(Names_FourCH));
pol_basSept_dia = NaN(mean_len_dia,length(Names_FourCH));
pol_midLat_dia = NaN(mean_len_dia,length(Names_FourCH));
pol_midSept_dia = NaN(mean_len_dia,length(Names_FourCH));

for n = 1:N

x = linspace(0,1,len_sys(n));
v = Seg_apLat_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_apLat_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_apLat_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_apLat_dia(:,n) = interp1(x,v,xq);

% 2
x = linspace(0,1,len_sys(n));
v = Seg_apSept_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_apSept_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_apSept_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_apSept_dia(:,n) = interp1(x,v,xq);

%3
x = linspace(0,1,len_sys(n));
v = Seg_basLat_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_basLat_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_basLat_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_basLat_dia(:,n) = interp1(x,v,xq);

%4
x = linspace(0,1,len_sys(n));
v = Seg_basSept_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_basSept_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_basSept_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_basSept_dia(:,n) = interp1(x,v,xq);

%5
x = linspace(0,1,len_sys(n));
v = Seg_midLat_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_midLat_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_midLat_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_midLat_dia(:,n) = interp1(x,v,xq);

%6
x = linspace(0,1,len_sys(n));
v = Seg_midSept_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_midSept_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_midSept_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_midSept_dia(:,n) = interp1(x,v,xq);

end 

pol_apLat = [pol_apLat_sys(1:mean_len_sys-1,:);pol_apLat_dia];
pol_apSept = [pol_apSept_sys(1:mean_len_sys-1,:);pol_apSept_dia];
pol_basLat = [pol_basLat_sys(1:mean_len_sys-1,:);pol_basLat_dia];
pol_basSept = [pol_basSept_sys(1:mean_len_sys-1,:);pol_basSept_dia];
pol_midLat = [pol_midLat_sys(1:mean_len_sys-1,:);pol_midLat_dia];
pol_midSept = [pol_midSept_sys(1:mean_len_sys-1,:);pol_midSept_dia];

%% mean (average) of all segments
for i = 1:length(Names_FourCH)
for n = 1:(mean_len_cyc-1)
Av_strain(n,i) = mean([pol_apLat(n,i),pol_apSept(n,i),pol_basLat(n,i),pol_basSept(n,i),pol_midLat(n,i),pol_midSept(n,i)],'omitnan');
end 
end

FourCH_Av_strain_ML = Av_strain';

%% plot average strain
figure

for n = 1:N
plot(linspace(0,1,mean_len_cyc-1),Av_strain(1:mean_len_cyc-1,n),'Color',[0, 0, 0, 0.05])
hold on
end
title('4CH mean of all segments')
hold off

%%
FourCH_Av_strain_ML = array2table(FourCH_Av_strain_ML);
FourCH_Av_strain_ML.Properties.RowNames = Names_FourCH;

pol_apLat = pol_apLat';
pol_apSept = pol_apSept';
pol_basLat = pol_basLat';
pol_basSept = pol_basSept';
pol_midLat = pol_midLat';
pol_midSept = pol_midSept';

pol_apLat = array2table(pol_apLat);
pol_apSept = array2table(pol_apSept);
pol_basLat = array2table(pol_basLat);
pol_basSept = array2table(pol_basSept);
pol_midLat = array2table(pol_midLat);
pol_midSept = array2table(pol_midSept);

pol_apLat.Properties.RowNames = Names_FourCH;
pol_apSept.Properties.RowNames = Names_FourCH;
pol_basLat.Properties.RowNames = Names_FourCH;
pol_basSept.Properties.RowNames = Names_FourCH;
pol_midLat.Properties.RowNames = Names_FourCH;
pol_midSept.Properties.RowNames = Names_FourCH;

pol_basLat_dia = pol_basLat_dia';
pol_midLat_dia = pol_midLat_dia';
pol_apLat_dia = pol_apLat_dia';
pol_apSept_dia = pol_apSept_dia';
pol_midSept_dia = pol_midSept_dia';
pol_basSept_dia = pol_basSept_dia';

pol_basLat_dia = array2table(pol_basLat_dia);
pol_midLat_dia = array2table(pol_midLat_dia);
pol_apLat_dia = array2table(pol_apLat_dia);
pol_apSept_dia = array2table(pol_apSept_dia);
pol_midSept_dia = array2table(pol_midSept_dia);
pol_basSept_dia = array2table(pol_basSept_dia);

pol_basLat_dia.Properties.RowNames = Names_FourCH;
pol_midLat_dia.Properties.RowNames = Names_FourCH;
pol_apLat_dia.Properties.RowNames = Names_FourCH;
pol_apSept_dia.Properties.RowNames = Names_FourCH;
pol_midSept_dia.Properties.RowNames = Names_FourCH;
pol_basSept_dia.Properties.RowNames = Names_FourCH;

end
