function [pol_basPost,pol_midPost,pol_apPost,pol_apAntSept,pol_midAntSept,pol_basAntSept,APLAX_CH_Average_strain,Names_APLAX_CH,pol_basAntSept_dia,pol_midAntSept_dia,pol_apAntSept_dia,pol_apPost_dia,pol_midPost_dia,pol_basPost_dia] = Standard_Strain_project_Interpolation_APLAX_30_09_2025(strain_DATA,Segment_type,final_stata)
%% Info
% for details please refer to the similar function for 2CH: Standard_Strain_project_Interpolation_2CH_30_09_2025
%% create name-vector
Names = final_stata.Properties.RowNames;

%% Removing names without APLAX

%creates af vector where 1 is patients with APLAX and 0 is without
for n = 1:length(Names)
APLAX(n) = isfield(strain_DATA.(Names{n,1}),'APLAX');
end

Names_APLAX_CH = Names(APLAX == 1);
Names_APLAX_CH_table = table('Size',[length(Names_APLAX_CH),1],'VariableTypes',{'cell'});
Names_APLAX_CH_table.Properties.RowNames = Names_APLAX_CH;
A = innerjoin(Names_APLAX_CH_table,final_stata,'Keys','Row');

%% Times: T0, T_MVO, T1
for n = 1:length(Names_APLAX_CH)
    t0(n) = strain_DATA.(Names_APLAX_CH{n}).APLAX.left_marker;
end

for n = 1:length(Names_APLAX_CH)
    if ischar(strain_DATA.(Names_APLAX_CH{n}).APLAX.right_marker) == 1;
        strain_DATA.(Names_APLAX_CH{n}).APLAX.right_marker = str2double(strain_DATA.(Names_APLAX_CH{n}).APLAX.right_marker); 
    end
end

for n = 1:length(Names_APLAX_CH)
    t1(n) = strain_DATA.(Names_APLAX_CH{n}).APLAX.right_marker;
end

t_MVO = A.MVO';
%% finding wrong time interals
t_false = find(t0 < t_MVO & t_MVO < t1);
t0 = t0(t_false);
t_MVO = t_MVO(t_false);
t1 = t1(t_false);
A = A(t_false,:);
Names_APLAX_CH = Names_APLAX_CH(t_false);
N = length(Names_APLAX_CH);

for n = 1:length(Names_APLAX_CH)
    [~,I_t0(n)] = min(abs(strain_DATA.(Names_APLAX_CH{n}).APLAX.Time-t0(n)));
end

for n = 1:length(Names_APLAX_CH)
    [~,I_MVOtemp(n)] = min(abs(strain_DATA.(Names_APLAX_CH{n}).APLAX.Time-t_MVO(n)));
    I_MVO(n) = I_MVOtemp(n) + I_t0(n);
end 

% Finder for hver patient den værdi af tiden som er tættest på t0.
for n = 1:length(Names_APLAX_CH)
    [d,I_t1(n)] = min(abs(strain_DATA.(Names_APLAX_CH{n}).APLAX.Time-t1(n)));
end 

%% plot times intervals
figure
cmap = copper(1);
for n = 16:16
    if isfield(strain_DATA.(Names_APLAX_CH{n}).APLAX,Segment_type{13}) == 1
    figure()
    plot(strain_DATA.(Names_APLAX_CH{n}).APLAX.Time(1:I_t1(n)-I_t0(n)+1),strain_DATA.(Names_APLAX_CH{n}).APLAX.(Segment_type{13})(I_t0(n):I_t1(n)),'color', cmap(1,:),'LineWidth',4)
    hold on
    
    AVO = table2array(final_stata(Names_APLAX_CH{n},'AVO'));
    xline(AVO,':',{'AVO'},'color', cmap(1,:),'LineWidth',4)
    hold on   

    AVC = table2array(final_stata(Names_APLAX_CH{n},'AVC'));
    xline(AVC,'--',{'AVC'},'color', cmap(1,:),'LineWidth',4)
    hold on 
    
    MVO = table2array(final_stata(Names_APLAX_CH{n},'MVO'));
    xline(MVO,'-',{'MVO'},'color', cmap(1,:),'LineWidth',4)
    hold on  
    
    MVC = table2array(final_stata(Names_APLAX_CH{n},'MVC'));
    xline(MVC,'-.',{'MVC'},'color', cmap(1,:),'LineWidth',4)
    hold on  
    pause(0.8)
    end
    %close all
end

%% cycle length and sys/dia length

for n = 1:length(Names_APLAX_CH)
     len_cyc(n) = I_t1(n)-I_t0(n)+1;
     len_sys(n) = I_MVO(n)-I_t0(n)+1;
     len_dia(n) = I_t1(n)-I_MVO(n)+1;
end

mean_len_sys = 30;
mean_len_dia = 31; 
mean_len_cyc = 61; 
%% Segments 
Seg_apAntSept_sys = NaN(300,length(Names_APLAX_CH));
Seg_apPost_sys = NaN(300,length(Names_APLAX_CH));
Seg_basAntSept_sys = NaN(300,length(Names_APLAX_CH));
Seg_basPost_sys = NaN(300,length(Names_APLAX_CH));
Seg_midAntSept_sys = NaN(300,length(Names_APLAX_CH));
Seg_midPost_sys = NaN(300,length(Names_APLAX_CH));

Seg_apAntSept_dia = NaN(300,length(Names_APLAX_CH));
Seg_apPost_dia = NaN(300,length(Names_APLAX_CH));
Seg_basAntSept_dia = NaN(300,length(Names_APLAX_CH));
Seg_basPost_dia = NaN(300,length(Names_APLAX_CH));
Seg_midAntSept_dia = NaN(300,length(Names_APLAX_CH));
Seg_midPost_dia = NaN(300,length(Names_APLAX_CH));

for n = 1:N
    if isfield(strain_DATA.(Names_APLAX_CH{n}).APLAX,Segment_type{13}) == 1
        Seg_basPost_sys(1:len_sys(n),n) = strain_DATA.(Names_APLAX_CH{n}).APLAX.(Segment_type{13})(I_t0(n):I_MVO(n));
        Seg_basPost_dia(1:len_dia(n),n) = strain_DATA.(Names_APLAX_CH{n}).APLAX.(Segment_type{13})(I_MVO(n):I_t1(n));
    end
    
    if isfield(strain_DATA.(Names_APLAX_CH{n}).APLAX,Segment_type{14}) == 1
        Seg_midPost_sys(1:len_sys(n),n) = strain_DATA.(Names_APLAX_CH{n}).APLAX.(Segment_type{14})(I_t0(n):I_MVO(n));
        Seg_midPost_dia(1:len_dia(n),n) = strain_DATA.(Names_APLAX_CH{n}).APLAX.(Segment_type{14})(I_MVO(n):I_t1(n));
    end
    
    if isfield(strain_DATA.(Names_APLAX_CH{n}).APLAX,Segment_type{15}) == 1
        Seg_apPost_sys(1:len_sys(n),n) = strain_DATA.(Names_APLAX_CH{n}).APLAX.(Segment_type{15})(I_t0(n):I_MVO(n));
        Seg_apPost_dia(1:len_dia(n),n) = strain_DATA.(Names_APLAX_CH{n}).APLAX.(Segment_type{15})(I_MVO(n):I_t1(n));
    end
    
    if isfield(strain_DATA.(Names_APLAX_CH{n}).APLAX,Segment_type{16}) == 1
        Seg_apAntSept_sys(1:len_sys(n),n) = strain_DATA.(Names_APLAX_CH{n}).APLAX.(Segment_type{16})(I_t0(n):I_MVO(n));
        Seg_apAntSept_dia(1:len_dia(n),n) = strain_DATA.(Names_APLAX_CH{n}).APLAX.(Segment_type{16})(I_MVO(n):I_t1(n));
    end
    
    if isfield(strain_DATA.(Names_APLAX_CH{n}).APLAX,Segment_type{17}) == 1
        Seg_midAntSept_sys(1:len_sys(n),n) = strain_DATA.(Names_APLAX_CH{n}).APLAX.(Segment_type{17})(I_t0(n):I_MVO(n));
        Seg_midAntSept_dia(1:len_dia(n),n) = strain_DATA.(Names_APLAX_CH{n}).APLAX.(Segment_type{17})(I_MVO(n):I_t1(n));
    end
    
     if isfield(strain_DATA.(Names_APLAX_CH{n}).APLAX,Segment_type{18}) == 1
        Seg_basAntSept_sys(1:len_sys(n),n) = strain_DATA.(Names_APLAX_CH{n}).APLAX.(Segment_type{18})(I_t0(n):I_MVO(n));
        Seg_basAntSept_dia(1:len_dia(n),n) = strain_DATA.(Names_APLAX_CH{n}).APLAX.(Segment_type{18})(I_MVO(n):I_t1(n));
    end
end

%% Interpolate all segments to 61

pol_apAntSept_sys = NaN(mean_len_sys,length(Names_APLAX_CH));
pol_apPost_sys = NaN(mean_len_sys,length(Names_APLAX_CH));
pol_basAntSept_sys = NaN(mean_len_sys,length(Names_APLAX_CH));
pol_basPost_sys = NaN(mean_len_sys,length(Names_APLAX_CH));
pol_midAntSept_sys = NaN(mean_len_sys,length(Names_APLAX_CH));
pol_midPost_sys = NaN(mean_len_sys,length(Names_APLAX_CH));

pol_apAntSept_dia = NaN(mean_len_dia,length(Names_APLAX_CH));
pol_apPost_dia = NaN(mean_len_dia,length(Names_APLAX_CH));
pol_basAntSept_dia = NaN(mean_len_dia,length(Names_APLAX_CH));
pol_basPost_dia = NaN(mean_len_dia,length(Names_APLAX_CH));
pol_midAntSept_dia = NaN(mean_len_dia,length(Names_APLAX_CH));
pol_midPost_dia = NaN(mean_len_dia,length(Names_APLAX_CH));

for n = 1:N

x = linspace(0,1,len_sys(n));
v = Seg_apAntSept_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_apAntSept_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_apAntSept_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_apAntSept_dia(:,n) = interp1(x,v,xq);

% 2
x = linspace(0,1,len_sys(n));
v = Seg_apPost_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_apPost_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_apPost_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_apPost_dia(:,n) = interp1(x,v,xq);

%3
x = linspace(0,1,len_sys(n));
v = Seg_basAntSept_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_basAntSept_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_basAntSept_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_basAntSept_dia(:,n) = interp1(x,v,xq);

%4
x = linspace(0,1,len_sys(n));
v = Seg_basPost_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_basPost_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_basPost_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_basPost_dia(:,n) = interp1(x,v,xq);

%5
x = linspace(0,1,len_sys(n));
v = Seg_midAntSept_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_midAntSept_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_midAntSept_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_midAntSept_dia(:,n) = interp1(x,v,xq);

%6
x = linspace(0,1,len_sys(n));
v = Seg_midPost_sys(1:len_sys(n),n);
xq = linspace(0,1,mean_len_sys);
pol_midPost_sys(:,n) = interp1(x,v,xq);

x = linspace(0,1,len_dia(n));
v = Seg_midPost_dia(1:len_dia(n),n);
xq = linspace(0,1,mean_len_dia);
pol_midPost_dia(:,n) = interp1(x,v,xq);

end 

pol_apAntSept = [pol_apAntSept_sys(1:mean_len_sys-1,:);pol_apAntSept_dia];
pol_apPost = [pol_apPost_sys(1:mean_len_sys-1,:);pol_apPost_dia];
pol_basAntSept = [pol_basAntSept_sys(1:mean_len_sys-1,:);pol_basAntSept_dia];
pol_basPost = [pol_basPost_sys(1:mean_len_sys-1,:);pol_basPost_dia];
pol_midAntSept = [pol_midAntSept_sys(1:mean_len_sys-1,:);pol_midAntSept_dia];
pol_midPost = [pol_midPost_sys(1:mean_len_sys-1,:);pol_midPost_dia];

%% mean (average) of all segments
for i = 1:length(Names_APLAX_CH)
for n = 1:(mean_len_cyc-1)
Av_strain(n,i) = mean([pol_apAntSept(n,i),pol_apPost(n,i),pol_basAntSept(n,i),pol_basPost(n,i),pol_midAntSept(n,i),pol_midPost(n,i)],'omitnan');
end 
end

APLAX_CH_Average_strain = Av_strain';

%% plot average strain
figure

for n = 1:N
plot(linspace(0,1,mean_len_cyc-1),Av_strain(1:mean_len_cyc-1,n),'Color',[0, 0, 0, 0.05])
hold on
end
title('APLAX mean of all segments')
hold off

%%
APLAX_CH_Average_strain = array2table(APLAX_CH_Average_strain);
APLAX_CH_Average_strain.Properties.RowNames = Names_APLAX_CH;

pol_apAntSept = pol_apAntSept';
pol_apPost = pol_apPost';
pol_basAntSept = pol_basAntSept';
pol_basPost = pol_basPost';
pol_midAntSept = pol_midAntSept';
pol_midPost = pol_midPost';

pol_apAntSept = array2table(pol_apAntSept);
pol_apPost = array2table(pol_apPost);
pol_basAntSept = array2table(pol_basAntSept);
pol_basPost = array2table(pol_basPost);
pol_midAntSept = array2table(pol_midAntSept);
pol_midPost = array2table(pol_midPost);

pol_apAntSept.Properties.RowNames = Names_APLAX_CH;
pol_apPost.Properties.RowNames = Names_APLAX_CH;
pol_basAntSept.Properties.RowNames = Names_APLAX_CH;
pol_basPost.Properties.RowNames = Names_APLAX_CH;
pol_midAntSept.Properties.RowNames = Names_APLAX_CH;
pol_midPost.Properties.RowNames = Names_APLAX_CH;

pol_basAntSept_dia = pol_basAntSept_dia';
pol_midAntSept_dia = pol_midAntSept_dia';
pol_apAntSept_dia = pol_apAntSept_dia';
pol_apPost_dia = pol_apPost_dia';
pol_midPost_dia = pol_midPost_dia';
pol_basPost_dia = pol_basPost_dia';

pol_basAntSept_dia = array2table(pol_basAntSept_dia);
pol_midAntSept_dia = array2table(pol_midAntSept_dia);
pol_apAntSept_dia = array2table(pol_apAntSept_dia);
pol_apPost_dia = array2table(pol_apPost_dia);
pol_midPost_dia = array2table(pol_midPost_dia);
pol_basPost_dia = array2table(pol_basPost_dia);

pol_basAntSept_dia.Properties.RowNames = Names_APLAX_CH;
pol_midAntSept_dia.Properties.RowNames = Names_APLAX_CH;
pol_apAntSept_dia.Properties.RowNames = Names_APLAX_CH;
pol_apPost_dia.Properties.RowNames = Names_APLAX_CH;
pol_midPost_dia.Properties.RowNames = Names_APLAX_CH;
pol_basPost_dia.Properties.RowNames = Names_APLAX_CH;

end