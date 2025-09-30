
function [Seg_basInf, Seg_midInf, Seg_apInf,Seg_apAnt,Seg_midAnt,Seg_basAnt,Time_2CH] = Standard_strain_reorg_2CH_30_09_2025(strain_DATA,Segment_type,names)
%% Load pepdata
N = length(names);

%% Index for T0
I_t0 = nan(N,1);
t0 = nan(N,1);
I_t1 = nan(N,1);
t1 = nan(N,1);
len_cyc = nan(N,1);
Time_2CH = NaN(N, 300);


for n = 1:N
    name = names{n};
    if isfield(strain_DATA.(name),'x2CH')
        t0(n) = strain_DATA.(name).x2CH.left_marker;
        [~,I_t0(n)] = min(abs(strain_DATA.(name).x2CH.Time-t0(n)));
        if ischar(strain_DATA.(name).x2CH.right_marker) == 1
            strain_DATA.(name).x2CH.right_marker = str2double(strain_DATA.(name).x2CH.right_marker); 
        end
        t1(n) = strain_DATA.(name).x2CH.right_marker;
        [~,I_t1(n)] = min(abs(strain_DATA.(name).x2CH.Time-t1(n)));
        len_cyc(n) = I_t1(n)-I_t0(n)+1;
        Time_2CH(n,1:len_cyc(n)) = strain_DATA.(name).x2CH.Time(1:len_cyc(n));
    end
end

%% Segments 
Seg_basInf = NaN(N,300);
Seg_midInf = NaN(N,300);
Seg_apInf = NaN(N,300);
Seg_apAnt = NaN(N,300);
Seg_midAnt = NaN(N,300);
Seg_basAnt = NaN(N,300);


for m = 1:6
    for n = 1:N
        name = names{n};
        if isfield(strain_DATA.(name),'x2CH')
            if isfield(strain_DATA.(name).x2CH,Segment_type{m})
                strain_values_n = strain_DATA.(name).x2CH.(Segment_type{m})(I_t0(n):I_t1(n));
                switch Segment_type{m}
                    case 'basInf'
                        Seg_basInf(n, 1:len_cyc(n)) = strain_values_n;
                    case 'midInf'
                        Seg_midInf(n, 1:len_cyc(n)) = strain_values_n;
                    case 'apInf'
                        Seg_apInf(n, 1:len_cyc(n)) = strain_values_n;
                    case 'apAnt'
                        Seg_apAnt(n, 1:len_cyc(n)) = strain_values_n;
                    case 'midAnt'
                        Seg_midAnt(n, 1:len_cyc(n)) = strain_values_n;
                    case 'basAnt'
                        Seg_basAnt(n, 1:len_cyc(n)) = strain_values_n;
                    otherwise
                        disp('unknown segment type')
                end
            end
        end
    end
end

%%
Seg_basInf = array2table(Seg_basInf);
Seg_midInf = array2table(Seg_midInf);
Seg_apInf = array2table(Seg_apInf);
Seg_apAnt = array2table(Seg_apAnt);
Seg_midAnt = array2table(Seg_midAnt);
Seg_basAnt = array2table(Seg_basAnt);

Seg_basInf.Properties.RowNames = names;
Seg_midInf.Properties.RowNames = names;
Seg_apInf.Properties.RowNames = names;
Seg_apAnt.Properties.RowNames = names;
Seg_midAnt.Properties.RowNames = names;
Seg_basAnt.Properties.RowNames = names;

Time_2CH = array2table(Time_2CH);
Time_2CH.Properties.RowNames = names;

end
