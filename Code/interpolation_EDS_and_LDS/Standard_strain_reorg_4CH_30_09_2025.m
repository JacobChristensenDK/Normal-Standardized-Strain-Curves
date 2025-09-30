function [Seg_basSept,Seg_midSept,Seg_apSept,Seg_apLat,Seg_midLat,Seg_basLat,Time_4CH] = Standard_strain_reorg_4CH_30_09_2025(strain_DATA,Segment_type,names)
N = length(names);
%% Index for T0
I_t0 = nan(N,1);
t0 = nan(N,1);
I_t1 = nan(N,1);
t1 = nan(N,1);
len_cyc = nan(N,1);
Time_4CH = NaN(N, 300);

for n = 1:N
    name = names{n};
    if isfield(strain_DATA.(name),'x4CH')
        t0(n) = strain_DATA.(name).x4CH.left_marker;
        [~,I_t0(n)] = min(abs(strain_DATA.(name).x4CH.Time-t0(n)));
        if ischar(strain_DATA.(name).x4CH.right_marker) == 1
            strain_DATA.(name).x4CH.right_marker = str2double(strain_DATA.(name).x4CH.right_marker); 
        end
        t1(n) = strain_DATA.(name).x4CH.right_marker;
        [~,I_t1(n)] = min(abs(strain_DATA.(name).x4CH.Time-t1(n)));
        len_cyc(n) = I_t1(n)-I_t0(n)+1;
        Time_4CH(n,1:len_cyc(n)) = strain_DATA.(name).x4CH.Time(1:len_cyc(n));
    end
end

%% Segments 
Seg_basSept = NaN(N,300);
Seg_midSept = NaN(N,300);
Seg_apSept = NaN(N,300);
Seg_apLat = NaN(N,300);
Seg_midLat = NaN(N,300);
Seg_basLat = NaN(N,300);

for m = 1:6
    for n = 1:N
        name = names{n};
        if isfield(strain_DATA.(name),'x4CH') && isfield(strain_DATA.(name).x4CH,Segment_type{m}) == 1
            strain_values_n = strain_DATA.(name).x4CH.(Segment_type{m})(I_t0(n):I_t1(n));
            switch Segment_type{m}
                case 'basSept'
                    Seg_basSept(n, 1:len_cyc(n)) = strain_values_n;
                case 'midSept'
                    Seg_midSept(n, 1:len_cyc(n)) = strain_values_n;
                case 'apSept'
                    Seg_apSept(n, 1:len_cyc(n)) = strain_values_n;
                case 'apLat'
                    Seg_apLat(n, 1:len_cyc(n)) = strain_values_n;
                case 'midLat'
                    Seg_midLat(n, 1:len_cyc(n)) = strain_values_n;
                case 'basLat'
                    Seg_basLat(n, 1:len_cyc(n)) = strain_values_n;
                otherwise
                    disp('unknown segment type')
            end
        end
    end
end


%%
Seg_basSept = array2table(Seg_basSept);
Seg_midSept = array2table(Seg_midSept);
Seg_apSept = array2table(Seg_apSept);
Seg_apLat = array2table(Seg_apLat);
Seg_midLat = array2table(Seg_midLat);
Seg_basLat = array2table(Seg_basLat);

Seg_basSept.Properties.RowNames = names;
Seg_midSept.Properties.RowNames = names;
Seg_apSept.Properties.RowNames = names;
Seg_apLat.Properties.RowNames = names;
Seg_midLat.Properties.RowNames = names;
Seg_basLat.Properties.RowNames = names;

Time_4CH = array2table(Time_4CH);
Time_4CH.Properties.RowNames = names;

end
