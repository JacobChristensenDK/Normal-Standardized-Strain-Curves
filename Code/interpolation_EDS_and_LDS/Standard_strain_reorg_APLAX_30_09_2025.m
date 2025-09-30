function [Seg_basPost,Seg_midPost,Seg_apPost,Seg_apAntSept,Seg_midAntSept,Seg_basAntSept,Time_APLAX] = Standard_strain_reorg_APLAX_30_09_2025(strain_DATA,Segment_type,names)

N = length(names);

%% Index for T0
I_t0 = nan(N,1);
t0 = nan(N,1);
I_t1 = nan(N,1);
t1 = nan(N,1);
len_cyc = nan(N,1);
Time_APLAX = NaN(N, 300);

for n = 1:N
    name = names{n};
    if isfield(strain_DATA.(name),'APLAX')
        t0(n) = strain_DATA.(name).APLAX.left_marker;
        [~,I_t0(n)] = min(abs(strain_DATA.(name).APLAX.Time-t0(n)));
        if ischar(strain_DATA.(name).APLAX.right_marker) == 1
            strain_DATA.(name).APLAX.right_marker = str2double(strain_DATA.(name).APLAX.right_marker); 
        end
        t1(n) = strain_DATA.(name).APLAX.right_marker;
        [~,I_t1(n)] = min(abs(strain_DATA.(name).APLAX.Time-t1(n)));
        len_cyc(n) = I_t1(n)-I_t0(n)+1;
        Time_APLAX(n,1:len_cyc(n)) = strain_DATA.(name).APLAX.Time(1:len_cyc(n)); 

    end
end

%% Segments 
Seg_basPost = NaN(N,300);
Seg_midPost = NaN(N,300);
Seg_apPost = NaN(N,300);
Seg_apAntSept = NaN(N,300);
Seg_midAntSept = NaN(N,300);
Seg_basAntSept = NaN(N,300);

for m = 1:6
    for n = 1:N
        name = names{n};
        if isfield(strain_DATA.(name),'APLAX')
            if isfield(strain_DATA.(name).APLAX,Segment_type{m})
                strain_values_n = strain_DATA.(name).APLAX.(Segment_type{m})(I_t0(n):I_t1(n));
                switch Segment_type{m}
                    case 'basPost'
                        Seg_basPost(n, 1:len_cyc(n)) = strain_values_n;
                    case 'midPost'
                        Seg_midPost(n, 1:len_cyc(n)) = strain_values_n;
                    case 'apPost'
                        Seg_apPost(n, 1:len_cyc(n)) = strain_values_n;
                    case 'apAntSept'
                        Seg_apAntSept(n, 1:len_cyc(n)) = strain_values_n;
                    case 'midAntSept'
                        Seg_midAntSept(n, 1:len_cyc(n)) = strain_values_n;
                    case 'basAntSept'
                        Seg_basAntSept(n, 1:len_cyc(n)) = strain_values_n;
                    otherwise
                        disp('unknown segment type')
                end
            end
        end
    end
end


%%
Seg_basPost = array2table(Seg_basPost);
Seg_midPost = array2table(Seg_midPost);
Seg_apPost = array2table(Seg_apPost);
Seg_apAntSept = array2table(Seg_apAntSept);
Seg_midAntSept = array2table(Seg_midAntSept);
Seg_basAntSept = array2table(Seg_basAntSept);

Seg_basPost.Properties.RowNames = names;
Seg_midPost.Properties.RowNames = names;
Seg_apPost.Properties.RowNames = names;
Seg_apAntSept.Properties.RowNames = names;
Seg_midAntSept.Properties.RowNames = names;
Seg_basAntSept.Properties.RowNames = names;

Time_APLAX = array2table(Time_APLAX);
Time_APLAX.Properties.RowNames = names;

end
