
%%

K = cell(1,num_embryos);

for i = 1:num_embryos
    cellsOI = cells_raw.get_embryoID(i);
    n = 2^nextpow2( size(cellsOI(1).area,1) );
    for j = 1:numel(cellsOI)
        
        m = central_diff( cellsOI(j).myosin_sm );
        
        if numel(nonans(m)) > 2
            m = interp_and_truncate_nan(m);
            K{i}(i,:) = fft(m,n);
        end
    end
end

%%

for i = 1:2
    power = squeeze( abs(K{i}) );
    semilogy((0:n/2-1)/n/input_gap(3).dt, ...
        power(:,1:n/2) , ...
        'Color', C(i,:)); hold on;
end

%%

