
%%

K = cell(1,num_embryos);

for i = 1:num_embryos

    cellsOI = cells_raw.get_embryoID(i);
    total_frames = size(cellsOI(1).area,1);
    numFFT = 2^nextpow2( total_frames );

    for j = 1:numel(cellsOI)
        
        m = central_diff( cellsOI(j).myosin_sm, ...
                (1:total_frames)*input(i).dt );
%         m = cellsOI(j).myosin_sm;
        
%         keyboard
        if numel(nonans(m)) > 2
            m = interp_and_truncate_nan(m);
            K{i}(j,:) = fft(m,numFFT);
        end
        
    end
end

%%

for i = 6:10
    power = mean( squeeze( abs(K{i}) ),1);
    semilogy((0:numFFT/2-1)/numFFT/input(i).dt, ...
        power(1:numFFT/2)/max(power) , ...
        'Color', C(i,:)); hold on;
end

%%

