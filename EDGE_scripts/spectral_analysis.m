
%%

K = cell(1,num_embryos);

for i = 1:num_embryos

    cellsOI = cells_raw.get_embryoID(i);
    total_frames = size(cellsOI(1).area,1);
    numFFT = 2^nextpow2( total_frames );
    K{i} = zeros(numel(cellsOI),numFFT);

    for j = 1:numel(cellsOI)
        
        m = central_diff( cellsOI(j).myosin_sm, ...
<<<<<<< HEAD
                (1:total_frames)*input(i).dt );
=======
                (1:total_frames)*input_rok(i).dt );
>>>>>>> FETCH_HEAD
%         m = cellsOI(j).myosin_sm;
        
%         keyboard
        if numel(nonans(m)) > 2
            m = interp_and_truncate_nan(m);
            K{i}(j,:) = fft(m,numFFT);
        end
        
    end
end

%%

for i = 1:2
    power = mean( squeeze( abs(K{i}) ),1);
<<<<<<< HEAD
    plot((0:numFFT/2-1)/numFFT/input_gap(i).dt, ...
=======
    semilogy((0:numFFT/2-1)/numFFT/input_rok(i).dt, ...
>>>>>>> FETCH_HEAD
        power(1:numFFT/2)/max(power) , ...
        'Color', C(i,:)); hold on;
end

%%

