function [dec_bits] = LDPC_harddecode(enc_bits,H,n_iterations)
    %Hard decoding using the Tanner graph
    [sizeCheck, sizeVar] = size(H);
    rate = sizeVar/sizeCheck;
    n_bits = length(enc_bits);
    dec_bits = zeros(1, n_bits/rate); 
    
    for idx = 1:sizeVar:n_bits
        y = enc_bits(idx:idx+sizeVar-1);
        m = y;
        iteration = 0;
        syndrome = mod(m*H',2)';
        
        %loop while syndrome != zero vector and iteration hasn't reached
        %max iteration. Syndrome not zero means that there is an error
        while any(syndrome) && (iteration < n_iterations)
            iteration = iteration + 1;
            m = tanner(y, syndrome, m, H);
            syndrome = mod(H*m',2);
        end
        
        %remove parity bits
        dec_bits(1+(idx-1)/rate:(idx-1)/rate+sizeCheck) = m(end-sizeCheck+1:end);
    end
end
