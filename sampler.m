function [outp] = sampler(inp, M, direction)
 % UPSAMPLING IS PUTTING ZEROS INBETWEEN SYMBOLS
 % DOWNSAMPLING IS REMOVING THE ZEROS
out = [];

switch direction
    
    case 'up'
        for i=1:length(inp)
                out = [out inp(i) zeros(1,M)];
        end

    case 'down'

         out = inp(1:M+1:end);
end

outp = out;
end

