function [ condition_number, c11 ] = c11_condition_number...
    ( params, paramsP, grid, qdeg, mu, c11_min, c11_max, c11_num_interval )
%C11_CONDITION_NUMBER Summary of this function goes here
%   Detailed explanation goes here

increment = (c11_max-c11_min)/(c11_num_interval);
c11 = c11_min:increment:c11_max;
condition_number = zeros(length(c11),1);

for i=1:1:length(c11)
    
    [ params, paramsP, rhs, stifness_matrix] = assemble_stifness_matrix...
        ( params, paramsP, grid, qdeg, mu, c11(i) );
    condition_number(i,1) = cond(full(stifness_matrix));
    
end

figure()
plot(c11,condition_number);
title('condition numer vs c11 graph')