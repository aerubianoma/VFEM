function [nm,ne] = format_own(n)
ne = floor(log10(n));  % EXPONENT
nm = n/10^ne;          % MANTISSA