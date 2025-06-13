function [output] = SCHWEFEL(x1, x2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

output = 418.9829*2 - ( x1*sin(sqrt(abs(x1))) + x2*sin(sqrt(abs(x2))) );
%/180*pi
end