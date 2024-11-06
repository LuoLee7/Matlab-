function [p_l,z_l] = is_conjugate(p,z)
%   返回开环零极点其中的共轭

mask_p = real(p) == 0 | imag(p) == 0;
p(mask_p) = [];
p_l = p;

mask_z = real(z) == 0 | imag(z) == 0;
z(mask_z) = [];
z_l = z;
end

