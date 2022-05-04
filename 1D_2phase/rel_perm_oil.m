function kro = rel_perm_oil(Sw)
% Relative permeability function for oil

if (Sw>=0.25 && Sw<=1)
    kro=1.5*(1-Sw)^2;
elseif (Sw<0.25 && Sw>=0)
    kro=1;
else
    msg = 'Error occurred.';
    error(msg)
end
end

