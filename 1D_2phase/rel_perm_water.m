function krw = rel_perm_water(Sw)
% Relative permeability function for water

if (Sw>=0.25 && Sw<=1)
    krw=(Sw-0.25)^2;
elseif (Sw<0.25 && Sw>=0)
    krw=0;
else
    msg = 'Error occurred.';
    error(msg)
end
end

