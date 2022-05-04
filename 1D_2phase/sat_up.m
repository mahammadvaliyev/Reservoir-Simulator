function Swup = sat_up(p_i_1,p_i,Sw_i_1,Sw_i)
% Relative permeability function for water

if (p_i_1>p_i)
    Swup=Sw_i_1;
else
    Swup=Sw_i;
end
end
