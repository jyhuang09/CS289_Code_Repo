% Takes values and mimics Warmflash et. al's initial states
%
function [return_arr] = initial_cell(x_cell, y_cell, square_dim)

    % For all cells near the periphery. Have an increased probability state of OCT4 (warmflash) 
    % If OCT4 is present in higher doses, SOX2 is suppressed: Oct4 priority
    % OCT4 Quantitative dynamic range from .7 to .3, 
    % NANOG Quantitative dynamic range from .4 to .1,
    % SOX2 Quantitative dynamic range from .5 to .2,
    % Overall distinction increases for larger colonies.
    % Slow creep for averaged intensity for 1000um colony
    % [BMP4 FGF CHD OCT4 NANOG SOX2 SOX17 CDX2 BRA]

    center_coord = square_dim / 2;
    radial_distance = sqrt((center_coord-x_cell)^2 + (center_coord-y_cell)^2);
    rad_norm = radial_distance/center_coord;
    
    oct4_sto = (rand-.5)/6;
    nanog_sto = (rand-.5)/6;
    sox2_sto = (rand-.5)/6;
    chd_sto = (rand-.5)/6;
    cdx2_sto = (rand-.5)/8;
    
    oct4_ini = .7;
    nanog_ini = .61;
    sox2_ini = .59;
    chd_ini = .75;
    cdx2_ini = .1;
    
    bmp4 = .7 + (rand-.5)/10;
    fgf = .5 + (rand-.5)/10;
    bra = .5 + (rand-.5)/10;
    
    oct4 = rad_norm*(.95-oct4_ini)+oct4_sto + oct4_ini;
    nanog = rad_norm*(.95-nanog_ini)+nanog_sto + nanog_ini;
    sox2 = rad_norm*(.95-sox2_ini)+sox2_sto + sox2_ini; - bmp4/5;  % Suppression of Sox2
    chd =  chd_ini - rad_norm*(.95-chd_ini)+chd_sto;
    cdx2 = rad_norm*(.4-cdx2_ini)+cdx2_sto + cdx2_ini;
    sox17 = 1 - sox2 + (rand-.5)/4;
    %End is just boolean with True Value i.e. not edge case
    return_arr = double([bmp4 fgf chd oct4 nanog sox2 sox17 cdx2 bra, 1]);

end