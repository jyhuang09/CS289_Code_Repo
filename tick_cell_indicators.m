% Updates the state of every non-boundary cell.
% All boundary cells are of 

function [return_arr] = tick_cell_indicators(patterned_cell, square_dim)

% Cell_arr organization is as follows starting from (1):
% [BMP4 FGF CHD OCT4 NANOG SOX2 SOX17 CDX2 BRA]

% Diffusion
% BMP4 has diffusion with a known slow rate 
% FGF signaling with slower rate ? 
% CHD diffusion of inhibitor since is mainly on the outside has fast
% All TFs should have just nominal diffusion
    vis = [.2, .1, 10, .01, .01, .01, .01, .01, .01];
    diff_arr = cell(1,9);
    for i = 1:9 
        diff_arr{i} = dirich_diff(patterned_cell, i,square_dim , vis(i));
    end
% Production/Suppresion
% BMP4 should always have a strong pull towards mid-point but inhibited
    bmp4_midpoint = .8;
    diff_bmp4 = diff_arr{1};
    diff_fgf = diff_arr{2};
    diff_chd = diff_arr{3};
    radius = sqrt(square_dim/2);
    center_coord = square_dim / 2;
  
    bmp4 = .5.*(bmp4_midpoint-diff_bmp4) - .1*(diff_chd) + diff_bmp4 ;


% FGF production should be regulated by  feedback loop

    fgf_midpoint = .7;
    fgf = .1*(fgf_midpoint - diff_fgf) + diff_fgf;

% CHD production is delayed CHD response directly impacted by FGF

    chd = .05*(fgf_midpoint-diff_fgf) + diff_chd;

% Empirical data suggests that SOX2 is significantly inhibited by the
% presence of BMP4 as the overalls of the SMAD2 with SOX2 are reflective
% BRA exists most highly at the interface of SOX2 and BMP4
% FGF impacts nanog
% Depletion of Oct4 mRNA by 50% is sufficient
% to result in the formation of trophectoderm cells, while
% Oct4 overexpression by 50% will promote mesodermal
% and endodermal differentiation

% CDX2 
% Traditional Class of TFs
    diff_oct4 = diff_arr{4};
    diff_nanog = diff_arr{5};
    diff_sox2 = diff_arr{6};

% Highly impacted then by BMP4 (currently runaway)
    oct4 = diff_oct4 - bmp4*.05 - .1*diff_sox2 + .5*(.7-diff_oct4);
    sox2 = diff_sox2 - bmp4*.1 + .5*(.6-diff_sox2) - oct4*.1;
    nanog = -.2*sox2 + diff_nanog + .2*(.59-diff_nanog);
    

% Other class sox17, BRA, CDX2
    sox17 = diff_arr{7} - sox2 *.1 + bmp4*.1 + .3*(.4-diff_arr{7});
    cdx2 = diff_arr{8} - sox2*.05 + .2*(.25-diff_arr{8});
    bra = diff_arr{9} - bmp4*.1 - sox2*.2  + .4*(.5-diff_arr{9});
% All increases have diminshing returns toward 1 (will implement soon)

    temp_arr = cat(3, bmp4, fgf, chd, oct4, nanog, sox2, sox17, cdx2, bra);
    
    for x=1:9
        for i=1:square_dim
            for j=1:square_dim
                if (sqrt((i-center_coord)^2 + (j-center_coord)^2)<(radius-.5)^2)
                else
                    temp_arr(i,j,x)=0;
                end
            end
        end
    end
    
% Collapse all ticked arrays into patterned_cell format; 
    
    temp_arr(:,:,:) = (temp_arr(:,:,:) > 0).*temp_arr(:,:,:);
    
    return_arr = cell(square_dim);
    
    for i = 1:square_dim
        for j = 1:square_dim
            return_arr{i,j} = reshape(temp_arr(i,j,:), 1, 9);
        end
    end

end