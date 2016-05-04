% Takes the state of the cell and builds a snapshot
% [BMP4 FGF CHD OCT4 NANOG SOX2 SOX17 CDX2 BRA]
function [] = pattern_vis(patterned_cell, square_dim)
    
    title_var = {'BMP4'; 'FGF'; 'CHD'; 'OCT4'; 'NANOG'; 'SOX2'; 'SOX17'; 'CDX2'; 'BRA'};
    
    figure
    for i=1:9
    
        subplot(3,3,i);
        A = cellfun(@(v) v(i), patterned_cell(:,:));
        lowestValue = min(A(A(:)>0));
        highestValue = max(A(:));
        imagesc(A);
        cmap = jet(256);
        caxis(gca, [lowestValue-2/256, highestValue]);
        cmap(1,:) = [0,0,0];
        colormap(cmap);
        title(title_var(i));
        colorbar;
    
    end

end