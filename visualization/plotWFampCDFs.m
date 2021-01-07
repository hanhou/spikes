
function f = plotWFampCDFs(pdfs, cdfs, ampBins, depthBins, lfpSurfaceCh)
f = figure; 

depthX = depthBins(1:end-1)+mean(diff(depthBins))/2;
depthX = lfpSurfaceCh * 10 - depthX;   % Use LFP surface chan

ampX = ampBins(1:end-1)+mean(diff(ampBins))/2;

subplot(1,2,1); 
imagesc(ampX, depthX, pdfs)
xlabel('spike amplitude (�V)');
% ylabel('depth on probe (�m)');
ylabel('Distance from pia (LFP surface, �m)')

title('pdf');
% set(gca, 'YDir', 'normal');
set(gca, 'YDir', 'reverse');

hold on; plot(xlim(), [0 0], 'b--', 'linew', 2);
text(min(xlim())+100, -100, sprintf('LFP surface: ch #%g', lfpSurfaceCh), 'color', 'b');

makepretty

subplot(1,2,2); 
imagesc(ampX, depthX, cdfs)
xlabel('spike amplitude (�V)');
% ylabel('depth on probe (�m)');
ylabel('Distance from pia (LFP surface, �m)')

title('inverse cdf');
% set(gca, 'YDir', 'normal');
set(gca, 'YDir', 'reverse');

colorbar
colormap(colormap_greyZero_blackred)
caxis([0 20]);
makepretty

ch = get(f, 'Children');
chTypes = get(ch, 'Type');
cbar = ch(strcmp(chTypes, 'colorbar'));
cbar.Label.String = 'firing rate (sp/s)';

hold on; plot(xlim(), [0 0], 'b--', 'linew', 2);
text(min(xlim())+100, -100, sprintf('LFP surface: ch #%g', lfpSurfaceCh), 'color', 'b');




