% --- Subfunction for plotting and saving ---
function plot_and_save(field, title_str, filename, Yq, Zq, output_folder)
    figure;
    imagesc(unique(Zq(1,:)), flipud(unique(Yq(:,1))), field);
    axis xy
    % colormap('viridis');
    colorbar('Ticks',[min(field(:)) max(field(:))], 'TickLabels',{num2str(min(field(:))), num2str(max(field(:)))});
    xlabel('Z');
    ylabel('Y');
    title(['Velocity Component: ', title_str], 'Interpreter', 'latex');
    axis image
    set(gca, 'YDir', 'normal'); % Ensure correct orientation
    saveas(gcf, fullfile(output_folder, filename));
    close(gcf);
end