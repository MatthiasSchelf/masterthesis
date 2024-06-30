% Function to visualize each individual Interaction Information matrix 
function visualize_II(II_data, participant, data_type)
    figure;
    imagesc(II_data);
    colorbar;
    title(sprintf('II %s for %s', data_type, participant));
    colormap('jet'); 
    clim([-max(abs(II_data(:))) max(abs(II_data(:)))]); 
    axis([1 size(II_data, 2) 1 size(II_data, 1)]);
end