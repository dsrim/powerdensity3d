function colorbarLogscale(fig_handle)

current_labels = fig_handle.Children(1).TickLabels;
K = length(current_labels);
for k = 1:K
    fig_handle.Children(1).TickLabels{k} = ['10^{' current_labels{k} '}' ];
end


end
