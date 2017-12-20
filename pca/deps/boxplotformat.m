function boxplotformat(colorV)
set(findobj(gca,'tag','Upper Adjacent Value'), 'Color', 'w'); 
set(findobj(gca,'tag','Lower Adjacent Value'), 'Color', 'w'); 
set(findobj(gca,'tag','Upper Whisker'), 'LineStyle', '-'); 
set(findobj(gca,'tag','Lower Whisker'), 'LineStyle', '-'); 
set(findobj(gca,'tag','Outliers'), 'Marker', 'o'); 
h = findobj(gca,'tag','Outliers');
if ~isempty(h)
    for i = 1:length(h)
        if nargin > 0
            set(h(i), 'MarkerEdgeColor', colorV(length(h)-i+1));
        else
            set(h(i), 'MarkerEdgeColor', 'r');
        end
    end
end
%set(findobj(gca,'tag','Outliers'), 'MarkerEdgeColor', 'k'); 
end
