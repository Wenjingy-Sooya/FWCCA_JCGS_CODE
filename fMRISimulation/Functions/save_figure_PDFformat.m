function save_figure_PDFformat(fig, output_dir)

    fig_name = get(fig,'Name');
    safe_name = regexprep(fig_name,'[^\w\-+]','_');

    set(fig,'PaperPositionMode','auto');

    exportgraphics(fig,fullfile(output_dir,[safe_name '.pdf']), ...
        'ContentType','vector');

    fprintf('Saved %s to %s\n',safe_name,output_dir);
end