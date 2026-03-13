function salva_fig(nome, folder_name, SALVATAGGIO)

    if SALVATAGGIO == 1
        saveas(gcf, fullfile(folder_name, nome));
    end

end