function val_struct = getDataFromWebplotDigitizer(filename)
    data_matrix = readmatrix(filename, 'NumHeaderLines', 2);
    
    opts = detectImportOptions(filename);
    opts.DataLines = [1 1];
    header_matrix = readmatrix(filename, opts);
    
    val_struct.headers = header_matrix(header_matrix==header_matrix);
    
    numAxes = size(data_matrix,2)/2;
    val_struct.values = cell(1, numAxes);
    for ax_sel = 1:numAxes
        tableColSel = (ax_sel-1)*2 + 1;
        tableCols = data_matrix(:, [tableColSel, tableColSel+1]);
        val_struct.values{ax_sel} = tableCols(tableCols(:,1)==tableCols(:,1), :);
    
    end


end
