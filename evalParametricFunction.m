
function y = evalParametricFunction(val_struct, x, par, fail_if_param_out)
    
    [header_id_low, header_id_high] = getBoundingIds(val_struct.headers, par, fail_if_param_out);
    if isnan(header_id_low) && fail_if_param_out
        error('Parameter is out of bound')
    end
        

    vals_low = val_struct.values{header_id_low};
    vals_high = val_struct.values{header_id_high};

    y_low = interp1(vals_low(:,1), vals_low(:,2), x, 'linear', NaN);
    y_high = interp1(vals_high(:,1), vals_high(:,2), x, 'linear', NaN);

    coeff_y = ((par - val_struct.headers(header_id_low)))/(val_struct.headers(header_id_high) - val_struct.headers(header_id_low));
    coeff_y = min([1, max([0, coeff_y])]);
    y = y_low + coeff_y * (y_high - y_low);


end


function [id_low, id_high] = getBoundingIds(vals, par, fail_if_outofbound)

    id_low = NaN;
    id_high = NaN;
    if par<min(vals)
        if fail_if_outofbound
            return
        else
            id_low = 1;
            id_high = 2;
        end
    elseif par>max(vals)
        if fail_if_outofbound
            return
        else
            id_high = numel(vals);
            id_low = id_high - 1;
        end
    else
        id_high = find(vals>=par, 1, 'first');
        id_low = id_high - 1;
    end

end
