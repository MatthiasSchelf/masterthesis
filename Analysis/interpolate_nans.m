function data_interpolated = interpolate_nans(data)
    data_interpolated = fillmissing(data, 'linear', 2);
end