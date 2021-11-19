function [err] = fit_currentspread(p)


[I, ~, Vq] = create_currentspread(p);
if isnan(Vq)
    Vq = max(max(I(:, :, round(size(I, 3)/2))));
end

err = (Vq(1)-p.thresh).^2;

end