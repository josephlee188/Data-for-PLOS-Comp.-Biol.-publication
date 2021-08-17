function C = pearson_corr(a, b)

ma = nanmean(a);
mb = nanmean(b);

C = nansum((a - ma) .* (b - mb));

C = C / sqrt(nansum((a - ma).^2)) / sqrt(nansum((b - mb).^2));

end