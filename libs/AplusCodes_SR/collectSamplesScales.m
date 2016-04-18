function [lores hires V_pca] = collectSamplesScales(conf, ohires, numscales, scalefactor)

lores = [];
hires = [];

for scale = 1:numscales
    sfactor = scalefactor^(scale-1);
    chires = resize(ohires, sfactor, 'bicubic');
    
    chires = modcrop(chires, conf.scale); % crop a bit (to simplify scaling issues)
    % Scale down images
    clores = resize(chires, 1/conf.scale, conf.interpolate_kernel);
    midres = resize(clores, conf.upsample_factor, conf.interpolate_kernel);
    features = collect(conf, midres, conf.upsample_factor, conf.filters);
    clear midres

    interpolated = resize(clores, conf.scale, conf.interpolate_kernel);
    clear clores
    patches = cell(size(chires));
    for i = 1:numel(patches) % Remove low frequencies
        patches{i} = chires{i} - interpolated{i};
    end
    clear chires interpolated

    hires = [hires collect(conf, patches, conf.scale, {})];
    
    lores = [lores  features];
end


% PCA dimensionality reduction
C = double(lores * lores');
[V, D] = eig(C);
D = diag(D); % perform PCA on features matrix 
D = cumsum(D) / sum(D);
k = find(D >= 1e-3, 1); % ignore 0.1% energy
V_pca = V(:, k:end); % choose the largest eigenvectors' projection
lores = V_pca' * lores;