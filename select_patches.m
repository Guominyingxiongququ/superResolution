%% f and g are two image sets 
%  this function is used to select the patches 
%  which has the worst performance based on PSNR


function patch_list = select_patches(conf)

    for j = 2:numel(conf.results{1})
        results_truth = cell(numel(conf.results), 1);
        results_test = cell(numel(conf.results), 1);

        for i = 1:numel(conf.results)
            results_truth{i} = conf.results{i}{1};
            results_test{i} = conf.results{i}{j};
            
%             score = calc_performance(conf.results{i}{1}, conf.results{i}{j});
%             scores(i,j) = score;
%             fprintf(fid, '<TD><A HREF=%s>%.1f</A>(%.2f)</TD>\n', ...
%                 esc(conf.results{i}{j}), score, conf.countedtime(j-1,i));
%             fprintf(' : %.1f dB', score)
        end
    end

    % how to select patches from a image
    % f -> global truth
    % g -> experiment result
    windows = load_images_patches(glob('training_set', pattern_str),conf);
    patches = collect(conf, windows, conf.scale, {});

end
% 
% %% collect file names from result 
% function result  = collect_fileName(c)
% 
% end
% 
% function res = calc_PeakSNR(f, g)
% F = im2double(imread(f)); % original
% G = im2double(imread(g)); % distorted
% E = F - G; % error signal
% N = numel(E); % Assume the original signal is at peak (|F|=1)
% res = 10*log10( N / sum(E(:).^2) );
% 
