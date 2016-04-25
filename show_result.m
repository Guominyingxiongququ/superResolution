function show_result(conf_set)
    set_num = size(conf_set,2);
    qmkdir(strcat('comparison_', conf_set{10}.result_dir));
    fid = fopen(fullfile(conf_set{1}.result_dir, 'index.html'), 'wt');
    fprintf(fid, ...
        '<HTML><HEAD><TITLE>Super-Resolution Summary</TITLE></HEAD><BODY>');
    fprintf(fid, '<H1>Simulation results</H1>\n');
    fprintf(fid, '<TABLE border="1">\n');
    fprintf(fid, '<TR>');
    fprintf(fid, '<TD>%s</TD>', ' ');
%     for i = 1:conf_set{1}.setNum
    for i = 1:10
        s1 = num2str(i);
        s2 = 'set_'
        s3 = strcat(s2,s1);
        fprintf(fid, '<TD>%s</TD>', s3);
    end
    fprintf(fid, '<TD>PSNR</TD>', s3);
    fprintf(fid, '</TR>\n');
    fprintf(fid, '<TR>');
    fprintf(fid, '<TD>%s</TD>', 'average PSNR');
    sum_PSNR = 0.0;
    for i = 1:10
        mean_PSNR = mean(conf_set{i}.scores(:,3))
        fprintf(fid, '<TD>%6.2f</TD>',mean_PSNR);
        sum_PSNR = sum_PSNR + mean_PSNR;
    end
    fprintf(fid, '<TD>%6.2f</TD>',sum_PSNR/10.0);
    fprintf(fid, '</TR>\n');
    fprintf(fid, '<TR>');
    fprintf(fid, '<TD>%s</TD>', 'COV PSNR');
    
    sum_COV = 0.0;
    for i = 1:10
        mean_COV = cov(conf_set{i}.scores(:,3))
        fprintf(fid, '<TD>%6.2f</TD>', mean_COV);
        sum_COV = sum_COV + mean_COV;
    end
    fprintf(fid, '<TD>%6.2f</TD>',sum_COV/10.0);
    fprintf(fid, '</TR>\n');
    fprintf(fid, '<TR>');
    fprintf(fid, '<TD>%s</TD>', 'patch_num');
    for i = 1:10
        fprintf(fid, '<TD>%5d</TD>', conf_set{i}.patch_num);
    end
    fprintf(fid, '</TR>\n');
    
%     fprintf(fid, '<H1>Simulation parameters</H1>\n<TABLE border="1">\n');
    fprintf(fid, sprintf('<TR><TD>Scaling factor<TD>x%d</TR>\n', conf_set{1}.scale));
    fprintf(fid, sprintf('<TR><TD>window num<TD>x%d</TR>\n', conf_set{1}.window_num));
    fprintf(fid, '</TABLE>\n');
    
    fprintf(fid, '%s\n', datestr(now));
    fprintf(fid, '</BODY></HTML>\n');
    fclose(fid);
    fprintf('\n');
end