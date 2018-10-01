tcode= array2table(bptcodevec);
long = cell2table(bplabvec_long);
short = cell2table(bplabvec_short);
bpraw = array2table(bpdata_raw);
bptrans = array2table(bpdata);

writetable(bpraw, '/Users/stanza/Desktop/lasso_r/bpraw.csv')
writetable(tcode, '/Users/stanza/Desktop/lasso_r/tcode.csv')
writetable(short, '/Users/stanza/Desktop/lasso_r/short.csv','Encoding','UTF-8')
writetable(long, '/Users/stanza/Desktop/lasso_r/long.csv','Encoding','UTF-8')
writetable(bptrans, '/Users/stanza/Desktop/lasso_r/bptrans.csv')