tcode= array2table(bptcodevec);
long = cell2table(bplabvec_long);
short = cell2table(bplabvec_short);
bpraw = array2table(bpdata_raw);
bptrans = array2table(bpdata);

writetable(bpraw, '/Users/stanza/Desktop/bpraw.csv')
writetable(tcode, '/Users/stanza/Desktop/tcode.csv')
writetable(short, '/Users/stanza/Desktop/short.csv','Encoding','UTF-8')
writetable(long, '/Users/stanza/Desktop/long.csv','Encoding','UTF-8')
writetable(bptrans, '/Users/stanza/Desktop/bptrans.csv')