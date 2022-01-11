function savebln( norm_error,outfile ) 
n = length(norm_error);
fp=fopen(outfile,'wt');
fprintf(fp,'%d\t',n);
fprintf(fp,'%d\n',0);
for i=1:n
    fprintf(fp,'%d\t%f\n',i,norm_error(i));
end
fclose(fp);