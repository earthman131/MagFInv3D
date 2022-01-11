% 保存三维网格模型数据

function savemod(model,nx,ny,nz,outfile4)
fp=fopen(outfile4,'wt');
for i=1:ny
    for j=1:nx
        for k=1:nz
            fprintf(fp,'%f\n',model(i,j,k));
        end
    end
end
fclose(fp);
