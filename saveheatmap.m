% Save matrix data so that it can be rendered as heatmap in pgfplots

function saveheatmap(data,labx,laby,filename)    
    fid = fopen(filename,'w');
    for i = 1:length(labx)
        for j = 1:length(laby)
            fprintf(fid, '%d %d %d \n',labx(i),laby(j),data(i,j));
        end
        fprintf(fid,'\n');
    end
    
    fclose(fid);
end