function [im] = reconstructFace(featVec,meanFace,eigenfaces,dims)
    recon = meanFace;
    for i = 1:size(eigenfaces,2)
        recon = recon+(eigenfaces(:,i).*featVec(i));
    end
    im = reshape(recon,dims);
end