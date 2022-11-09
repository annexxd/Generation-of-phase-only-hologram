function [PSNR]=PSNR(img,img_processed)
    [m,n]=size(img);
    img = im2double(img);
    img_processed = im2double(img_processed);
    img = 255*img;
    img_processed = 255*img_processed;
    num_MSE=0;
    %MSE
    for i=1:m
        for j=1:n
            diff = (img(i,j) - img_processed(i,j))^2;
            num_MSE = num_MSE + diff;
        end
    end
    MSE = num_MSE / (m*n) ;   
    PSNR = 10*log10(255*255/MSE);
end
