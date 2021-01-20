function im_noise = addNoiseToImage(im,ron,gain,nPhotonBackground,QE)

%% ADD NOISE TO IMAGE
if ismatrix(im)
    if license('checkout','statistics_toolbox')        
        im_noise = poissrnd(im + nPhotonBackground) - nPhotonBackground;
        im_noise = QE*im_noise;
        if ron > 0
            im_noise = im_noise + randn(size(im_noise)).*ron/gain;
        end
    else
        buffer    = im + nPhotonBackground;
        image = im + randn(size(im)).*(im + nPhotonBackground);
        index     = image<0;
        image(index) = buffer(index);
        im_noise = QE*image;
        
        if obj.readOutNoise>0
            im_noise = im_noise + randn(size(im_noise)).*ron/gain;
        end
    end          
end
                
