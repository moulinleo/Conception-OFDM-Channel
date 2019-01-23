function plotImage(phi, teta, angularSpectrum, limitAmp)

%limitAmp: max and min on first step
    crsz = get(0,'ScreenSize');
    fig = figure;
    panoramicPicture = imread(sprintf('pictureLab2.JPG'));
    im = image(phi,teta,panoramicPicture);
    hold on;
    axis image
    bb=pcolor(phi,teta,angularSpectrum); colorbar;
    set(bb,'FaceAlpha',0.5,'EdgeAlpha',0);
    caxis([limitAmp(2) limitAmp(1)]);
    set(gca, 'DataAspectRatio', [1 2 1])
    set(gca,'YDir','reverse')
    set(gca,'XDir','reverse')

end