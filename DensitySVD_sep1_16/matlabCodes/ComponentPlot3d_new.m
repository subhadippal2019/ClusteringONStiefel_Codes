
function p1=ComponentPlot3d_new(Component,NiiFile,cutoff,smooting1,FaceColor,EdgeColor,FaceAlpha,gaussSmooth,boxSmooth)
%NiiFile='C:\Users\subha\Box Sync\DTIANALYSIS\selected\600085654611\CN_23FA_0.3Baseline_18\out\12__sub01_component_ica_s1_';
BB=load_nii([NiiFile,'.nii.gz'])
BB=BB.img;
size(BB)
CC=BB(:,:,:);
cdata=CC;
        CC_reduced=reducevolume(CC,smooting1); %%%%  choice of 1 or 2 in this step is crucial
        cdata=CC_reduced;
         cdata = smooth3(cdata,'gaussian',3);
         cdata = smooth3(cdata,'box',boxSmooth);
        cdata = smooth3(cdata,'gaussian',gaussSmooth);
        cdata = smooth3(cdata,'box',boxSmooth+2);
          
       %cdata=reducevolume(cdata,smooting1);
       %cdata = smooth3(cdata,'box',3);
       fv_red = isosurface(cdata,cutoff);
        %fv_red = smooth3(fv_red,'box',3);
        % [faces,verts,colors] = isosurface(cdata,.3,colData);
        p1 = patch(fv_red,'FaceColor',FaceColor,'EdgeColor',EdgeColor,'EdgeAlpha',.31,'FaceAlpha', FaceAlpha);
        %isonormals(cdata,p1)
        %p1 = patch(fv_red,'EdgeColor',EdgeColor,'EdgeAlpha',.31,'FaceAlpha', .3);
      %  p1=surf(fv_red)
          axis tight;
          axis equal;
          camlight
        %camlight(-100,-10)
        lighting gouraud