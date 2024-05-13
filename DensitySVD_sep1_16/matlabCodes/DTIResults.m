addpath(genpath('C:\Users\subha\Dropbox\Subhadip Pal\MatlabToolbox'))
cd 'C:\Users\subha\Dropbox\Subhadip Pal\Matlab Codes\';



FA_file='C:\Users\subha\Box Sync\DTIANALYSIS\600085654611\data_FA.nii';
mask_file='C:\Users\subha\Desktop\DTIDATA1\CN_20FA_03Baseline_13\mask';

NiiFile='C:\Users\subha\Box Sync\DTIANALYSIS\selected\600085654611\CN_23FA_0.3Baseline_18\P'
FA_file='C:\Users\subha\Box Sync\DTIANALYSIS\600085654611\data_FA.nii';

%[h, colors]=PlotBrainImageformFAnii(FA_file,.2,'Color2','b','b',.2,2);
[h, colors] = PlotBrainImageformFAnii(FA_file,.2,'Color','b','b',.1,2,.1,'box');
ComponentPlot3d(3,NiiFile,.1,2,'b','b',1,3,3);


custerNUm=9
NiiFile1=['C:\Users\subha\Desktop\data_stiefel\new_gen_L1_0.31_',num2str(custerNUm)];
resolution='-r200';
[h, colors]=PlotBrainImageformFAnii(FA_file,.12,'Color2','b','b',.1,2);
set(gca, 'XTick', [], 'YTick', [], 'ZTick',[]);set(gca,'Visible','off');
colormap('copper')

ComponentPlot3d_new(1,NiiFile1,.04,2,'b','none',.4,3,3);
ComponentPlot3d_new(1,NiiFile1,.06,2,'b','none',.9,3,3);
box off;
%Component,NiiFile,cutoff,smooting1,FaceColor,EdgeColor,FaceAl
%save_path='C:\Users\subha\Dropbox\projects\ClusteringDTIonStiefel\RealData\DTIPLOTS\Cluster_';
%print('C:\Users\subha\Dropbox\projects\ClusteringDTIonStiefel\RealData\DTIPLOTS\TOP_Cluster3.png','-dpng','-r400')
view([0 0 1])  
print([save_path,num2str(custerNUm),'Top.png'],'-dpng','-r200');

   az = 90;el = 10;view(az, el);
  %print('C:\Users\subha\Dropbox\projects\ClusteringDTIonStiefel\RealData\DTIPLOTS\Back_Cluster3.png','-dpng','-r400')
   print([save_path,num2str(custerNUm),'Front.png'],'-dpng','-r200')
  
  az = 180;el = 10;view(az, el);
  print([save_path,num2str(custerNUm),'Side.png'],'-dpng','-r200')
  
  close(gcf)
  % Y Axis View
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
custerNUm=11
NiiFile1=['C:\Users\subha\Desktop\data_stiefel\new_gen_L1_0.31_',num2str(custerNUm)];
save_path='C:\Users\subha\Dropbox\projects\ClusteringDTIonStiefel\AOS\figures\DTI_results\Cluster_';
%print('C:\Users\s

[h, colors]=PlotBrainImageformFAnii(FA_file,.12,'Color2','b','b',.1,2);
set(gca, 'XTick', [], 'YTick', [], 'ZTick',[])
colormap('copper');set(gca,'Visible','off');
ComponentPlot3d_new(1,NiiFile1,.08,2,'b','none',.4,3,3);
ComponentPlot3d_new(1,NiiFile1,.1,2,'b','none',.9,3,3);
box off;
%Component,NiiFile,cutoff,smooting1,FaceColor,EdgeColor,FaceAl
%save_path='C:\Users\subha\Dropbox\projects\ClusteringDTIonStiefel\RealData\DTIPLOTS\Cluster_';
%print('C:\Users\subha\Dropbox\projects\ClusteringDTIonStiefel\RealData\DTIPLOTS\TOP_Cluster3.png','-dpng','-r400')
view([0 0 1])  
print([save_path,num2str(custerNUm),'Top.png'],'-dpng','-r200')

 
  az = 90;el = 10;view(az, el);
  %print('C:\Users\subha\Dropbox\projects\ClusteringDTIonStiefel\RealData\DTIPLOTS\Back_Cluster3.png','-dpng','-r400')
   print([save_path,num2str(custerNUm),'Front.png'],'-dpng','-r200')
  
  az = 180;el = 10;view(az, el);
  print([save_path,num2str(custerNUm),'Side.png'],'-dpng','-r200')
  
  close(gcf)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  custerNUm=3
NiiFile1=['C:\Users\subha\Desktop\data_stiefel\new_gen_L1_0.31_',num2str(custerNUm)];

[h, colors]=PlotBrainImageformFAnii(FA_file,.12,'Color2','b','b',.1,2);
set(gca, 'XTick', [], 'YTick', [], 'ZTick',[])
colormap('copper');set(gca,'Visible','off');
ComponentPlot3d_new(1,NiiFile1,.1,2,'b','none',.4,3,3);
ComponentPlot3d_new(1,NiiFile1,.12,2,'b','none',.9,3,3);
box off;
%Component,NiiFile,cutoff,smooting1,FaceColor,EdgeColor,FaceAl
save_path='C:\Users\subha\Dropbox\projects\ClusteringDTIonStiefel\RealData\DTIPLOTS\Cluster_';
%print('C:\Users\subha\Dropbox\projects\ClusteringDTIonStiefel\RealData\DTIPLOTS\TOP_Cluster3.png','-dpng','-r400')
view([0 0 1])  
print([save_path,num2str(custerNUm),'Top.png'],'-dpng','-r200')

 
  az = 90;el = 10;view(az, el);
  %print('C:\Users\subha\Dropbox\projects\ClusteringDTIonStiefel\RealData\DTIPLOTS\Back_Cluster3.png','-dpng','-r400')
   print([save_path,num2str(custerNUm),'Front.png'],'-dpng','-r200')
  
  az = 180;el = 10;view(az, el);
  print([save_path,num2str(custerNUm),'Side.png'],'-dpng','-r200')
  
  close(gcf)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
custerNUm=2
NiiFile1=['C:\Users\subha\Desktop\data_stiefel\new_gen_L1_0.31_',num2str(custerNUm)];

[h, colors]=PlotBrainImageformFAnii(FA_file,.12,'Color2','b','b',.1,2);
set(gca, 'XTick', [], 'YTick', [], 'ZTick',[])
colormap('copper');set(gca,'Visible','off');
ComponentPlot3d_new(1,NiiFile1,.06,2,'b','none',.4,3,3);
ComponentPlot3d_new(1,NiiFile1,.08,2,'b','none',.9,3,3);
box off;
%Component,NiiFile,cutoff,smooting1,FaceColor,EdgeColor,FaceAl
%save_path='C:\Users\subha\Dropbox\projects\ClusteringDTIonStiefel\RealData\DTIPLOTS\Cluster_';
%print('C:\Users\subha\Dropbox\projects\ClusteringDTIonStiefel\RealData\DTIPLOTS\TOP_Cluster3.png','-dpng','-r400')
view([0 0 1])  
print([save_path,num2str(custerNUm),'Top.png'],'-dpng','-r200')

 
  az = 90;el = 10;view(az, el);
  %print('C:\Users\subha\Dropbox\projects\ClusteringDTIonStiefel\RealData\DTIPLOTS\Back_Cluster3.png','-dpng','-r400')
   print([save_path,num2str(custerNUm),'Fron.png'],'-dpng','-r200')
  
  az = 180;el = 10;view(az, el);
  print([save_path,num2str(custerNUm),'Side.png'],'-dpng','-r200')
  
  close(gcf)
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
custerNUm=4
NiiFile1=['C:\Users\subha\Desktop\data_stiefel\new_gen_L1_0.31_',num2str(custerNUm)];
    
[h, colors]=PlotBrainImageformFAnii(FA_file,.12,'Color2','b','b',.1,2);
 set(gca,   'Box', 'off');
set(gca,'color','red');set(gca,'Visible','off');
%set(gca,'XColor',[0,0,0],'YColor',[0,0,0],'TickDir','out')
set(gca, 'Box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[],'ZTickLabel',[],'ZTick',[])
colormap('copper')
ComponentPlot3d_new(1,NiiFile1,.1,2,'b','none',.4,3,3);
ComponentPlot3d_new(1,NiiFile1,.12,2,'b','none',.9,3,3);
%Component,NiiFile,cutoff,smooting1,FaceColor,EdgeColor,FaceAl
save_path='C:\Users\subha\Dropbox\projects\ClusteringDTIonStiefel\AOS\figures\DTI_results\Cluster_';
%print('C:\Users\subha\Dropbox\projects\ClusteringDTIonStiefel\RealData\DTIPLOTS\TOP_Cluster3.png','-dpng','-r400')
view([0 0 1])  
print([save_path,num2str(custerNUm),'Top.png'],'-dpng','-r200')
    
 
  az = 90;el = 10;view(az, el);
  %print('C:\Users\subha\Dropbox\projects\ClusteringDTIonStiefel\RealData\DTIPLOTS\Back_Cluster3.png','-dpng','-r400')
   print([save_path,num2str(custerNUm),'Front.png'],'-dpng','-r200')
  
  az = 180;el = 10;view(az, el);
  print([save_path,num2str(custerNUm),'Side.png'],'-dpng','-r200')
  
  close(gcf)
  