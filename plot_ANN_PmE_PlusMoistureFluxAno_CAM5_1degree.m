clear

    infile1 =  './Control_SSTs_CAM5_0.9x1.25_gx1v6.cam.h0.0021-0120._ANN_climatology.nc';
 
    [Moisture_Div,UQBar_div,VQBar_div,VQBar_2D,UQBar_2D,TimeMean_Moisture_Div,UBQB_div,VBQB_div,VBQB_2D,UBQB_2D,Eddy_Moisture_Div,UpQp_div,VpQp_div,VpQp_2D,UpQp_2D,P,E,lat,lon]=get_vertinteg_moistureflux_withdivergence(infile1);
    % P and E in units of kg/m^2/s & moisture div is in m/s 
    clear infile1  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PmE_control=squeeze(P(:,:)-E(:,:)).*86400;       %convert from kg/m^2/s to mm/day -> ./1000.*1000.*86400

Moist_Div_control=squeeze(Moisture_Div(:,:,1)).*86400.*1000;       %convert from m/s to mm/day -> .*1000.*86400

[Z_lon,Z_lat] = meshgrid(lon,lat);

% Define map limits
% latlim = [0 50];
% lonlim = [40 160];
 latlim = [-70 70];
 lonlim = [0 360];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure('Units','centimeters',...
%       'Position',[1 1 18 25],...
%       'PaperPosition',[1 1 18 25],...
%       'PaperUnits','centimeters')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PP(7,:) = [ 0.1300-0.1   0.1100-0.04     0.3347+0.14    0.1577+0.06 ];

   % subplot('position',[PP(j,:)]);

    axesm('robinson','MapLatLimit',latlim,'MapLonLimit',lonlim, ...
   'Frame','on','Grid','on','MeridianLabel','on', ...
   'ParallelLabel','on')

    setm(gca,'mlinelocation',60,'MLabelLocation',60);
    setm(gca,'plinelocation',30,'PLabelLocation',30);
    setm(gca,'frame','off');
    setm(gca,'fontsize',10);
    setm(gca,'FEdgeColor',[0 1 0]);


    [cc,AA1]=contourfm(Z_lat,Z_lon,PmE_control',[-8:1:7],'LineStyle', 'none');
    hold on;
   %contourm(Z_lat,Z_lon,(Z_plio-Z_control)',[0],'k');
    
   % scatterm(sig_lat(1:1:end),sig_lon(1:1:end),0.5,'m','filled');

   % h = quiverm(Z_lat(1:2:end,1:2:end),Z_lon(1:2:end,1:2:end),VQBar_2D(1:2:end,1:2:end,1)',UQBar_2D(1:2:end,1:2:end,1)',1);
   % set(h,'color','m');
   
   
    title('Control P minus E');
   
    coast = load('coast.mat');
    plotm(coast.lat,coast.long,'k','linewidth',1)
    
    i=gca;
    set(i,'Visible','off')

     caxis([-8 8]);

    i=gca;
    set(i,'Visible','off')

    C2=colorbar('hori'); 
    set(C2,'xtick',[-8:2:8],'fontsize',15);

    [C]=text(2,-1.52,'mm/day');
    set(C,'fontsize',10);
    
    clin = clmap(39);
    colormap(clin([1 2 4 6 7 8 9 10 12 13 14 15 16 18 20 21],:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure_name=strcat('Cnt_PmE_years21to120_AnnualMean');
 %figure_name=strcat('PlioceneAnominusCnt_PmE_years21to120_AnnualMean_Global');
  H=gcf;
  set(H,'color','white');
  hold off;
 
 % filename = strcat(figure_name,'.eps');
 % print('-depsc','-tiff',filename);
 % export_fig(figure_name,'-pdf')
  
  handle = getframe(gcf);
  filename = strcat(figure_name,'.png');
  imwrite(handle.cdata,filename,'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %close all
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure('Units','centimeters',...
%       'Position',[1 1 18 25],...
%       'PaperPosition',[1 1 18 25],...
%       'PaperUnits','centimeters')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PP(7,:) = [ 0.1300-0.1   0.1100-0.04     0.3347+0.14    0.1577+0.06 ];

   % subplot('position',[PP(j,:)]);

    axesm('robinson','MapLatLimit',latlim,'MapLonLimit',lonlim, ...
   'Frame','on','Grid','on','MeridianLabel','on', ...
   'ParallelLabel','on')

    setm(gca,'mlinelocation',60,'MLabelLocation',60);
    setm(gca,'plinelocation',30,'PLabelLocation',30);
    setm(gca,'frame','off');
    setm(gca,'fontsize',10);
    setm(gca,'FEdgeColor',[0 1 0]);


    [cc,AA1]=contourfm(Z_lat,Z_lon,-Moist_Div_control',[-8:1:7],'LineStyle', 'none');
    hold on;
   %contourm(Z_lat,Z_lon,(Z_plio-Z_control)',[0],'k');
    
   % scatterm(sig_lat(1:1:end),sig_lon(1:1:end),0.5,'m','filled');

  %  h = quiverm(Z_lat(1:2:end,1:2:end),Z_lon(1:2:end,1:2:end),VQBar_2D(1:2:end,1:2:end,1)',UQBar_2D(1:2:end,1:2:end,1)',1);
  %  set(h,'color','m');
   
   
    title('Control Moisture Divergence');
   
    coast = load('coast.mat');
    plotm(coast.lat,coast.long,'k','linewidth',1)
    
    i=gca;
    set(i,'Visible','off')

     caxis([-8 8]);

    i=gca;
    set(i,'Visible','off')

    C2=colorbar('hori'); 
    set(C2,'xtick',[-8:2:8],'fontsize',15);

    [C]=text(2,-1.52,'mm/day');
    set(C,'fontsize',10);
    
    clin = clmap(39);
    colormap(clin([1 2 4 6 7 8 9 10 12 13 14 15 16 18 20 21],:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure_name=strcat('Cnt_MoistDiv_years21to120_AnnualMean');
 %figure_name=strcat('PlioceneAnominusCnt_PmE_years21to120_AnnualMean_Global');
  H=gcf;
  set(H,'color','white');
  hold off;
 
 % filename = strcat(figure_name,'.eps');
 % print('-depsc','-tiff',filename);
 % export_fig(figure_name,'-pdf')
  
  handle = getframe(gcf);
  filename = strcat(figure_name,'.png');
  imwrite(handle.cdata,filename,'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %close all

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure('Units','centimeters',...
%       'Position',[1 1 18 25],...
%       'PaperPosition',[1 1 18 25],...
%       'PaperUnits','centimeters')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PP(7,:) = [ 0.1300-0.1   0.1100-0.04     0.3347+0.14    0.1577+0.06 ];

   % subplot('position',[PP(j,:)]);

    axesm('robinson','MapLatLimit',latlim,'MapLonLimit',lonlim, ...
   'Frame','on','Grid','on','MeridianLabel','on', ...
   'ParallelLabel','on')

    setm(gca,'mlinelocation',60,'MLabelLocation',60);
    setm(gca,'plinelocation',30,'PLabelLocation',30);
    setm(gca,'frame','off');
    setm(gca,'fontsize',10);
    setm(gca,'FEdgeColor',[0 1 0]);


    [cc,AA1]=contourfm(Z_lat,Z_lon,PmE_control'+Moist_Div_control',[-8:1:7],'LineStyle', 'none');
    hold on;
   %contourm(Z_lat,Z_lon,(Z_plio-Z_control)',[0],'k');
    
   % scatterm(sig_lat(1:1:end),sig_lon(1:1:end),0.5,'m','filled');

  %  h = quiverm(Z_lat(1:2:end,1:2:end),Z_lon(1:2:end,1:2:end),VQBar_2D(1:2:end,1:2:end,1)',UQBar_2D(1:2:end,1:2:end,1)',1);
  %  set(h,'color','m');
   
   
    title('Control P-E minus Moisture Divergence');
   
    coast = load('coast.mat');
    plotm(coast.lat,coast.long,'k','linewidth',1)
    
    i=gca;
    set(i,'Visible','off')

     caxis([-8 8]);

    i=gca;
    set(i,'Visible','off')

    C2=colorbar('hori'); 
    set(C2,'xtick',[-8:2:8],'fontsize',15);

    [C]=text(2,-1.52,'mm/day');
    set(C,'fontsize',10);
    
    clin = clmap(39);
    colormap(clin([1 2 4 6 7 8 9 10 12 13 14 15 16 18 20 21],:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure_name=strcat('Cnt_PmE_minus_MoistDiv_years21to120_AnnualMean');
 %figure_name=strcat('PlioceneAnominusCnt_PmE_years21to120_AnnualMean_Global');
  H=gcf;
  set(H,'color','white');
  hold off;
 
 % filename = strcat(figure_name,'.eps');
 % print('-depsc','-tiff',filename);
 % export_fig(figure_name,'-pdf')
  
  handle = getframe(gcf);
  filename = strcat(figure_name,'.png');
  imwrite(handle.cdata,filename,'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TimeMean_Moist_Div_control=squeeze(TimeMean_Moisture_Div(:,:,1)).*86400.*1000;       %convert from m/s to mm/day -> .*1000.*86400

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure('Units','centimeters',...
%       'Position',[1 1 18 25],...
%       'PaperPosition',[1 1 18 25],...
%       'PaperUnits','centimeters')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PP(7,:) = [ 0.1300-0.1   0.1100-0.04     0.3347+0.14    0.1577+0.06 ];

   % subplot('position',[PP(j,:)]);

    axesm('robinson','MapLatLimit',latlim,'MapLonLimit',lonlim, ...
   'Frame','on','Grid','on','MeridianLabel','on', ...
   'ParallelLabel','on')

    setm(gca,'mlinelocation',60,'MLabelLocation',60);
    setm(gca,'plinelocation',30,'PLabelLocation',30);
    setm(gca,'frame','off');
    setm(gca,'fontsize',10);
    setm(gca,'FEdgeColor',[0 1 0]);


    [cc,AA1]=contourfm(Z_lat,Z_lon,-TimeMean_Moist_Div_control',[-8:1:7],'LineStyle', 'none');
    hold on;
   %contourm(Z_lat,Z_lon,(Z_plio-Z_control)',[0],'k');
    
   % scatterm(sig_lat(1:1:end),sig_lon(1:1:end),0.5,'m','filled');

  %  h = quiverm(Z_lat(1:2:end,1:2:end),Z_lon(1:2:end,1:2:end),VQBar_2D(1:2:end,1:2:end,1)',UQBar_2D(1:2:end,1:2:end,1)',1);
  %  set(h,'color','m');
   
   
    title('Control Time Mean Moisture Divergence');
   
    coast = load('coast.mat');
    plotm(coast.lat,coast.long,'k','linewidth',1)
    
    i=gca;
    set(i,'Visible','off')

     caxis([-8 8]);

    i=gca;
    set(i,'Visible','off')

    C2=colorbar('hori'); 
    set(C2,'xtick',[-8:2:8],'fontsize',15);

    [C]=text(2,-1.52,'mm/day');
    set(C,'fontsize',10);
    
    clin = clmap(39);
    colormap(clin([1 2 4 6 7 8 9 10 12 13 14 15 16 18 20 21],:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure_name=strcat('Cnt_TimeMean_MoistDiv_years21to120_AnnualMean');
 %figure_name=strcat('PlioceneAnominusCnt_PmE_years21to120_AnnualMean_Global');
  H=gcf;
  set(H,'color','white');
  hold off;
 
 % filename = strcat(figure_name,'.eps');
 % print('-depsc','-tiff',filename);
 % export_fig(figure_name,'-pdf')
  
  handle = getframe(gcf);
  filename = strcat(figure_name,'.png');
  imwrite(handle.cdata,filename,'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %close all

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure('Units','centimeters',...
%       'Position',[1 1 18 25],...
%       'PaperPosition',[1 1 18 25],...
%       'PaperUnits','centimeters')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PP(7,:) = [ 0.1300-0.1   0.1100-0.04     0.3347+0.14    0.1577+0.06 ];

   % subplot('position',[PP(j,:)]);

    axesm('robinson','MapLatLimit',latlim,'MapLonLimit',lonlim, ...
   'Frame','on','Grid','on','MeridianLabel','on', ...
   'ParallelLabel','on')

    setm(gca,'mlinelocation',60,'MLabelLocation',60);
    setm(gca,'plinelocation',30,'PLabelLocation',30);
    setm(gca,'frame','off');
    setm(gca,'fontsize',10);
    setm(gca,'FEdgeColor',[0 1 0]);


    [cc,AA1]=contourfm(Z_lat,Z_lon,PmE_control'+TimeMean_Moist_Div_control',[-8:1:7],'LineStyle', 'none');
    hold on;
   %contourm(Z_lat,Z_lon,(Z_plio-Z_control)',[0],'k');
    
   % scatterm(sig_lat(1:1:end),sig_lon(1:1:end),0.5,'m','filled');

  %  h = quiverm(Z_lat(1:2:end,1:2:end),Z_lon(1:2:end,1:2:end),VQBar_2D(1:2:end,1:2:end,1)',UQBar_2D(1:2:end,1:2:end,1)',1);
  %  set(h,'color','m');
   
   
    title('Control P-E minus Time Mean Moisture Divergence');
   
    coast = load('coast.mat');
    plotm(coast.lat,coast.long,'k','linewidth',1)
    
    i=gca;
    set(i,'Visible','off')

     caxis([-8 8]);

    i=gca;
    set(i,'Visible','off')

    C2=colorbar('hori'); 
    set(C2,'xtick',[-8:2:8],'fontsize',15);

    [C]=text(2,-1.52,'mm/day');
    set(C,'fontsize',10);
    
    clin = clmap(39);
    colormap(clin([1 2 4 6 7 8 9 10 12 13 14 15 16 18 20 21],:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure_name=strcat('Cnt_PmE_minus_TimeMean_MoistDiv_years21to120_AnnualMean');
 %figure_name=strcat('PlioceneAnominusCnt_PmE_years21to120_AnnualMean_Global');
  H=gcf;
  set(H,'color','white');
  hold off;
 
 % filename = strcat(figure_name,'.eps');
 % print('-depsc','-tiff',filename);
 % export_fig(figure_name,'-pdf')
  
  handle = getframe(gcf);
  filename = strcat(figure_name,'.png');
  imwrite(handle.cdata,filename,'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Eddy_Moist_Div_control=squeeze(Eddy_Moisture_Div(:,:,1)).*86400.*1000;       %convert from m/s to mm/day -> .*1000.*86400

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure('Units','centimeters',...
%       'Position',[1 1 18 25],...
%       'PaperPosition',[1 1 18 25],...
%       'PaperUnits','centimeters')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PP(7,:) = [ 0.1300-0.1   0.1100-0.04     0.3347+0.14    0.1577+0.06 ];

   % subplot('position',[PP(j,:)]);

    axesm('robinson','MapLatLimit',latlim,'MapLonLimit',lonlim, ...
   'Frame','on','Grid','on','MeridianLabel','on', ...
   'ParallelLabel','on')

    setm(gca,'mlinelocation',60,'MLabelLocation',60);
    setm(gca,'plinelocation',30,'PLabelLocation',30);
    setm(gca,'frame','off');
    setm(gca,'fontsize',10);
    setm(gca,'FEdgeColor',[0 1 0]);


    [cc,AA1]=contourfm(Z_lat,Z_lon,-Eddy_Moist_Div_control',[-8:1:7],'LineStyle', 'none');
    hold on;
   %contourm(Z_lat,Z_lon,(Z_plio-Z_control)',[0],'k');
    
   % scatterm(sig_lat(1:1:end),sig_lon(1:1:end),0.5,'m','filled');

  %  h = quiverm(Z_lat(1:2:end,1:2:end),Z_lon(1:2:end,1:2:end),VQBar_2D(1:2:end,1:2:end,1)',UQBar_2D(1:2:end,1:2:end,1)',1);
  %  set(h,'color','m');
   
   
    title('Control Eddy Moisture Divergence');
   
    coast = load('coast.mat');
    plotm(coast.lat,coast.long,'k','linewidth',1)
    
    i=gca;
    set(i,'Visible','off')

     caxis([-8 8]);

    i=gca;
    set(i,'Visible','off')

    C2=colorbar('hori'); 
    set(C2,'xtick',[-8:2:8],'fontsize',15);

    [C]=text(2,-1.52,'mm/day');
    set(C,'fontsize',10);
    
    clin = clmap(39);
    colormap(clin([1 2 4 6 7 8 9 10 12 13 14 15 16 18 20 21],:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure_name=strcat('Cnt_Eddy_MoistDiv_years21to120_AnnualMean');
 %figure_name=strcat('PlioceneAnominusCnt_PmE_years21to120_AnnualMean_Global');
  H=gcf;
  set(H,'color','white');
  hold off;
 
 % filename = strcat(figure_name,'.eps');
 % print('-depsc','-tiff',filename);
 % export_fig(figure_name,'-pdf')
  
  handle = getframe(gcf);
  filename = strcat(figure_name,'.png');
  imwrite(handle.cdata,filename,'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %close all

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure('Units','centimeters',...
%       'Position',[1 1 18 25],...
%       'PaperPosition',[1 1 18 25],...
%       'PaperUnits','centimeters')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PP(7,:) = [ 0.1300-0.1   0.1100-0.04     0.3347+0.14    0.1577+0.06 ];

   % subplot('position',[PP(j,:)]);

    axesm('robinson','MapLatLimit',latlim,'MapLonLimit',lonlim, ...
   'Frame','on','Grid','on','MeridianLabel','on', ...
   'ParallelLabel','on')

    setm(gca,'mlinelocation',60,'MLabelLocation',60);
    setm(gca,'plinelocation',30,'PLabelLocation',30);
    setm(gca,'frame','off');
    setm(gca,'fontsize',10);
    setm(gca,'FEdgeColor',[0 1 0]);


    [cc,AA1]=contourfm(Z_lat,Z_lon,PmE_control'+Eddy_Moist_Div_control',[-8:1:7],'LineStyle', 'none');
    hold on;
   %contourm(Z_lat,Z_lon,(Z_plio-Z_control)',[0],'k');
    
   % scatterm(sig_lat(1:1:end),sig_lon(1:1:end),0.5,'m','filled');

  %  h = quiverm(Z_lat(1:2:end,1:2:end),Z_lon(1:2:end,1:2:end),VQBar_2D(1:2:end,1:2:end,1)',UQBar_2D(1:2:end,1:2:end,1)',1);
  %  set(h,'color','m');
   
   
    title('Control P-E minus Eddy Moisture Divergence');
   
    coast = load('coast.mat');
    plotm(coast.lat,coast.long,'k','linewidth',1)
    
    i=gca;
    set(i,'Visible','off')

     caxis([-8 8]);

    i=gca;
    set(i,'Visible','off')

    C2=colorbar('hori'); 
    set(C2,'xtick',[-8:2:8],'fontsize',15);

    [C]=text(2,-1.52,'mm/day');
    set(C,'fontsize',10);
    
    clin = clmap(39);
    colormap(clin([1 2 4 6 7 8 9 10 12 13 14 15 16 18 20 21],:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure_name=strcat('Cnt_PmE_minus_Eddy_MoistDiv_years21to120_AnnualMean');
 %figure_name=strcat('PlioceneAnominusCnt_PmE_years21to120_AnnualMean_Global');
  H=gcf;
  set(H,'color','white');
  hold off;
 
 % filename = strcat(figure_name,'.eps');
 % print('-depsc','-tiff',filename);
 % export_fig(figure_name,'-pdf')
  
  handle = getframe(gcf);
  filename = strcat(figure_name,'.png');
  imwrite(handle.cdata,filename,'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %close all



