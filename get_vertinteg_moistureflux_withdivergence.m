function [Moisture_Div,UQBar_div,VQBar_div,VQBar_2D,UQBar_2D,TimeMean_Moisture_Div,UBQB_div,VBQB_div,VBQB_2D,UBQB_2D,Eddy_Moisture_Div,UpQpBar_div,VpQpBar_div,VpQpBar_2D,UpQpBar_2D,P,E,lat,lon]=get_vertinteg_moistureflux_withdivergence(infile)

%
% Reads in:
%
%	float VQ(time, lev, lat, lon) ;
%		VQ:mdims = 1 ;
%		VQ:units = "m/skg/kg" ;
%		VQ:long_name = "Meridional water transport" ;
%		VQ:cell_methods = "time: mean" ;
%   float UQ(time, lev, lat, lon) ;
%       UQ:units = "m/skg/kg" ;
%       UQ:mdims = 1 ;
%       UQ:cell_methods = "time: mean time: mean" ;
%       UQ:long_name = "Zonal water transport" ;
%
%
% Vertically interpolated onto uniform pressure levels
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % From constants.F90 
   rho_air   = 1.2;                   % ambient air density (kg/m^3)
   rho_sw    = 4.1/3.996;             % density of salt water (g/cm^3)
   rho_sw_SIunits = rho_sw*1000;      % density of salt water (kg/cm^3)
   rho_fw = 1.0;                      % density of fresh water (g/cm^3)
   rho_fw_SIunits = rho_fw*1000;      % density of fresh water (kg/cm^3)
   
   g=9.81 %gravitational constant
 
   [Q,V,VQBar,U,UQBar,P,E,p,dp,dx,dy,lat,lon]=get_moistureflux_fields(infile);

   [nlon,nlat,nlvl]=size(Q);
   
   dx_midy = squeeze(dx(:,1:end-1,1)+dx(:,2:end,1))./2;
   dy_midx = squeeze(dy(1:end-1,:,1)+dy(2:end,:,1))./2; 
   
   dx_mid = squeeze(dx(1:end-1,:,1)+dx(2:end,:,1))./2;
   dy_mid = squeeze(dy(:,1:end-1,1)+dy(:,2:end,1))./2; 
   
   area=squeeze(dx(:,:,1).*dy(:,:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VQBar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vertically integrated

  VQBardp = VQBar;
  VQBardp = VQBar.*dp;  % kg/kg m/s Pa

   VQBardp_noNaN=VQBardp;
   I=find(isnan(VQBardp_noNaN)==1);
   VQBardp_noNaN(I)=0; clear I

   VQBar_2D = zeros(nlon,nlat);
   VQBar_2D = sum(VQBardp_noNaN,3)./sum(dp,3).*1000; %convert to g/kg m/s
   
   VQBar_vertinterp = (1./g).*sum(VQBardp_noNaN,3); % in kg/m/s (kg/kg m/s Pa s^2/m -> kg/kg m/s kg/m/s^2 s^2/m) 

% VQBar_vertinterp has the units of kg/m/s 
% convert to m^2/s by dividing by water density and then multiply by y or x
% distance (m) to get m^3/s

   VQBar_vertinterp_mid = VQBar_vertinterp(:,1:end-1)+((VQBar_vertinterp(:,2:end)-VQBar_vertinterp(:,1:end-1))./dy_mid).*0.5.*squeeze(dy(:,1:end-1,1));

   VQBar_transport = zeros(nlon,nlat+1);
   VQBar_transport(:,2:end-1)=VQBar_vertinterp_mid./rho_fw_SIunits.*dx_midy; % m^3/s

   VQBar_div=VQBar_transport(:,2:end)-VQBar_transport(:,1:end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UQBar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vertically integrated

  UQBardp = UQBar;
  UQBardp = UQBar.*dp;  % kg/kg m/s Pa

   UQBardp_noNaN=UQBardp;
   I=find(isnan(UQBardp_noNaN)==1);
   UQBardp_noNaN(I)=0; clear I

   UQBar_2D = zeros(nlon,nlat);
   UQBar_2D = sum(UQBardp_noNaN,3)./sum(dp,3).*1000; %convert to g/kg m/s
   
   UQBar_vertinterp = (1./g).*sum(UQBardp_noNaN,3); % in kg/m/s (kg/kg m/s Pa s^2/m -> kg/kg m/s kg/m/s^2 s^2/m) 
% UQBar_vertinterp has the units of kg/m/s 
% convert to m^2/s by dividing by water density and then multiply by y or x
% distance (m) to get m^3/s

   UQBar_vertinterp_mid = UQBar_vertinterp(1:end-1,:)+((UQBar_vertinterp(2:end,:)-UQBar_vertinterp(1:end-1,:))./dx_mid).*0.5.*squeeze(dx(1:end-1,:,1));
   UQBar_vertinterp_wrap = UQBar_vertinterp(end,:)+((UQBar_vertinterp(1,:)-UQBar_vertinterp(end,:))./(squeeze(dx(1,:,1)+dx(end,:,1))./2)).*0.5.*squeeze(dx(end,:,1));

   UQBar_transport = zeros(nlon+1,nlat);
   UQBar_transport(2:end-1,:)=UQBar_vertinterp_mid./rho_fw_SIunits.*dy_midx; % m^3/s
      UQBar_transport(1,:)=UQBar_vertinterp_wrap./rho_fw_SIunits.*(squeeze(dy(1,:,1)+dy(end,:,1))./2); % m^3/s
      UQBar_transport(end,:)=UQBar_vertinterp_wrap./rho_fw_SIunits.*(squeeze(dy(1,:,1)+dy(end,:,1))./2); % m^3/s

   UQBar_div=UQBar_transport(2:end,:)-UQBar_transport(1:end-1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Moisture_Div = (UQBar_div+VQBar_div)./area; %m/s 
Moisture_Div(:,1) = 0; Moisture_Div(:,end) = 0; % set moisture div to zero at poles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VBarQBar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% VBarQBar
%
   VBQB_3D = V.*Q;

% Vertically integrated
   VBQBdp = VBQB_3D;
   VBQBdp = VBQB_3D.*dp;  % kg/kg m/s Pa OR kg/kg m/s 

   VBQBdp_noNaN=VBQBdp;
   I=find(isnan(VBQBdp_noNaN)==1);
   VBQBdp_noNaN(I)=0; clear I

   VBQB_2D = zeros(nlon,nlat);
   VBQB_2D = sum(VBQBdp_noNaN,3)./sum(dp,3).*1000; %convert to g/kg m/s
% Vertically integrated
   VBQBdp = VBQB_3D;
   VBQBdp = VBQB_3D.*dp;  % kg/kg m/s Pa

   VBQBdp_noNaN=VBQBdp;
   I=find(isnan(VBQBdp_noNaN)==1);
   VBQBdp_noNaN(I)=0; clear I

   VBQB_2D = zeros(nlon,nlat);
   VBQB_2D = sum(VBQBdp_noNaN,3)./sum(dp,3).*1000; %convert to g/kg m/s
   
   VBQB_vertinterp = (1./g).*sum(VBQBdp_noNaN,3); % in kg/m/s (kg/kg m/s Pa s^2/m -> kg/kg m/s kg/m/s^2 s^2/m) 

% VQBar_vertinterp has the units of kg/m/s 
% convert to m^2/s by dividing by water density and then multiply by y or x
% distance (m) to get m^3/s

   VBQB_vertinterp_mid = VBQB_vertinterp(:,1:end-1)+((VBQB_vertinterp(:,2:end)-VBQB_vertinterp(:,1:end-1))./dy_mid).*0.5.*squeeze(dy(:,1:end-1,1));

   VBQB_transport = zeros(nlon,nlat+1);
   VBQB_transport(:,2:end-1)=VBQB_vertinterp_mid./rho_fw_SIunits.*dx_midy; % m^3/s

   VBQB_div=VBQB_transport(:,2:end)-VBQB_transport(:,1:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UBarQBar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% UBarQBar
%
   UBQB_3D = U.*Q;

% Uertically integrated
   UBQBdp = UBQB_3D;
   UBQBdp = UBQB_3D.*dp;  % kg/kg m/s Pa OR kg/kg m/s 

   UBQBdp_noNaN=UBQBdp;
   I=find(isnan(UBQBdp_noNaN)==1);
   UBQBdp_noNaN(I)=0; clear I

   UBQB_2D = zeros(nlon,nlat);
   UBQB_2D = sum(UBQBdp_noNaN,3)./sum(dp,3).*1000; %conUert to g/kg m/s
% Uertically integrated
   UBQBdp = UBQB_3D;
   UBQBdp = UBQB_3D.*dp;  % kg/kg m/s Pa

   UBQBdp_noNaN=UBQBdp;
   I=find(isnan(UBQBdp_noNaN)==1);
   UBQBdp_noNaN(I)=0; clear I

   UBQB_2D = zeros(nlon,nlat);
   UBQB_2D = sum(UBQBdp_noNaN,3)./sum(dp,3).*1000; %conUert to g/kg m/s
   
   UBQB_vertinterp = (1./g).*sum(UBQBdp_noNaN,3); % in kg/m/s (kg/kg m/s Pa s^2/m -> kg/kg m/s kg/m/s^2 s^2/m) 
% UBQB_vertinterp has the units of kg/m/s 
% convert to m^2/s by dividing by water density and then multiply by y or x
% distance (m) to get m^3/s

   UBQB_vertinterp_mid = UBQB_vertinterp(1:end-1,:)+((UBQB_vertinterp(2:end,:)-UBQB_vertinterp(1:end-1,:))./dx_mid).*0.5.*squeeze(dx(1:end-1,:,1));
   UBQB_vertinterp_wrap = UBQB_vertinterp(end,:)+((UBQB_vertinterp(1,:)-UBQB_vertinterp(end,:))./(squeeze(dx(1,:,1)+dx(end,:,1))./2)).*0.5.*squeeze(dx(end,:,1));

   UBQB_transport = zeros(nlon+1,nlat);
   UBQB_transport(2:end-1,:)=UBQB_vertinterp_mid./rho_fw_SIunits.*dy_midx; % m^3/s
   UBQB_transport(1,:)=UBQB_vertinterp_wrap./rho_fw_SIunits.*(squeeze(dy(1,:,1)+dy(end,:,1))./2); % m^3/s
   UBQB_transport(end,:)=UBQB_vertinterp_wrap./rho_fw_SIunits.*(squeeze(dy(1,:,1)+dy(end,:,1))./2); % m^3/s

   UBQB_div=UBQB_transport(2:end,:)-UBQB_transport(1:end-1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TimeMean_Moisture_Div = (UBQB_div+VBQB_div)./area; %m/s 
TimeMean_Moisture_Div(:,1) = 0; TimeMean_Moisture_Div(:,end) = 0; % set moisture div to zero at poles
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Eddy components of the vertically integrated moister flux and moister flux divergence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VpQpBar_2D = VQBar_2D-VBQB_2D;  
VpQpBar_div=VQBar_div-VBQB_div;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Eddy components of the Uertically integrated moister flux and moister flux diUergence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UpQpBar_2D = UQBar_2D-UBQB_2D;
UpQpBar_div=UQBar_div-UBQB_div;

Eddy_Moisture_Div=Moisture_Div-TimeMean_Moisture_Div