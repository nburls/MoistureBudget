function [Q,V,VQ,U,UQ,P,E,p,dp_3D,dx_rho_3D,dy_rho_3D,lat,lon]=get_moistureflux_fields_interp(infile);



% Reads in:
%
%	float Q(time, lev, lat, lon) ;
%		Q:mdims = 1 ;
%		Q:units = "kg/kg" ;
%		Q:long_name = "Specific humidity" ;
%		Q:cell_methods = "time: mean" ;
%	float V(time, lev, lat, lon) ;
%		V:mdims = 1 ;
%		V:units = "m/s" ;
%		V:long_name = "Meridional wind" ;
%		V:cell_methods = "time: mean" ;
%	float VQ(time, lev, lat, lon) ;
%		VQ:mdims = 1 ;
%		VQ:units = "m/skg/kg" ;
%		VQ:long_name = "Meridional water transport" ;
%		VQ:cell_methods = "time: mean" ;
%	float U(time, lev, lat, lon) ;
%		U:mdims = 1 ;
%		U:units = "m/s" ;
%		U:long_name = "Zonal wind" ;
%		U:cell_methods = "time: mean" ;
%   float UQ(time, lev, lat, lon) ;
%       UQ:units = "m/skg/kg" ;
%       UQ:mdims = 1 ;
%       UQ:cell_methods = "time: mean time: mean" ;
%       UQ:long_name = "Zonal water transport" ;
%
% PLUS
%
%       p: One-dimensional array of pressure level values of vertical
%          dimension ORDERED TOP TO BOTTOM. Units must be Pa. Type float.
%          The first value must be greater than 500 Pa (5mb), and the
%          last value must be less than 100500 Pa (1005mb).
%
%       ps: Two-dimensional (lat,lon) array of surface pressures. Units
%           must be Pa. Type float.
%       READ
%	    float PS(time, lat, lon) ;
%		PS:units = "Pa" ;
%		PS:long_name = "Surface pressure" ;
%		PS:cell_methods = "time: mean" ;
% 	double hyam(lev) ;
%		hyam:long_name = "hybrid A coefficient at layer midpoints" ;
%	double hybm(lev) ;
%		hybm:long_name = "hybrid B coefficient at layer midpoints"

      ncid = netcdf.open(infile,'NOWRITE');
      
      varid = netcdf.inqVarID(ncid,'Q'); 
      Q = netcdf.getVar(ncid,varid,'double'); clear  varid   

      varid = netcdf.inqVarID(ncid,'V'); 
      V = netcdf.getVar(ncid,varid,'double'); clear  varid       
   
      varid = netcdf.inqVarID(ncid,'VQ'); 
      VQ = netcdf.getVar(ncid,varid,'double'); clear  varid      

      varid = netcdf.inqVarID(ncid,'U'); 
      U = netcdf.getVar(ncid,varid,'double'); clear  varid   
      
      varid = netcdf.inqVarID(ncid,'UQ'); 
      UQ = netcdf.getVar(ncid,varid,'double'); clear  varid  

      varid = netcdf.inqVarID(ncid,'PS'); 
      PS = netcdf.getVar(ncid,varid,'double'); clear  varid 

      lat=ncread(infile,'lat');   
      lon=ncread(infile,'lon');
      A=ncread(infile,'hyam');
      B=ncread(infile,'hybm');

      P0 = 100000 ; % reference pressure in PA

      [I J K]=size(V);
      for i=1:I 
        for j=1:J
          for k=1:K
           p(i,j,k)=A(k).*P0+B(k).*PS(i,j);
          end
        end  
      end

      [nlon,nlat,nlvl]=size(V);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % get interface pressures
      p_int_1 = (p(:,:,1:end-1)+p(:,:,2:end))./2;

      % according to CAM4 documention top of the model pressure is 2.917 mb
      p_top = zeros(nlon,nlat)+291.7; 
            
      p_int = zeros(nlon,nlat,nlvl+1);
      p_int(:,:,1)=p_top;
      p_int(:,:,2:end-1)=p_int_1; clear p_int_1
      p_int(:,:,end)=PS;

      dp_3D = zeros(nlon,nlat,nlvl);
      dp_3D =p_int(:,:,2:end)-p_int(:,:,1:end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reads in the average precip field
% 
%	float PRECC(time, lat, lon) ;
%		PRECC:units = "m/s" ;
%		PRECC:long_name = "Convective precipitation rate (liq + ice)" ;
%		PRECC:cell_methods = "time: mean" ;
%	float PRECL(time, lat, lon) ;
%		PRECL:units = "m/s" ;
%		PRECL:long_name = "Large-scale (stable) precipitation rate (liq + ice)" ;
%		PRECL:cell_methods = "time: mean" ;
%	float QFLX(time, lat, lon) ;
%		QFLX:units = "kg/m2/s" ;
%		QFLX:long_name = "Surface water flux" ;
%		QFLX:cell_methods = "time: mean" ;

      ncid = netcdf.open(infile,'NOWRITE');
      
      varid = netcdf.inqVarID(ncid,'PRECC'); 
      precc = netcdf.getVar(ncid,varid,'double'); clear  varid 
      
      varid = netcdf.inqVarID(ncid,'PRECL'); 
      precl = netcdf.getVar(ncid,varid,'double'); clear  varid       
      
      precip_field =  precc + precl;
      
      varid = netcdf.inqVarID(ncid,'QFLX'); 
      evap = netcdf.getVar(ncid,varid,'double'); clear  varid  
      
      netcdf.close(ncid);

% From amwg_diag5.2 > functions_surfaces.ncl
%      ep@long_name = "evap-precip"
%      ep@units = "mm/day"
%      ep@derive_op = "EP=QFLX*8.64e4-(PRECC+PRECL)*8.64e7"        

% convert from m/s to kg/m2/s using water density of 1000kg/m3
P = precip_field.*1000;

% In kg/m2/s 
E = evap; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lat_rho,lon_rho]=meshgrid(lat,lon); 

[I,J]=size(lon_rho); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get longitudes and latitudes of C-Grid u and v points
lon_u=zeros(I+1,J);
lon_u(2:end-1,1:end)=0.5*(lon_rho(1:end-1,1:end)+lon_rho(2:end,1:end));
lon_u(1,1:end)=0.5*(lon_rho(end,1:end)+lon_rho(1,1:end)+360);
lon_u(end,1:end)=0.5*(lon_rho(end,1:end)+lon_rho(1,1:end)+360);

lat_u=zeros(I+1,J);
lat_u(2:end-1,1:end)=0.5*(lat_rho(1:end-1,1:end)+lat_rho(2:end,1:end));
lat_u(1,1:end)=0.5*(lat_rho(end,1:end)+lat_rho(1,1:end));
lat_u(end,1:end)=0.5*(lat_rho(end,1:end)+lat_rho(1,1:end));

lon_v=zeros(I,J+1);
lon_v(1:end,2:end-1)=0.5*(lon_rho(1:end,1:end-1)+lon_rho(1:end,2:end));
lon_v(1:end,1)=0.5*(lon_rho(1:end,end)+lon_rho(1:end,1));
lon_v(1:end,end)=0.5*(lon_rho(1:end,end)+lon_rho(1:end,1));

lat_v=zeros(I,J+1);
lat_v(1:end,2:end-1)=0.5*(lat_rho(1:end,1:end-1)+lat_rho(1:end,2:end));
lat_v(1:end,1)=0.5*(-90+lat_rho(1:end,1));
lat_v(1:end,end)=0.5*(90+lat_rho(1:end,end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx_rho=zeros(I,J);
dx_rho=spheric_dist(lat_u(1:end-1,:),lat_u(2:end,:),...
                       lon_u(1:end-1,:),lon_u(2:end,:));
dy_rho=zeros(I,J);
dy_rho=spheric_dist(lat_v(:,1:end-1,:),lat_v(:,2:end,:),...
                       lon_v(:,1:end-1),lon_v(:,2:end));

area=dy_rho.*dx_rho;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx_rho_3D=repmat(dx_rho,[1 1 nlvl]);
dy_rho_3D=repmat(dy_rho,[1 1 nlvl]);


