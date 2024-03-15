datadir=['/home/weng/Works/Softwares/sem2dpack/EXAMPLES/Vw_test_zhu'];
grid = sem2d_read_specgrid(datadir);

% Read fault data
data = sem2d_read_fault('Flt01');
% Plot fault slipda
figure(1)
for i=400:200:data.nt
    tt=i*data.dt;
    formatSpec = "The current time is: %d";
    str = sprintf(formatSpec,tt);
     plot(data.x/1e3,(data.d(:,i)'+ abs(min(data.d(:,i)))),'DisplayName',str);
%     plot(data.x/1e3,(data.d(:,i)'),'DisplayName',str);  % non-dimension 
    xlim([0 inf])
    ylim([-inf inf])
    xlabel('Along strike distance (km)')
    ylabel('Slip (m)')
    set(gca,'FontSize',18,'Fontname', 'Times New Roman')
    hold on
end
hold off
% % Plot fault slip rate
figure(2)
% % for i=400:1600:data.nt
% %     tt=i*data.dt;
% %     formatSpec = "The current time is: %d";
% %     str = sprintf(formatSpec,tt);
% %     plot(data.x/10,(data.v(:,i)'/2),'DisplayName',str);
% %     xlim([0 inf])
% %     ylim([0 6])
% %     hold on
% % end

plot(data.x/1000,(data.v(:,3600)'));
     xlim([-inf inf])
     ylim([-inf inf])
     xlabel('Along strike distance (km)')
     ylabel('Slip rate (m/s)')
     set(gca,'FontSize',18,'Fontname', 'Times New Roman')
figure (3)     
plot(data.x/1000,(data.mu(:,3600)'));
     xlim([-inf inf])
     ylim([-inf inf])
     xlabel('Along strike distance (km)')
     ylabel('Friction coefficient')
     set(gca,'FontSize',18,'Fontname', 'Times New Roman')
% 
% %Plot slip rate
figure(4)
mesh([0:data.nt-1]*data.dt,data.x/1e3,data.v);
xlabel('Time (s)')
ylabel('Along strike distance (km)')
zlabel('Slip rate (m/s)')
set(gca,'FontSize',18,'Fontname', 'Times New Roman')


%plots a spectral element grid
%h=sem2d_plot_grid(grid);

s11 = sem2d_snapshot_read('vx',7,datadir);


% % Plot velocity a horizontal fault located at x=0:
%  [xu,fu] = sem2d_extract_line(s,grid,0.0001,'h');
%   plot(xu,fu);
% Extracts field value at an arbitrary point, by SEM interpolation
% xz=[100,0];
% F = sem2d_extract_point(s11,grid,xz);


%Plot slip rate
%  mesh([0:data.nt-2]*data.dt,data.x/1e3,data.v);
%  xlabel('Time (s)')
%  ylabel('Along strike distance (km)')
%  zlabel('Slip rate (m/s)')



%    sem2d_velocity_plot(2500)
   figure(5)
   sem2d_snapshot_plot(s11,grid,[-inf,inf]);
   axis equal
   title('vx (m/s)')
