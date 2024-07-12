% This script is used for the evaluation of the Monte-Carlo data 
close all; clear;

%noise-free data:
ddsmpc = load('results_noisefree.mat');
%noisy data:
ddsmpc_noisy = load('results.mat');


%the following parameters are assumed to be equal for all methods:
load("MPC_init.mat", "params");
load("MPC_init.mat", "sys");
noise = load("noise_MC.mat");

sim_maxiter = size(ddsmpc.t_ocp,1);
sim_maxMCruns = size(ddsmpc.t_ocp,2);

%number of simulation steps per run:
t = 1:sim_maxiter;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  TIME PER ITERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%computation times for optimization: ddsmpc
tmean_ddsmpc = mean([ddsmpc.t_ocp],2);
tstd_ddsmpc = std([ddsmpc.t_ocp],1,2);
tminall_ddsmpc = min([ddsmpc.t_ocp],[],'all');
tmeanall_ddsmpc = mean([ddsmpc.t_ocp],'all');
tmaxall_ddsmpc = max([ddsmpc.t_ocp],[],'all');

%computation times for optimization: ddsmpc_noisy
tmean_ddsmpc_noisy = mean([ddsmpc_noisy.t_ocp],2);
tstd_ddsmpc_noisy = std([ddsmpc_noisy.t_ocp],1,2);
tminall_ddsmpc_noisy = min([ddsmpc_noisy.t_ocp],[],'all');
tmeanall_ddsmpc_noisy = mean([ddsmpc_noisy.t_ocp],'all');
tmaxall_ddsmpc_noisy = max([ddsmpc_noisy.t_ocp],[],'all');


figure;
%ddsmpc:
plot(t,1000*tmean_ddsmpc,'g','LineWidth',1.5)
hold on
set(gca,'Fontsize',12);
plot(t,1000*ones(size(t))*tmeanall_ddsmpc,'--g','LineWidth',1.5)
%ddsmpc_noisy:
plot(t,1000*tmean_ddsmpc_noisy,'b','LineWidth',1.5)
plot(t,1000*ones(size(t))*tmeanall_ddsmpc_noisy,'--b','LineWidth',1.5)
hold off
grid on
legend('noise-free','mean','noisy','mean', 'Interpreter', 'Latex', 'FontSize', 12)
xlim([t(2),t(end)])
%ylim([0,4])
xlabel('simulation step', 'Interpreter', 'Latex', 'FontSize', 16)
ylabel('time per iteration (ms)', 'Interpreter', 'Latex', 'FontSize', 16)
set(gcf,'units','normalized','outerposition',[0.25 0.25 0.5 0.5])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  TRAJECTORY COST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q = params.Q;
R = params.R;

cost_ddsmpc = zeros(sim_maxMCruns,sim_maxiter);
cost_ddsmpc_noisy = zeros(sim_maxMCruns,sim_maxiter);

for j=1:sim_maxMCruns
    %ddsmpc:
    for i = 1:sim_maxiter
        cost_ddsmpc(j,i) = (ddsmpc.X(:,i,j))'*Q*(ddsmpc.X(:,i,j)) + (ddsmpc.U(:,i,j))'*R*(ddsmpc.U(:,i,j));
    end
    
    %ddsmpc_noisy:
    for i = 1:sim_maxiter
        cost_ddsmpc_noisy(j,i) = (ddsmpc_noisy.X(:,i,j))'*Q*(ddsmpc_noisy.X(:,i,j)) + (ddsmpc_noisy.U(:,i,j))'*R*(ddsmpc_noisy.U(:,i,j));
    end
    
end

cost_ddsmpc_total = sum(cost_ddsmpc,2);
cost_ddsmpc_noisy_total = sum(cost_ddsmpc_noisy,2);
cost_ddsmpc_noisy_total = cost_ddsmpc_noisy_total/median(cost_ddsmpc_total);
cost_ddsmpc_total = cost_ddsmpc_total/median(cost_ddsmpc_total);

figure;
%subplot(1,2,1)
X = categorical([repmat("noise-free data",sim_maxMCruns,1);
                repmat("noisy data",sim_maxMCruns,1)],{'noise-free data','noisy data'});
Y = [cost_ddsmpc_total;
     cost_ddsmpc_noisy_total];
boxchart(X,Y,'MarkerStyle','none');
ylabel('Total trajectory cost','FontSize',14,'Interpreter','Latex')
%ylim([4000,11000])
set(gca,'TickLabelInterpreter','Latex')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  SCENARIO PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time axis of simulations:
tplot = 0:sys.T:sys.T*(sim_maxiter-1);

c1 = [0.2 0.5 0.84]; %color 1
c2 = [0.84 0.34 0.34]; %color 2
lw = 1.5; %LineWidth of state, input, reference
lw_c = 2; %Linewidth of constraints
c_c = [0.2 0.2 0.2]; %Color of constraints
fs = 14; %fontsize
%% plot all states:
%Plot states - noise-free data
f1 = figure('visible','off');
for k_plot = 1:10
    x = ddsmpc.X(:,:,k_plot);

    subplot(2,2,1)
    hold on
    p1 = plot(tplot,x(1,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
    p2 = plot(tplot,x(2,:),'LineWidth',lw,'Color',c2,'LineStyle','-.');


    subplot(2,2,3)
    hold on
    p3 = plot(tplot,x(3,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
    p4 = plot(tplot,x(4,:),'LineWidth',lw,'Color',c2,'LineStyle','-.');
end
subplot(2,2,1)
title('Noise-free data','FontSize',fs,'Interpreter','latex')
p5 = plot(tplot,zeros(size(tplot)),'-k','LineWidth',1.2); %Reference
ylim(ylim)
ylim("manual")
legend([p1,p2,p5(1)],{'$\theta_1$','$\theta_2$','Reference'},'AutoUpdate','Off', 'Interpreter', 'Latex','FontSize',10);
grid on;
yline(sys.g_X(1),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
yline(-sys.g_X(2),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
xlabel('time $t$ in s', 'Interpreter', 'Latex', 'FontSize', fs)
ylabel('Angles in rad', 'Interpreter', 'Latex', 'FontSize', fs)
xlim([0,sys.T*(sim_maxiter-1)]);
ylim([-0.25,1.6])
%set(gca,'YTick',-sys.g_X(2)-pi/2:pi/2:sys.g_X(1)+pi/2)
%set(gca,'YTickLabel',{'-2.5{\pi}','-2{\pi}','-1.5{\pi}','-{\pi}','-0.5{\pi}','0','0.5{\pi}','{\pi}','1.5{\pi}','2{\pi}','2.5{\pi}'})
grid on

subplot(2,2,3)
p5 = plot(tplot,zeros(size(tplot)),'-k','LineWidth',1.2); %Reference
legend([p3,p4,p5(1)],{'$\omega_1$','$\omega_2$','Reference'}, 'AutoUpdate','Off','Interpreter', 'Latex', 'Location', 'southeast','FontSize',10);
ylim(ylim)
ylim("manual") %These two commands make it such that all plots form now on are ignored for ylim
yline(sys.g_X(3),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
yline(-sys.g_X(4),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
xlabel('time $t$ in s', 'Interpreter', 'Latex', 'FontSize', fs)
ylabel('Velocities in rad/s', 'Interpreter', 'Latex', 'FontSize', fs)
xlim([0,sys.T*(sim_maxiter-1)]);
ylim([-1.8,0.8])
grid on

%set(f1,'visible','on')

%Plot states - noisy data
%f2 = figure('visible','off');
for k_plot = 1:10
    x = ddsmpc_noisy.X(:,:,k_plot);

    subplot(2,2,2)
    hold on
    p1 = plot(tplot,x(1,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
    p2 = plot(tplot,x(2,:),'LineWidth',lw,'Color',c2,'LineStyle','-.');


    subplot(2,2,4)
    hold on
    p3 = plot(tplot,x(3,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
    p4 = plot(tplot,x(4,:),'LineWidth',lw,'Color',c2,'LineStyle','-.');
end
subplot(2,2,2)
title('Noisy data','FontSize',fs,'Interpreter','latex')
p5 = plot(tplot,zeros(size(tplot)),'-k','LineWidth',1.2); %Reference
ylim(ylim)
ylim("manual")
legend([p1,p2,p5(1)],{'$\theta_1$','$\theta_2$','Reference'},'AutoUpdate','Off', 'Interpreter', 'Latex','FontSize',10);
grid on;
yline(sys.g_X(1),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
yline(-sys.g_X(2),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
xlabel('time $t$ in s', 'Interpreter', 'Latex', 'FontSize', fs)
ylabel('Angles in rad', 'Interpreter', 'Latex', 'FontSize', fs)
xlim([0,sys.T*(sim_maxiter-1)]);
ylim([-0.25,1.6])
%set(gca,'YTick',-sys.g_X(2)-pi/2:pi/2:sys.g_X(1)+pi/2)
%set(gca,'YTickLabel',{'-2.5{\pi}','-2{\pi}','-1.5{\pi}','-{\pi}','-0.5{\pi}','0','0.5{\pi}','{\pi}','1.5{\pi}','2{\pi}','2.5{\pi}'})
grid on

subplot(2,2,4)
p5 = plot(tplot,zeros(size(tplot)),'-k','LineWidth',1.2); %Reference
legend([p3,p4,p5(1)],{'$\omega_1$','$\omega_2$','Reference'}, 'AutoUpdate','Off','Interpreter', 'Latex', 'Location', 'southeast','FontSize',10);
ylim(ylim)
ylim("manual") %These two commands make it such that all plots form now on are ignored for ylim
yline(sys.g_X(3),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
yline(-sys.g_X(4),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
xlabel('time $t$ in s', 'Interpreter', 'Latex', 'FontSize', fs)
ylabel('Velocities in rad/s', 'Interpreter', 'Latex', 'FontSize', fs)
xlim([0,sys.T*(sim_maxiter-1)]);
ylim([-1.8,0.8])
grid on
set(gcf,'units','normalized','outerposition',[0.2 0.2 0.3 0.5])

set(f1,'visible','on')

%% plot states 1 and 3:
%Plot states - noise-free data
f1 = figure('visible','off');
for k_plot = 1:10
    x = ddsmpc.X(:,:,k_plot);

    subplot(2,2,1)
    hold on
    p1 = plot(tplot,x(1,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
    
    subplot(2,2,3)
    hold on
    p3 = plot(tplot,x(3,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
end
subplot(2,2,1)
title('Noise-free data','FontSize',fs,'Interpreter','latex')
p5 = plot(tplot,zeros(size(tplot)),'-k','LineWidth',1.2); %Reference
ylim(ylim)
ylim("manual")
%legend([p1,p2,p5(1)],{'$\theta_1$','Reference'},'AutoUpdate','Off', 'Interpreter', 'Latex');
grid on;
yline(sys.g_X(1),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
yline(-sys.g_X(2),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
xlabel('time $t$ in s', 'Interpreter', 'Latex', 'FontSize', fs)
ylabel('Angle $\theta_1$ in rad', 'Interpreter', 'Latex', 'FontSize', fs)
xlim([0,sys.T*(sim_maxiter-1)]);
ylim([-0.25,1.6])
%set(gca,'YTick',-sys.g_X(2)-pi/2:pi/2:sys.g_X(1)+pi/2)
%set(gca,'YTickLabel',{'-2.5{\pi}','-2{\pi}','-1.5{\pi}','-{\pi}','-0.5{\pi}','0','0.5{\pi}','{\pi}','1.5{\pi}','2{\pi}','2.5{\pi}'})
grid on

subplot(2,2,3)
p5 = plot(tplot,zeros(size(tplot)),'-k','LineWidth',1.2); %Reference
%legend([p3,p4,p5(1)],{'$\omega_1$','Reference'}, 'AutoUpdate','Off','Interpreter', 'Latex', 'Location', 'southeast');
ylim(ylim)
ylim("manual") %These two commands make it such that all plots form now on are ignored for ylim
yline(sys.g_X(3),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
yline(-sys.g_X(4),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
xlabel('time $t$ in s', 'Interpreter', 'Latex', 'FontSize', fs)
ylabel('Velocity $\omega_1$ in rad/s', 'Interpreter', 'Latex', 'FontSize', fs)
xlim([0,sys.T*(sim_maxiter-1)]);
ylim([-1.8,0.8])
grid on

%set(f1,'visible','on')

%Plot states - noisy data
%f2 = figure('visible','off');
for k_plot = 1:10
    x = ddsmpc_noisy.X(:,:,k_plot);

    subplot(2,2,2)
    hold on
    p1 = plot(tplot,x(1,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
    
    subplot(2,2,4)
    hold on
    p3 = plot(tplot,x(3,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
end
subplot(2,2,2)
title('Noisy data','FontSize',fs,'Interpreter','latex')
p5 = plot(tplot,zeros(size(tplot)),'-k','LineWidth',1.2); %Reference
ylim(ylim)
ylim("manual")
%legend([p1,p2,p5(1)],{'$\theta_1$','Reference'},'AutoUpdate','Off', 'Interpreter', 'Latex');
grid on;
yline(sys.g_X(1),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
yline(-sys.g_X(2),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
xlabel('time $t$ in s', 'Interpreter', 'Latex', 'FontSize', fs)
ylabel('Angle $\theta_1$ in rad', 'Interpreter', 'Latex', 'FontSize', fs)
xlim([0,sys.T*(sim_maxiter-1)]);
ylim([-0.25,1.6])
%set(gca,'YTick',-sys.g_X(2)-pi/2:pi/2:sys.g_X(1)+pi/2)
%set(gca,'YTickLabel',{'-2.5{\pi}','-2{\pi}','-1.5{\pi}','-{\pi}','-0.5{\pi}','0','0.5{\pi}','{\pi}','1.5{\pi}','2{\pi}','2.5{\pi}'})
grid on

subplot(2,2,4)
p5 = plot(tplot,zeros(size(tplot)),'-k','LineWidth',1.2); %Reference
%legend([p3,p4,p5(1)],{'$\omega_1$','Reference'}, 'AutoUpdate','Off','Interpreter', 'Latex', 'Location', 'southeast');
ylim(ylim)
ylim("manual") %These two commands make it such that all plots form now on are ignored for ylim
yline(sys.g_X(3),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
yline(-sys.g_X(4),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
xlabel('time $t$ in s', 'Interpreter', 'Latex', 'FontSize', fs)
ylabel('Velocity $\omega_1$ in rad/s', 'Interpreter', 'Latex', 'FontSize', fs)
xlim([0,sys.T*(sim_maxiter-1)]);
ylim([-1.8,0.8])
grid on
set(gcf,'units','normalized','outerposition',[0.2 0.2 0.3 0.5])

set(f1,'visible','on')

%% plot states 2 and 4:
%Plot states - noise-free data
f1 = figure('visible','off');
for k_plot = 1:10
    x = ddsmpc.X(:,:,k_plot);

    subplot(2,2,1)
    hold on
    p1 = plot(tplot,x(2,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
    
    subplot(2,2,3)
    hold on
    p3 = plot(tplot,x(4,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
end
subplot(2,2,1)
title('Noise-free data','FontSize',fs,'Interpreter','latex')
p5 = plot(tplot,zeros(size(tplot)),'-k','LineWidth',1.2); %Reference
ylim(ylim)
ylim("manual")
%legend([p1,p2,p5(1)],{'$\theta_2$','Reference'},'AutoUpdate','Off', 'Interpreter', 'Latex');
grid on;
yline(sys.g_X(1),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
yline(-sys.g_X(2),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
xlabel('time $t$ in s', 'Interpreter', 'Latex', 'FontSize', fs)
ylabel('Angle $\theta_2$ in rad', 'Interpreter', 'Latex', 'FontSize', fs)
xlim([0,sys.T*(sim_maxiter-1)]);
ylim([-0.25,1.6])
%set(gca,'YTick',-sys.g_X(2)-pi/2:pi/2:sys.g_X(1)+pi/2)
%set(gca,'YTickLabel',{'-2.5{\pi}','-2{\pi}','-1.5{\pi}','-{\pi}','-0.5{\pi}','0','0.5{\pi}','{\pi}','1.5{\pi}','2{\pi}','2.5{\pi}'})
grid on

subplot(2,2,3)
p5 = plot(tplot,zeros(size(tplot)),'-k','LineWidth',1.2); %Reference
%legend([p3,p4,p5(1)],{'$\omega_2$','Reference'}, 'AutoUpdate','Off','Interpreter', 'Latex', 'Location', 'southeast');
ylim(ylim)
ylim("manual") %These two commands make it such that all plots form now on are ignored for ylim
yline(sys.g_X(3),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
yline(-sys.g_X(4),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
xlabel('time $t$ in s', 'Interpreter', 'Latex', 'FontSize', fs)
ylabel('Velocity $\omega_2$ in rad/s', 'Interpreter', 'Latex', 'FontSize', fs)
xlim([0,sys.T*(sim_maxiter-1)]);
ylim([-1.8,0.8])
grid on

%set(f1,'visible','on')

%Plot states - noisy data
%f2 = figure('visible','off');
for k_plot = 1:10
    x = ddsmpc_noisy.X(:,:,k_plot);

    subplot(2,2,2)
    hold on
    p1 = plot(tplot,x(2,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
    
    subplot(2,2,4)
    hold on
    p3 = plot(tplot,x(4,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
end
subplot(2,2,2)
title('Noisy data','FontSize',fs,'Interpreter','latex')
p5 = plot(tplot,zeros(size(tplot)),'-k','LineWidth',1.2); %Reference
ylim(ylim)
ylim("manual")
%legend([p1,p2,p5(1)],{'$\theta_2$','Reference'},'AutoUpdate','Off', 'Interpreter', 'Latex');
grid on;
yline(sys.g_X(1),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
yline(-sys.g_X(2),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
xlabel('time $t$ in s', 'Interpreter', 'Latex', 'FontSize', fs)
ylabel('Angle $\theta_2$ in rad', 'Interpreter', 'Latex', 'FontSize', fs)
xlim([0,sys.T*(sim_maxiter-1)]);
ylim([-0.25,1.6])
%set(gca,'YTick',-sys.g_X(2)-pi/2:pi/2:sys.g_X(1)+pi/2)
%set(gca,'YTickLabel',{'-2.5{\pi}','-2{\pi}','-1.5{\pi}','-{\pi}','-0.5{\pi}','0','0.5{\pi}','{\pi}','1.5{\pi}','2{\pi}','2.5{\pi}'})
grid on

subplot(2,2,4)
p5 = plot(tplot,zeros(size(tplot)),'-k','LineWidth',1.2); %Reference
%legend([p3,p4,p5(1)],{'$\omega_2$','Reference'}, 'AutoUpdate','Off','Interpreter', 'Latex', 'Location', 'southeast');
ylim(ylim)
ylim("manual") %These two commands make it such that all plots form now on are ignored for ylim
yline(sys.g_X(3),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
yline(-sys.g_X(4),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
xlabel('time $t$ in s', 'Interpreter', 'Latex', 'FontSize', fs)
ylabel('Velocity $\omega_2$ in rad/s', 'Interpreter', 'Latex', 'FontSize', fs)
xlim([0,sys.T*(sim_maxiter-1)]);
ylim([-1.8,0.8])
grid on
set(gcf,'units','normalized','outerposition',[0.2 0.2 0.3 0.5])

set(f1,'visible','on')
