%Run main_init first.
%Alternatively, load the following structs: data, sys, controller,
load("MPC_init.mat")
%noise = load("noise_MC.mat");

rng('shuffle'); % random seed
s = rng;

params.lambda = 0; %Regularization of alpha

x0 = [2;3];%[2.5;3];
%x0 = [0;0;0;0]; %Initial state

%Specify state reference x_ref_user, length also specifies length of sim
%x_ref_user = [zeros(sys.n,10),pi*[ones(2,70);zeros(2,70)],zeros(sys.n,20)];
x_ref_user = zeros(sys.n, 20);

u_ref = zeros(params.L+1,1); %Input reference value is set to zero

N_runs = 500; %Number of MPC Runs

plot_predictions = 0; %Specify wether to plot predictions in each step, only works for n=2

options = optimset('Display','off', 'Algorithm', 'active-set','TolCon',1e-12);

%Get BigQ, BigR, H etc.
mpc_init = init_mpc(sys,params, controller, data, x_ref_user);

X = zeros(sys.n,mpc_init.sim_length,N_runs); U = zeros(sys.m,mpc_init.sim_length,N_runs); XH = zeros(sys.n,mpc_init.sim_length,N_runs);
D = zeros(sys.md,mpc_init.sim_length,N_runs); t_ocp = zeros(mpc_init.sim_length,N_runs); V_ocp = zeros(mpc_init.sim_length,N_runs);
alpha_norm = zeros(mpc_init.sim_length,N_runs);

for i = 1:N_runs
    disp("MPC run "+num2str(i))

    %Generate measurement noise and disturbance sequences
    eps_sequence = 1*random(sys.noise.eps_pdf,[sys.n,mpc_init.sim_length]);
    d_sequence = 1*random(sys.disturbance.d_pdf,[sys.md,mpc_init.sim_length]);
    %eps_sequence = noise.eps_MC{i};
    %d_sequence = noise.d_MC{i};

    for k = 1:mpc_init.sim_length
        if k == 1
            xk = x0;
        end
        xk_measured = xk + eps_sequence(:,k); %Measured state

        X(:,k,i) = xk; %Store
        XH(:,k,i) = xk_measured; %Store

        %Prepare quadprog OCP
        [A, b, Aeq, beq, lb, ub] = linconstraints(sys, data, constraints, controller.K, xk_measured);
        f = [-mpc_init.BigQ*mpc_init.x_ref(sys.n*(k-1)+(1:sys.n*(params.L+1))) - mpc_init.BigK'*mpc_init.BigR*u_ref;...
            - mpc_init.BigR*u_ref; zeros(size(data.Hx,2),1)];

        %Warmstart
        if k>1
            xshift = xualpha(1:sys.n*(params.L+1));
            ushift = xualpha(sys.n*(params.L+1)+(1:sys.m*(params.L+1))); %xualpha(sys.n*(params.L+1)+(sys.m*params.L+1:sys.m*(params.L+1)));
            alphashift = pinv([data.Hu-mpc_init.KHx;data.Hd;data.Hx])*[ushift;zeros(sys.md*(params.L+1),1);xshift];
            alpha0 = [xshift; ushift; alphashift];
        else
            alpha0 = zeros(size(mpc_init.H,1),1);
        end

        %Solve OCP:
        tic
        [xualpha, V, exitflag] = quadprog(mpc_init.H, f, A, b, Aeq, beq, lb, ub, alpha0, options);
        t_ocp(k,i) = toc; %Store elapsed time

        %Check Feasibility, stop if not feasible
        if(exitflag == -2)
            %Infeasible, only use state feedback
            disp("No feasible solution found in time step "+num2str(k));
            V = -1;
            %u = ushift(2)+controller.K*xk_measured;
        else
            if plot_predictions
                figure(77)
                xpred = xualpha(1:sys.n*(params.L+1));
                xpred = reshape(xpred,[sys.n, params.L+1]);
                plot(xpred(1,1),xpred(2,1),'ok')
                plot(xpred(1,2:end),xpred(2,2:end),'*k')
            end

            
        end
        u =  xualpha(sys.n*(params.L+1)+(1:sys.m)) + controller.K*xk_measured;
        if abs(u)>sys.input.u_max
            disp("Input constraints violated in time step "+num2str(k)+", out of bounds by "+num2str(abs(u)-sys.input.u_max))
            u = sys.input.u_max*sign(u);
        end

        d = d_sequence(:,k);

        %Store:
        V_ocp(k,i) = V;
        U(:,k,i) = u;
        D(:,k,i) = d;
        alpha_norm(k,i) = norm(xualpha(end-size(data.Hx,2)+1:end));

        %Apply input and disturbance to system:
        xk = sys.A*xk+sys.B*u+sys.Bd*d;
    end
end

disp(['Mean computation time of quadprog: ', num2str(mean(t_ocp,'all')*1e3), ' ms'])
disp([num2str(length(find(V_ocp==-1))),' infeasibilities, ',num2str(100*length(find(V_ocp==-1))/(mpc_init.sim_length*N_runs)),'% of all OCPs.'])
%disp(['Constraints satisfied for ',num2str(100*sum(sys.X.contains(X))/(size(X,2)*size(X,3))),'% of steps, risk param was ',num2str(params.p)])
%% Plot
fs = 12; %fontsize

disp(['Constraints satisfied for ',num2str(100*sum(sys.X.contains(X))/(size(X,2)*size(X,3))),'% of total steps, risk param was ',num2str(params.p)])

if sys.n==4
    tplot=0:sys.T:(mpc_init.sim_length-1)*sys.T;

    c1 = [0.2 0.5 0.84]; %color 1
    c2 = [0.84 0.34 0.34]; %color 2
    lw = 1; %LineWidth of state, input, reference
    lw_c = 2; %Linewidth of constraints
    c_c = [0.2 0.2 0.2]; %Color of constraints


    %Plot states
    f1 = figure('visible','off');
    for k_plot = 1:N_runs
        x = X(:,:,k_plot);

        subplot(2,1,1)
        hold on
        p1 = plot(tplot,x(1,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
        p2 = plot(tplot,x(2,:),'LineWidth',lw,'Color',c2,'LineStyle','--');


        subplot(2,1,2)
        hold on
        p3 = plot(tplot,x(3,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
        p4 = plot(tplot,x(4,:),'LineWidth',lw,'Color',c2,'LineStyle','--');
    end
    subplot(2,1,1)
    p5 = plot(tplot,x_ref_user(1:2,:)','-k','LineWidth',1.2); %Reference
    ylim(ylim)
    ylim("manual")
    legend([p1,p2,p5(1)],{'$\theta_1$ in rad','$\theta_2$ in rad','Reference for $\theta_1$, $\theta_2$'},'AutoUpdate','Off', 'Interpreter', 'Latex');
    grid on;
    yline(sys.g_X(1),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
    yline(-sys.g_X(2),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
    xlabel('time $t$ (s)', 'Interpreter', 'Latex', 'FontSize', fs)
    ylabel('Angles', 'Interpreter', 'Latex', 'FontSize', fs)
    xlim([0,sys.T*(mpc_init.sim_length-1)]);
    %ylim([min([min(min(x_ref_user(1:2,:))),x0(1:2)'])-0.25,max([max(max(x_ref_user(1:2,:))),x0(1:2)'])+0.25])
    %set(gca,'YTick',-sys.g_X(2)-pi/2:pi/2:sys.g_X(1)+pi/2)
    %set(gca,'YTickLabel',{'-2.5{\pi}','-2{\pi}','-1.5{\pi}','-{\pi}','-0.5{\pi}','0','0.5{\pi}','{\pi}','1.5{\pi}','2{\pi}','2.5{\pi}'})
    set(gcf,'units','normalized','outerposition',[0.2 0.2 0.5 0.5])
    grid on

    subplot(2,1,2)
    p5 = plot(tplot,x_ref_user(3:4,:)','-k','LineWidth',1.2); %Reference
    legend([p3,p4,p5(1)],{'$\omega_1$ in rad/s','$\omega_2$ in rad/s','Reference for $\omega_1$, $\omega_2$'}, 'AutoUpdate','Off','Interpreter', 'Latex');
    ylim(ylim)
    ylim("manual") %These two commands make it such that all plots form now on are ignored for ylim
    yline(sys.g_X(3),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
    yline(-sys.g_X(4),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
    xlabel('time $t$ (s)', 'Interpreter', 'Latex', 'FontSize', fs)
    ylabel('Velocities', 'Interpreter', 'Latex', 'FontSize', fs)
    xlim([0,sys.T*(mpc_init.sim_length-1)]);
    %ylim([min([min(min(x(3:4,:))),x0(3:4)'])-0.25,max([max(max(x(3:4,:))),x0(3:4)'])+0.25])
    grid on
    set(gcf,'units','normalized','outerposition',[0.2 0.2 0.4 0.8])

    set(f1,'visible','on')

    %Plot input and disturbance
    f2 = figure('visible','off');
    hold on;
    for k_plot = 1:N_runs
        u = U(:,:,k_plot);
        d = D(:,:,k_plot);
        plot(tplot,u','LineWidth',lw,'Color',c1,'LineStyle','-')
        plot(tplot,d','LineWidth',lw,'Color',c2,'LineStyle','--')
    end
    grid on;
    yline(sys.input.u_max,'LineStyle',':','Color',c1,'LineWidth',lw_c);
    yline(-sys.input.u_max,'LineStyle',':','Color',c1,'LineWidth',lw_c);
    yline(params.d_max,'LineStyle',':','Color',c2,'LineWidth',lw_c);
    yline(-params.d_max,'LineStyle',':','Color',c2,'LineWidth',lw_c);
    legend('Input $M_M$','Disturbance $M_d$', 'AutoUpdate','Off','Interpreter', 'Latex');
    xlabel('time $t$ (s)', 'Interpreter', 'Latex', 'FontSize', fs)
    ylabel('Torque (Nm)', 'Interpreter', 'Latex', 'FontSize', fs)
    xlim([0,sys.T*(mpc_init.sim_length-1)]);
    set(gcf,'units','normalized','outerposition',[0.2 0.2 0.5 0.5])
    ylim([-sys.input.u_max-0.1, sys.input.u_max+0.1])
    set(f2,'visible','on')

    %Plot measured states
    f3 = figure('visible','off');
    for k_plot = 1:N_runs
        xh = XH(:,:,k_plot);

        subplot(2,1,1)
        hold on
        p1 = plot(tplot,xh(1,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
        p2 = plot(tplot,xh(2,:),'LineWidth',lw,'Color',c2,'LineStyle','--');


        subplot(2,1,2)
        hold on
        p3 = plot(tplot,xh(3,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
        p4 = plot(tplot,xh(4,:),'LineWidth',lw,'Color',c2,'LineStyle','--');
    end
    subplot(2,1,1)
    p5 = plot(tplot,x_ref_user(1:2,:)','--k','LineWidth',1.2); %Reference
    legend([p1,p2,p5(1)],{'measured $\theta_1$ in rad','measured $\theta_2$ in rad','Reference for $\theta_1$, $\theta_2$'},'AutoUpdate','Off', 'Interpreter', 'Latex');
    ylim(ylim)
    ylim("manual")
    grid on;
    yline(sys.g_X(1),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
    yline(-sys.g_X(2),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
    xlabel('time $t$ (s)', 'Interpreter', 'Latex', 'FontSize', fs)
    ylabel('Angles', 'Interpreter', 'Latex', 'FontSize', fs)
    xlim([0,sys.T*(mpc_init.sim_length-1)]);
    %ylim([min([min(min(x_ref_user(1:2,:))),x0(1:2)'])-0.25,max([max(max(x_ref_user(1:2,:))),x0(1:2)'])+0.25])
    %set(gca,'YTick',-sys.g_X(2)-pi/2:pi/2:sys.g_X(1)+pi/2)
    %set(gca,'YTickLabel',{'-2.5{\pi}','-2{\pi}','-1.5{\pi}','-{\pi}','-0.5{\pi}','0','0.5{\pi}','{\pi}','1.5{\pi}','2{\pi}','2.5{\pi}'})
    set(gcf,'units','normalized','outerposition',[0.2 0.2 0.5 0.5])
    grid on

    subplot(2,1,2)
    p5 = plot(tplot,x_ref_user(3:4,:)','--k','LineWidth',1.2); %Reference
    legend([p3,p4,p5(1)],{'measured $\omega_1$ in rad/s','measured $\omega_2$ in rad/s','Reference for $\omega_1$, $\omega_2$'}, 'AutoUpdate','Off','Interpreter', 'Latex');
    ylim(ylim)
    ylim("manual")
    yline(sys.g_X(3),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
    yline(-sys.g_X(4),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
    xlabel('time $t$ (s)', 'Interpreter', 'Latex', 'FontSize', fs)
    ylabel('Velocities', 'Interpreter', 'Latex', 'FontSize', fs)
    xlim([0,sys.T*(mpc_init.sim_length-1)]);
    grid on
    set(gcf,'units','normalized','outerposition',[0.2 0.2 0.4 0.8])
    %ylim([min([min(min(x(3:4,:))),x0(3:4)'])-0.25,max([max(max(x(3:4,:))),x0(3:4)'])+0.25])
    set(f3,'visible','on')
end

if sys.n==2
    tplot=0:sys.T:(mpc_init.sim_length-1)*sys.T;

    c1 = [0.2 0.5 0.84]; %color 1
    c2 = [0.84 0.34 0.34]; %color 2
    lw = 1; %LineWidth of state, input, reference
    lw_c = 2; %Linewidth of constraints
    c_c = [0.2 0.2 0.2]; %Color of constraints
    fs = 11; %fontsize

    %Plot states
    f1 = figure('visible','off');
    for k_plot = 1:N_runs
        x = X(:,:,k_plot);
        hold on
        p1 = plot(tplot,x(1,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
        p2 = plot(tplot,x(2,:),'LineWidth',lw,'Color',c2,'LineStyle','--');
    end
    p5 = yline(sys.g_X(1),'LineStyle','-','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
    plot(tplot,x_ref_user(1:2,:)',':k','LineWidth',1.2); %Reference
    ylim(ylim)
    ylim("manual")
    legend([p1,p2,p5(1)],{'$\theta_1$ in rad','$\theta_2$ in rad','Constraint on $\theta_1$'},'AutoUpdate','Off', 'Interpreter', 'Latex');
    grid on;
    yline(sys.g_X(1),'LineStyle','-','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
    yline(-sys.g_X(2),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
    xlabel('time $t$ (s)', 'Interpreter', 'Latex', 'FontSize', fs)
    ylabel('Angles', 'Interpreter', 'Latex', 'FontSize', fs)
    xlim([0,sys.T*(mpc_init.sim_length-1)]);
    %ylim([min([min(min(x_ref_user(1:2,:))),x0(1:2)'])-0.25,max([max(max(x_ref_user(1:2,:))),x0(1:2)'])+0.25])
    %set(gca,'YTick',-sys.g_X(2)-pi/2:pi/2:sys.g_X(1)+pi/2)
    %set(gca,'YTickLabel',{'-2.5{\pi}','-2{\pi}','-1.5{\pi}','-{\pi}','-0.5{\pi}','0','0.5{\pi}','{\pi}','1.5{\pi}','2{\pi}','2.5{\pi}'})
    set(gcf,'units','normalized','outerposition',[0.2 0.2 0.5 0.5])
    grid on
    set(f1,'visible','on')

    %Plot input and disturbance
    f2 = figure('visible','off');
    hold on;
    for k_plot = 1:N_runs
        u = U(:,:,k_plot);
        d = D(:,:,k_plot);
        plot(tplot,u','LineWidth',lw,'Color',c1,'LineStyle','-')
        plot(tplot,d','LineWidth',lw,'Color',c2,'LineStyle','--')
    end
    grid on;
    yline(sys.input.u_max,'LineStyle',':','Color',c1,'LineWidth',lw_c);
    yline(-sys.input.u_max,'LineStyle',':','Color',c1,'LineWidth',lw_c);
    yline(params.d_max,'LineStyle',':','Color',c2,'LineWidth',lw_c);
    yline(-params.d_max,'LineStyle',':','Color',c2,'LineWidth',lw_c);
    legend('Input $M_M$','Disturbance $M_d$', 'AutoUpdate','Off','Interpreter', 'Latex');
    xlabel('time $t$ (s)', 'Interpreter', 'Latex', 'FontSize', fs)
    ylabel('Torque (Nm)', 'Interpreter', 'Latex', 'FontSize', fs)
    xlim([0,sys.T*(mpc_init.sim_length-1)]);
    set(gcf,'units','normalized','outerposition',[0.2 0.2 0.5 0.5])
    ylim([-sys.input.u_max-0.1, sys.input.u_max+0.1])
    set(f2,'visible','on')

    %Plot measured states
    f3 = figure('visible','off');
    for k_plot = 1:N_runs
        xh = XH(:,:,k_plot);
        hold on
        p1 = plot(tplot,xh(1,:),'LineWidth',lw,'Color',c1,'LineStyle','-');
        p2 = plot(tplot,xh(2,:),'LineWidth',lw,'Color',c2,'LineStyle','--');
    end
    p5 = plot(tplot,x_ref_user(1:2,:)','--k','LineWidth',1.2); %Reference
    legend([p1,p2,p5(1)],{'measured $\theta_1$ in rad','measured $\theta_2$ in rad','Reference for $\theta_1$, $\theta_2$'},'AutoUpdate','Off', 'Interpreter', 'Latex');
    ylim(ylim)
    ylim("manual")
    grid on;
    yline(sys.g_X(1),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c); %Constraints
    yline(-sys.g_X(2),'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',lw_c);
    xlabel('time $t$ (s)', 'Interpreter', 'Latex', 'FontSize', fs)
    ylabel('Angles', 'Interpreter', 'Latex', 'FontSize', fs)
    xlim([0,sys.T*(mpc_init.sim_length-1)]);
    %ylim([min([min(min(x_ref_user(1:2,:))),x0(1:2)'])-0.25,max([max(max(x_ref_user(1:2,:))),x0(1:2)'])+0.25])
    %set(gca,'YTick',-sys.g_X(2)-pi/2:pi/2:sys.g_X(1)+pi/2)
    %set(gca,'YTickLabel',{'-2.5{\pi}','-2{\pi}','-1.5{\pi}','-{\pi}','-0.5{\pi}','0','0.5{\pi}','{\pi}','1.5{\pi}','2{\pi}','2.5{\pi}'})
    set(gcf,'units','normalized','outerposition',[0.2 0.2 0.5 0.5])
    grid on
    set(f3,'visible','on')

    %Plot in state space
    figure(77)
    set(gcf,'Position',[ 680   589   543   389])
    s1 = subplot(2,2,1,'Position',[0.085 0.56 0.4 0.39]);
    %sys.X.plot('color','white','alpha',0.2)
    hold on
    % if isfield(constraints,'polyfs')
    %     constraints.polyfs.plot('Color','b','alpha',.1)
    %     constraints.polyinit.plot('Color','g','alpha',0.1)
    % end
    %constraints.Xf.plot('Color', 'yellow','alpha',.1)
    
    %title("Constraint sets")
    for k_plot = 1:N_runs
        x = X(:,:,k_plot);
        plot(x(1,1:end),x(2,1:end),'Color',c1,'LineWidth',0.2)
    end
    XsetPlot = sys.X.plot('color','white','alpha',0);
    l1 = legend(XsetPlot,'$\mathcal{X}$','AutoUpdate','off','FontSize', fs,'Interpreter','latex');%,'First Step Constraint','Feasible Initial set','Terminal set','AutoUpdate','off')
    set(l1,'Position',[0.33,0.62,0.115,0.059])
    set(gca,'FontSize',fs')
    xlabel('$x_1$','FontSize', fs+2,'Interpreter','latex')
    ylabel('$x_2$','FontSize', fs+2,'Interpreter','latex')
    plot(x0(1),x0(2),'*k','Linewidth',3)
    s2 = subplot(2,2,2,'Position',[0.565 0.56 0.38 0.39]);
    hold on
    % if isfield(constraints,'polyfs')
    %     constraints.polyfs.plot('Color','b','alpha',.1)
    %     constraints.polyinit.plot('Color','g','alpha',0.1)
    % end
    %constraints.Xf.plot('Color', 'yellow','alpha',.1)
    %legend('X','AutoUpdate','off','Location','Southeast','FontSize', fs);%,'First Step Constraint','Feasible Initial set','Terminal set','AutoUpdate','off')
    
    %text(x0(1)+0.05,x0(2)+0.05,'Start','FontSize',fs)
    %title("Constraint sets")
    for k_plot = 1:N_runs
        x = X(:,:,k_plot);
        plot(x(1,1:end),x(2,1:end),'Color',c1,'LineWidth',0.2)
    end
    sys.X.plot('color','white','alpha',0)
    xlim([1.5 2.2])
    ylim([0 3.1])
    set(gca,'FontSize',fs)
    xlabel('$x_1$','FontSize', fs+2,'Interpreter','latex')
    ylabel('$x_2$','FontSize', fs+2,'Interpreter','latex')
    plot(x0(1),x0(2),'*k','Linewidth',3)
    s3 = subplot(2,2,[3 4],'Position',[0.15 0.125 0.75 0.27]);
    constViolation = reshape(sys.X.contains(X),[mpc_init.sim_length,N_runs])';
    conSat_perStep = sum(constViolation)/N_runs;
    plot(100*conSat_perStep,'linewidth',lw)
    %title("Empirical probability of constraint violation")
    set(gca,'FontSize',fs)
    xlabel('step','FontSize', fs+2,'Interpreter','latex')
    ylabel('$\textnormal{Pr}(\bm{x}\in\mathcal{X})$','FontSize', fs+2,'Interpreter','latex')
    ylim([75 101])
    hold on 
    plot([0, 25], [params.p*100, params.p*100],'Color','k','LineWidth',0.8,'LineStyle','--')
    grid on
end
