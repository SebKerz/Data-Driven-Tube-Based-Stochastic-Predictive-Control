%Run main_DDSMPC one before, to have all necessary variables

%Is origin feasible point for quadprog?
xk_measured = x0;
[A, b, Aeq, beq, lb, ub] = linconstraints(sys, data, constraints, controller.K, xk_measured);

alpha0 = zeros(size(mpc_init.H,1),1);

if all(A*alpha0<=b)
    disp("Inequality consraint okay")
else
    idx_problem = find(A*alpha0>b);
    disp("Inequality consraint does not hold for xualpha0, problem with constraint "+num2str(idx_problem))
    %
    for k = 1:length(idx_problem)
        xualpha_idx_problem = find(A(idx_problem(k),:));
        disp("Problem "+num2str(idx_problem(k))+" concerns element "+num2str(xualpha_idx_problem)+" of xualpha")
    end
        disp("Up to "+num2str(sys.n*(params.L+1))+" are state constraints, up to "+num2str(sys.n*(params.L+1)+sys.m*(params.L))+" are input constraints")
end
%% Check First step Constraint
%Show allowed 1st steps
figure
sys.X.plot()
hold on
constraints.polyfs.plot('Color','b','alpha',1)
constraints.polyinit.plot('Color','g','alpha',0.1)
constraints.Xf.plot('Color', 'yellow','alpha',1)
xlabel('x1')
ylabel('x2')
title("Constraint sets")
legend('X','First Step Constraint','Feasible Initial set','Terminal set')

% %Input set
%BU = (-sys.B)*sys.input.U;
% figure
% subplot(2,1,1)
% BU.projection(1:2).plot()
% xlabel('x1')
% ylabel('x2')
% title("Input set BU")
% subplot(2,1,2)
% BU.projection(3:4).plot()
% xlabel('x3')
% ylabel('x4')
% title("Input set BU")
% 
% %Minus
%minusInput = constraints.polyfs+(BU);
% figure
% subplot(2,1,1)
% minusInput.projection(1:2).plot()
% xlabel('x1')
% ylabel('x2')
% title("FS minus InputSet")
% subplot(2,1,2)
% minusInput.projection(3:4).plot()
% xlabel('x3')
% ylabel('x4')
% title("FS minus InputSet")

%Check where we would need to start to reach that
%constraints.polyinit = inv(sys.A)*minusInput;
figure
subplot(2,1,1)
constraints.polyinit.plot()
xlabel('x1')
ylabel('x2')
title("Initial states from which FirstStep-Set reachable, projected")
subplot(2,1,2)
constraints.polyinit.projection(3:4).plot()
xlabel('x3')
ylabel('x4')
title("Initial states from which FirstStep-Set reachable, projected")

clear BU minusInput

%% Check Terminal set
figure
subplot(2,1,1)
constraints.Xf.projection(1:2).plot()
xlabel('x1')
ylabel('x2')
title("Terminal set, projected")
subplot(2,1,2)
constraints.Xf.projection(3:4).plot()
xlabel('x3')
ylabel('x4')
title("Terminal set, projected")