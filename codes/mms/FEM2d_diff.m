% FEM code for order of verification test in 2D diffusion equation
% Modified from the template code as a part of the course "Introduction to
% Finite element methods in geosciences", Spring 2019, ETH Zurich.
% (c) Sthavishtha, 2019
% -------------------------------------------------------------------------
% Features
% - implemented bilinear & biquadratic shape functions 
% - 2x2 & 3x3 Gauss-legendre quadrature numerical integration schemes
% - just comment/uncomment the necessary part to transfrom from bilinear to
%   biquadratic shape functions
% -------------------------------------------------------------------------

%% GENERAL STUFF
    clear                                                                   % clear the current Workspace
    close all                                                               % close all figure windows
    clc                                                                     % clear the Command Window
    
% GEOMETRICAL PARAMETERS
    lx          =	01;                                                     % length of model [m]
    ly          =   01;                                                     % computation time    
    
% PHYSICAL PARAMETERS
    kappa_x     =	1;                                                      % thermal diffusivity components [m2/s]
    kappa_y     =	1;                                  
    kappa       =	[kappa_x 0 ; 0 kappa_y];                                % thermal diffusivity matrix [m2/s]
    
% NUMERICAL PARAMETERS
    el_x        =   16;                                                     % #elements in x direction
    el_y        =   16;
    ndof        =   2;
    el_tot      =	el_x*el_y;                                              % #elements total
    n_x         =	ndof*el_x + 1;                                          % total #nodes along x direction
    n_y         =	ndof*el_y + 1;                                          % total #nodes along y direction
    n_per_el    =   9;                                                      % #nodes per element (bilinear/biquadratic)
    n_tot       =	n_x*n_y;                                                % #nodes total
    
% QUADRATURE PARAMETERS
    % 2 point GLg
    n_int       =   4;                                                      % #integration points  
    int_pts   	=   [-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3) ; ...        % integration points in [(-1,-1) ; (1,1)]
                    -1/sqrt(3)  1/sqrt(3) 1/sqrt(3) -1/sqrt(3)];                        
    weight      =   [1  1  1  1];                                           % weights of integration points  

    % 3 point GLg
%     n_int      =   9;                                                       % #integration points
%     int_pts    =   [-0.774596669241483  -0.774596669241483 -0.774596669241483 ...
%         0.0 0.0 0.0 ...
%         0.774596669241483 0.774596669241483 0.774596669241483 ;
%                     -0.774596669241483  0.0 0.774596669241483 ...
%         -0.774596669241483  0.0 0.774596669241483 ...
%         -0.774596669241483 0.0 0.774596669241483];                          % integration points in [-1,1]
%     weight     =   [0.3086 0.4938 0.3086 0.4938 0.7901 0.4938 0.3086 0.4938 ...
%         0.3086];                                                            % weights of integration points
                                                                         
% CREATE NUMERICAL GRID
    dx          =   lx/(n_x - 1);                                           % distance between two nodes
    dy          =   ly/(n_y - 1);        
    GCOORDS     =   zeros(2,n_x*n_y);                                       % global node coordinates
    GLINDEX     =   zeros(n_x,n_y);                                         % global node index   
    XC          =   zeros(n_x,n_y);                                         % x,y coordinates in 2D form
    YC          =   zeros(n_x,n_y);    
    deform      =   0.0;                                                    % perturb the global coords for debugging

% COMPUTING & PERTURBING GLOBAL COORDINATES 
    for i = 1: n_x
        for j = 1 : n_y
            xn = (i-1)*dx;
            yn = (j-1)*dy;
            nid = (i-1)+(j-1)*n_x + 1;
            GCOORDS(1,nid) = xn;
            GCOORDS(2,nid) = yn;
            
            if i == n_x/2 
                GCOORDS(1,nid) = xn + deform;
                if j == n_y/2                  
                    GCOORDS(2,nid) = yn + deform;
                end
                if j == n_y/2 - 1
                    GCOORDS(2,nid) = yn - deform;
                end
                if j == n_y/2 + 1
                    GCOORDS(2,nid) = yn + deform;
                end                
            end            
            GLINDEX(i,j) = nid;
        end
    end
        
% LOCAL-TO-GLOBAL MAPPING  
    EL_N = zeros(n_per_el,el_tot);
    el = 0;                                                                 % relates local to global node numbers per element       
    for i = 1 : ndof : n_y - ndof
        for j = 1 : ndof : n_x - ndof        
            el = el + 1;
            EL_N(1,el)   =   GLINDEX(j,i);
            EL_N(2,el)   =   GLINDEX(j,i + 1);  
            EL_N(3,el)   =   GLINDEX(j,i + 2);
            EL_N(4,el)   =   GLINDEX(j + 1,i + 2);
            EL_N(5,el)   =   GLINDEX(j + 1,i + 1);  
            EL_N(6,el)   =   GLINDEX(j + 1,i);  
            EL_N(7,el)   =   GLINDEX(j + 2,i);
            EL_N(8,el)   =   GLINDEX(j + 2,i + 1);  
            EL_N(9,el)   =   GLINDEX(j + 2,i + 2);
%             EL_N(1,el)   =   GLINDEX(j,i);
%             EL_N(2,el)   =   GLINDEX(j,i + 1);  
%             EL_N(3,el)   =   GLINDEX(j + 1,i + 1);
%             EL_N(4,el)   =   GLINDEX(j + 1,i);
        end
    end
	
% BOUNDARY CONDITIONS   
    i_b = 0;
    i_t = 0;
    i_l = 0;
    i_r = 0;    
    for i = 1:n_x
        for j = 1:n_y
            if GCOORDS(2,GLINDEX(i,j)) == 0 % bottom wall
                i_b = i_b + 1;
                bc_dof_b(i_b) = GLINDEX(i,j);
                bc_val_b(i_b) = MMS_Solution(GCOORDS(1,GLINDEX(i,j)),...
                    GCOORDS(2,GLINDEX(i,j)),lx,ly);
            elseif GCOORDS(2,GLINDEX(i,j)) == ly % top wall
                i_t = i_t + 1;
                bc_dof_t(i_t) = GLINDEX(i,j);
                bc_val_t(i_t) = MMS_Solution(GCOORDS(1,GLINDEX(i,j)),...
                    GCOORDS(2,GLINDEX(i,j)),lx,ly);
            elseif GCOORDS(1,GLINDEX(i,j)) == 0 % left wall
                i_l = i_l + 1;
                bc_dof_l(i_l) = GLINDEX(i,j);
                bc_val_l(i_l) = MMS_Solution(GCOORDS(1,GLINDEX(i,j)),...
                    GCOORDS(2,GLINDEX(i,j)),lx,ly);                
            elseif GCOORDS(1,GLINDEX(i,j)) == lx % right wall
                i_r = i_r + 1;
                bc_dof_r(i_r) = GLINDEX(i,j);
                bc_val_r(i_r) = MMS_Solution(GCOORDS(1,GLINDEX(i,j)),...
                    GCOORDS(2,GLINDEX(i,j)),lx,ly);                  
            end
        end
    end

    bc_dof      =     [bc_dof_b bc_dof_t bc_dof_l bc_dof_r];
    bc_val      =     [bc_val_b bc_val_t bc_val_l bc_val_r];        
    
% INITIALIZATION OF ALL KINDS OF STUFF    
    T_1         =   MMS_Solution(GCOORDS(1,:),GCOORDS(2,:),lx,ly);          % initial condition for temperature                                        
    T           =   reshape(T_1,n_tot,1);
    T_ana       =   zeros(n_x,n_y);                                         % analytical MMS solution
    T2d         =   zeros(n_x,n_y);                                         % temperature in 2D mesh grid format    
    KG          =	zeros(n_tot,n_tot);                                     % global stiffness matrix
    FG          =	zeros(n_tot,1);                                         % global force vector      
    Nloc        =   zeros(n_int,n_per_el);                                  % local shape function matrix
    dNloc       =   zeros(2,n_per_el,n_int);                                % local shape function derivative matrix
  
    for i_int = 1 : n_int
        Nloc(i_int,:)      =   Nbiquad2D(int_pts(1,i_int),int_pts(2,i_int));
        dNloc(:,:,i_int)   =   dNbiquad2D(int_pts(1,i_int),int_pts(2,i_int));    
%         Nloc(i_int,:)      =   Nbil2D(int_pts(1,i_int),int_pts(2,i_int));
%         dNloc(:,:,i_int)   =   dNbil2D(int_pts(1,i_int),int_pts(2,i_int));     
    end       		
           
    for i = 1: n_x
        for j = 1 : n_y
            nid     =   (i-1)+(j-1)*n_x + 1;
            XC(i,j) =   GCOORDS(1,nid);
            YC(i,j) =   GCOORDS(2,nid);
            T2d(i,j)=   T(nid);
        end
    end            

%% LOOPING
    for iel = 1:el_tot % ELEMENT LOOP
        n_now   =   EL_N(:,iel);                                            % which nodes are in the current element?
        
        Mloc    =   zeros(n_per_el,n_per_el);
        Kloc    =   zeros(n_per_el,n_per_el);                               % local stiffness matrix
        Floc    =   zeros(n_per_el,1);                                      % local force vector  

        for i_int = 1:n_int        
            [jacobdN,invjacobdN,det_jacobdN] = jacobDN_compute(dNloc(:,:,i_int), ...
                GCOORDS(:,n_now(:))');
            dNlocT  =   invjacobdN*dNloc(:,:,i_int);                        % derivative of shape fn. wrt global coords.
            Kloc    =   Kloc + dNlocT'*kappa*dNlocT*weight(i_int)*det_jacobdN;  
            x_int   =   Nloc(i_int,:)*GCOORDS(1,n_now)';                    % integration points in global coords.
            y_int   =   Nloc(i_int,:)*GCOORDS(2,n_now)';
            source  =   MMS_RHS(x_int,y_int,lx,ly);                         % source computed from RHS of MMS
            Floc    =   Floc + source*Nloc(i_int,:)'*weight(i_int)*det_jacobdN;         
        end

        % add local matrices to global ones
        KG(n_now,n_now)     =   KG(n_now,n_now) + Kloc;   
        FG(n_now)           =   FG(n_now) + Floc;
    end % ELEMENT LOOP        
    
    % APPlY BOUNDARY CONDITIONS
    for i = 1:length(bc_dof)
        KG(bc_dof(i),:)        =   0;                                       % set bc-rows to zero
        KG(bc_dof(i),bc_dof(i))=   1; 
        FG(bc_dof(i))          =   bc_val(i); 
    end

    % SOLVER
        T	= KG\FG;
           
    for i = 1: n_x
        for j = 1 : n_y 
            nid         =   (i-1)+(j-1)*n_x + 1;
            T2d(i,j)    =   T(nid);
            T_ana(i,j)  =   MMS_Solution(GCOORDS(1,nid),GCOORDS(2,nid),...
                lx,ly);
        end
    end

    % PLOTTING
    figure(1)
    surf(XC,YC,T2d);
    title('Steady state soln.');
    xlabel('X')
    ylabel('Y')
    zlabel('T')  
    
    figure(2)
    surf(XC,YC,T_ana);
    title('Analytical soln.');
    xlabel('X')
    ylabel('Y')
    zlabel('T')     
    
    fig = figure(3);
    surf(XC,YC,T2d);
    co = colorbar;
    ylabel(co,'$T[K]$','Interpreter','latex','FontSize',20)
    xlabel('$X[m]$')
    ylabel('$Y[m]$')
    zlabel('$T[K]$')  
    set(findall(gcf,'type','axes'),'FontSize',18,'LineWidth',2)
    set(findall(gcf,'type','line'),'LineWidth',2)
    set(findall(gcf,'type','text'),'Interpreter','latex','FontSize',20)
    set(findall(gcf,'tag','legend'),'Interpreter','latex','FontSize',18,'LineWidth',1)
    set(findall(gcf,'tag','annotation'),'LineStyle','none','Interpreter','latex','FontSize',20)    
    print(gcf,'-dpng','-r300', '2Ddiffeqn_ovs');    

%% COMPUTE L2 Error
    L2_error = 0.0;
    
    for iel = 1:el_tot
        n_now   = EL_N(:,iel);
        dV      = dx*dy;
        
        for i_int = 1:n_int
            T_int  = Nloc(i_int,:)*T(n_now);
            x_int  = Nloc(i_int,:)*GCOORDS(1,n_now)';
            y_int  = Nloc(i_int,:)*GCOORDS(2,n_now)';
            T_mms  = MMS_Solution(x_int,y_int,lx,ly);
            error  = T_int - T_mms;
            L2_error= L2_error + error*error*dV;
        end
    end
    
    L2_error = sqrt( L2_error );
    fprintf(1,'h =          E_2 =  \n' );
    fprintf(1,'%1.4e,  %1.6e ; \n', dx, L2_error );    
        
%% FUNCTIONS
    % compute bilinear shape function wrt local CS
    function val = Nbil2D(xi,eta)                                           
        val1 = (1 - xi).*(1 - eta)/4;
        val2 = (1 - xi).*(1 + eta)/4;
        val3 = (1 + xi).*(1 + eta)/4;
        val4 = (1 + xi).*(1 - eta)/4;
        val  = [val1 ; val2 ; val3 ; val4];
    end
    
    % compute bilinear shape function derivative wrt local CS
    function val = dNbil2D(xi_p,eta_p)                                      
        syms xi eta
        % shape function expressions in symbolic form
        f1          =   (1 - xi).*(1 - eta)/4;
        f2          =   (1 - xi).*(1 + eta)/4;
        f3          =   (1 + xi).*(1 + eta)/4;
        f4          =   (1 + xi).*(1 - eta)/4;
        
        % shape function derivatives in symbolic form
        df1_xi      =   diff(f1,xi);  
        df1_eta     =   diff(f1,eta);
        df2_xi      =   diff(f2,xi);
        df2_eta     =   diff(f2,eta);
        df3_xi      =   diff(f3,xi);
        df3_eta     =   diff(f3,eta);
        df4_xi      =   diff(f4,xi);
        df4_eta     =   diff(f4,eta);

        % compute shape function derivatives at GL quadrature points
        val1_xi     =   vpa(subs(df1_xi,eta,eta_p));
        val1_eta    =   vpa(subs(df1_eta,xi,xi_p));
        val2_xi     =   vpa(subs(df2_xi,eta,eta_p));
        val2_eta    =   vpa(subs(df2_eta,xi,xi_p));
        val3_xi     =   vpa(subs(df3_xi,eta,eta_p));
        val3_eta    =   vpa(subs(df3_eta,xi,xi_p));
        val4_xi     =   vpa(subs(df4_xi,eta,eta_p));
        val4_eta    =   vpa(subs(df4_eta,xi,xi_p));   
        
        % assemble the shape function derivative matrix
        val         =   [val1_xi  val2_xi  val3_xi  val4_xi ; ...
                         val1_eta val2_eta val3_eta val4_eta ];                     
    end    
    
    % compute biquadratic shape function wrt local CS
    function val = Nbiquad2D(xi,eta)                                           
        val1 = xi.*eta.*(xi - 1).*(eta - 1)/4;
        val2 = xi.*(xi - 1).*(1 - eta).*(1 + eta)/2;
        val3 = xi.*eta.*(xi - 1).*(eta + 1)/4;
        val4 = (1 - xi).*(1 + xi).*eta.*(eta + 1)/2;
        val5 = (1 - xi).*(1 + xi).*(1 - eta).*(1 + eta);
        val6 = (1 - xi).*(1 + xi).*eta.*(eta - 1)/2;
        val7 = xi.*(1 + xi).*eta.*(eta - 1)/4;
        val8 = xi.*(1 + xi).*(1 - eta).*(1 + eta)/2;
        val9 = xi.*eta.*(1 + xi).*(1 + eta)/4;
        val  = [val1 ; val2 ; val3 ; val4 ; val5 ; val6 ; val7 ; val8 ; ...
            val9];
    end
    
    % compute biquadratic shape function derivative wrt local CS
    function val = dNbiquad2D(xi_p,eta_p)                                      
        syms xi eta
        % shape function expressions in symbolic form
        f1 = xi.*eta.*(xi - 1).*(eta - 1)/4;
        f2 = xi.*(xi - 1).*(1 - eta).*(1 + eta)/2;
        f3 = xi.*eta.*(xi - 1).*(eta + 1)/4;
        f4 = (1 - xi).*(1 + xi).*eta.*(eta + 1)/2;
        f5 = (1 - xi).*(1 + xi).*(1 - eta).*(1 + eta);
        f6 = (1 - xi).*(1 + xi).*eta.*(eta - 1)/2;
        f7 = xi.*(1 + xi).*eta.*(eta - 1)/4;
        f8 = xi.*(1 + xi).*(1 - eta).*(1 + eta)/2;
        f9 = xi.*eta.*(1 + xi).*(1 + eta)/4;
        
        % shape function derivatives in symbolic form
        df1_xi      =   diff(f1,xi);  
        df1_eta     =   diff(f1,eta);
        df2_xi      =   diff(f2,xi);
        df2_eta     =   diff(f2,eta);
        df3_xi      =   diff(f3,xi);
        df3_eta     =   diff(f3,eta);
        df4_xi      =   diff(f4,xi);
        df4_eta     =   diff(f4,eta);
        df5_xi      =   diff(f5,xi);  
        df5_eta     =   diff(f5,eta);
        df6_xi      =   diff(f6,xi);
        df6_eta     =   diff(f6,eta);
        df7_xi      =   diff(f7,xi);
        df7_eta     =   diff(f7,eta);
        df8_xi      =   diff(f8,xi);
        df8_eta     =   diff(f8,eta);
        df9_xi      =   diff(f9,xi);
        df9_eta     =   diff(f9,eta);         

        % compute shape function derivatives at GL quadrature points
        val1_xi     =   vpa(subs(df1_xi,{xi,eta},{xi_p,eta_p}));
        val1_eta    =   vpa(subs(df1_eta,{xi,eta},{xi_p,eta_p}));
        val2_xi     =   vpa(subs(df2_xi,{xi,eta},{xi_p,eta_p}));
        val2_eta    =   vpa(subs(df2_eta,{xi,eta},{xi_p,eta_p}));
        val3_xi     =   vpa(subs(df3_xi,{xi,eta},{xi_p,eta_p}));
        val3_eta    =   vpa(subs(df3_eta,{xi,eta},{xi_p,eta_p}));
        val4_xi     =   vpa(subs(df4_xi,{xi,eta},{xi_p,eta_p}));
        val4_eta    =   vpa(subs(df4_eta,{xi,eta},{xi_p,eta_p}));
        val5_xi     =   vpa(subs(df5_xi,{xi,eta},{xi_p,eta_p}));
        val5_eta    =   vpa(subs(df5_eta,{xi,eta},{xi_p,eta_p}));
        val6_xi     =   vpa(subs(df6_xi,{xi,eta},{xi_p,eta_p}));
        val6_eta    =   vpa(subs(df6_eta,{xi,eta},{xi_p,eta_p}));
        val7_xi     =   vpa(subs(df7_xi,{xi,eta},{xi_p,eta_p}));
        val7_eta    =   vpa(subs(df7_eta,{xi,eta},{xi_p,eta_p}));
        val8_xi     =   vpa(subs(df8_xi,{xi,eta},{xi_p,eta_p}));
        val8_eta    =   vpa(subs(df8_eta,{xi,eta},{xi_p,eta_p}));         
        val9_xi     =   vpa(subs(df9_xi,{xi,eta},{xi_p,eta_p}));
        val9_eta    =   vpa(subs(df9_eta,{xi,eta},{xi_p,eta_p}));          
        
        % assemble the shape function derivative matrix
        val         =   [val1_xi  val2_xi  val3_xi  val4_xi val5_xi  ...
            val6_xi val7_xi val8_xi val9_xi; ...
                         val1_eta val2_eta val3_eta val4_eta val5_eta ...
            val6_eta val7_eta val8_eta val9_eta];                     
    end        
    
    % compute the jacobian and its determinant
    function [jacobdN,invjacobdN,det_jacobdN] = jacobDN_compute(dNquad,xy_p)
        jacobdN     =   dNquad*xy_p;
        invjacobdN  =   inv(jacobdN);
        det_jacobdN =   det(jacobdN);
    end    
    
    % Manufactured soln.
    function val = MMS_Solution(x,y,lx,ly)
        val = sin(pi*x/lx).*cos(pi*y/ly);
    end
    
    % RHS of manufactured soln. 
    function val = MMS_RHS(x,y,lx,ly)
        val = sin(pi*x/lx).*cos(pi*y/ly).*(pi*pi/(lx*lx) + ...
            (pi*pi)/(ly*ly));
    end
    
    
    
    
    
    
        