% FEM code for Assessment
% Modified from the code by Marcel Frehner, ETH Zurich, 2017
% (c) Sthavishtha, 2019
% -------------------------------------------------------------------------
% Features 
% - Transient solution of a 2D diffusion equation
% - constant grid spacing, source terms and diffusivity
% - bilinear shape functions with 4 nodes/element
% - 2nd & 3rd order Gauss-legendre quadrature numerical integration (just 
%   comment/uncomment the quadrature scheme you want to use)
% -------------------------------------------------------------------------

%% GENERAL STUFF
    clear                                                                   % clear the current Workspace
    close all                                                               % close all figure windows
    clc                                                                     % clear the Command Window
    
% GEOMETRICAL PARAMETERS
    Lx          =	01;                                                     % length of model [m]
    Ly          =   01;
    Time_tot    =   1;                                                      % computation time    
    
% PHYSICAL PARAMETERS
    kappa_x     =	0.01;                                                   % thermal diffusivity components [m2/s]
    kappa_y     =	0.01;                                  
    kappa       =	[kappa_x 0 ; 0 kappa_y];                                % thermal diffusivity matrix [m2/s]
    source      =	0;                                                      % heat source [K/s]
    
% NUMERICAL PARAMETERS
    el_x        =   29;                                                     % #elements in x direction
    el_y        =   29;
    ndof        =   1;
    el_tot      =	el_x*el_y;                                              % #elements total
    n_x         =	ndof*el_x + 1;                                          % total #nodes along x direction
    n_y         =	ndof*el_y + 1;                                          % total #nodes along y direction
    n_per_el    =   4;                                                      % #nodes per element (bilinear/biquadratic)
    n_tot       =	n_x*n_y;                                                % #nodes total
    
% QUADRATURE PARAMETERS
    % 2 point GLg
%     n_int       =   4;                                                      % #integration points  
%     int_pts   	=   [-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3) ; ...        % integration points in [(-1,-1) ; (1,1)]
%                     -1/sqrt(3)  1/sqrt(3) 1/sqrt(3) -1/sqrt(3)];                        
%     weight      =   [1  1  1  1];                                           % weights of integration points  

    % 3 point GLg
    n_int      =   9;                                                       % #integration points
    int_pts    =   [-0.774596669241483  -0.774596669241483 -0.774596669241483 ...
        0.0 0.0 0.0 ...
        0.774596669241483 0.774596669241483 0.774596669241483 ;
                    -0.774596669241483  0.0 0.774596669241483 ...
        -0.774596669241483  0.0 0.774596669241483 ...
        -0.774596669241483 0.0 0.774596669241483];                          % integration points in [-1,1]
    weight     =   [0.3086 0.4938 0.3086 0.4938 0.7901 0.4938 0.3086 0.4938 ...
        0.3086];                                                            % weights of integration points
                                                                         
% CREATE NUMERICAL GRID
    dx          =   Lx/(n_x - 1);                                           % distance between two nodes
    dy          =   Ly/(n_y - 1);
    dt          =   0.2;                                                    % incremental time step            
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
            GLINDEX(j,i) = nid;
        end
    end
        
% LOCAL-TO-GLOBAL MAPPING  
    EL_N = zeros(n_per_el,el_tot);
    el = 0;                                                                 % relates local to global node numbers per element       
    for i = 1 : ndof : n_y - ndof
        for j = 1 : ndof : n_x - ndof        
            el = el + 1;
            EL_N(1,el)   =   GLINDEX(i,j);
            EL_N(2,el)   =   GLINDEX(i + 1,j);  
            EL_N(3,el)   =   GLINDEX(i + 1,j + 1);
            EL_N(4,el)   =   GLINDEX(i,j + 1);            
        end
    end
	
% BOUNDARY CONDITIONS
    bc_dof_b    =   1:n_x;                                                  % bottom boundaries
    bc_dof_t    =   n_tot - n_x + 1 : n_tot;                                % top boundaries
    bc_dof_l    =   n_x + 1 : n_x : n_tot - 2*n_x + 1;                      % left boundaries
    bc_dof_r    =   2*n_x : n_x : n_tot - n_x;                              % right boundaries
    bc_dof      =   horzcat(bc_dof_b,bc_dof_t,bc_dof_l,bc_dof_r);           % dof's to which Dirichlet bc's are assigned                            
    bc_val      =   zeros(size(bc_dof));                                    % value for dirichlet dof's
    
% INITIALIZATION OF ALL KINDS OF STUFF       
    T2d         =   zeros(n_x,n_y);                                         % temperature in 2D mesh grid format
    LG          =	zeros(n_tot,n_tot);                                     % global stiffness matrix
    RG          =	zeros(n_tot,n_tot);                                     % global RHS vector (T^n terms in transient state)    
    FG          =	zeros(n_tot,1);                                         % global force vector      
    Nloc        =   zeros(n_int,n_per_el);                                  % local shape function matrix
    dNloc       =   zeros(2,n_per_el,n_int);                                % local shape function derivative matrix
    time        =   0;                                                      % time looping variable     
    for i_int = 1 : n_int
        Nloc(i_int,:)      =   Nbil2D(int_pts(1,i_int),int_pts(2,i_int));
        dNloc(:,:,i_int)   =   dNbil2D(int_pts(1,i_int),int_pts(2,i_int));        
    end       				 
     
     for i = 1: n_x
        for j = 1 : n_y
            nid     =   (i-1)+(j-1)*n_x + 1;
            XC(i,j) =   GCOORDS(1,nid);
            YC(i,j) =   GCOORDS(2,nid);
        end
     end    
     
     % initial condition for temperature
     T_new       =   100.0*exp(-(XC - 0.5*Lx).^2 - (YC - 0.5*Ly).^2);
     T           =   reshape(T_new,n_x*n_y,1);

%% LOOPING
while time < Time_tot
    
    for iel = 1:el_tot % ELEMENT LOOP
        n_now   =   EL_N(:,iel);                                            % which nodes are in the current element?
        
        Mloc    =   zeros(n_per_el,n_per_el);
        Kloc    =   zeros(n_per_el,n_per_el);
        Floc    =   zeros(n_per_el,1);

        for i_int = 1:n_int        
            [jacobdN,invjacobdN,det_jacobdN] = jacobDN_compute(dNloc(:,:,i_int), ...
                GCOORDS(:,n_now(:))');
            dNlocT  =   invjacobdN*dNloc(:,:,i_int);                        % derivative of shape fn. wrt global coords
            Mloc    =   Mloc + Nloc(i_int,:)'*Nloc(i_int,:)*weight(i_int)*det_jacobdN;
            Kloc    =   Kloc + dNlocT'*kappa*dNlocT*weight(i_int)*det_jacobdN;  % local stiffness matrix
            Floc    =   Floc + source*Nloc(i_int,:)'*weight(i_int)*det_jacobdN; % local force vector          
            Lloc    =   Mloc/dt + Kloc;
            Rloc    =   Mloc/dt;                                            % right hand side vector (transient state)
        end

        % add local matrices to global ones
        LG(n_now,n_now)     =   LG(n_now,n_now) + Lloc;
        RG(n_now,n_now)     =   RG(n_now,n_now) + Rloc;    
        FG(n_now)           =   FG(n_now) + Floc;
    end % ELEMENT LOOP
        b   =   RG*T + FG;

    % APPLY BOUNDARY CONDITIONS
    for i = 1:length(bc_dof)
        LG(bc_dof(i),:)        =   0;                                       % set bc-rows to zero
        LG(bc_dof(i),bc_dof(i))=   1;                                       % put 1 onto the main diagonal
        b(bc_dof(i))           =   bc_val(i);                               % set bc-value to the rhs
    end

    % SOLVER
        T	= LG\b; 

    % increment time
        time = time + dt;          
    
    time
    
    for i = 1: n_x
        for j = 1 : n_y 
            nid = (i-1)+(j-1)*n_x + 1;
            T2d(i,j) = T(nid);
        end
    end

    % PLOTTING        
%     surf(XC,YC,T2d);
%     xlabel('X')
%     ylabel('Y')
%     zlabel('T')  
%     drawnow
end          

% SAVING PLOTS
    figure(1)
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
    print(gcf,'-dpng','-r300', '2Ddiffeqn_surf');    
    
    figure(2)
    [c,h] = contour(XC,YC,T2d);
    clabel(c,h,'FontSize',10,'Color','red')
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
    print(gcf,'-dpng','-r300', '2Ddiffeqn_cont');  
        
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
    
    
    
    
    
    
        