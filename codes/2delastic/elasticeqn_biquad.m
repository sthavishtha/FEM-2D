% FEM code for Assessment
% Modified from the code by Marcel Frehner, ETH Zurich, 2017
% (c) Sthavishtha, 2019
% -------------------------------------------------------------------------
% Features : 
% - Solution of 2D elasticity equation in 1 time step
% - no time-dependece terms used, as its steady state
% - constant grid spacing and source terms
% - biquadratic shape functions with 9 nodes 
% - 2nd & 3rd order numerical integration to compute the element matrices
% -------------------------------------------------------------------------

%% GENERAL STUFF
    clear                                                                   % clear the current Workspace
    close all                                                               % close all figure windows
    clc                                                                     % clear the Command Window
    
% GEOMETRICAL PARAMETERS
    Lx          =	01;                                                     % length of model [m]
    Ly          =   01;   
    
% PHYSICAL PARAMETERS
    poisson     =	0.49;                                                   % poisson ratio
    E           =	10^5;                                                   % young's modulus of elasticity
    D           =	[1 - poisson    poisson     0 ; ...
                     poisson        1 - poisson 0 ; ...
                     0              0           (1 - 2*poisson)/2];         
    D           =   D*E/((1 + poisson)*(1 - 2*poisson));                    % constitutive formulation
    source_x    =	0;                                                      % source 
    source_y    =	0;
    
% NUMERICAL PARAMETERS
    el_x        =   32;                                                     % #elements in x direction
    el_y        =   32;
    el_tot      =	el_x*el_y;                                              % #elements total
    ndof_per_n  =   2;
    n_x         =	ndof_per_n*el_x + 1;                                    % total #nodes along x direction
    n_y         =	ndof_per_n*el_y + 1;                                    % total #nodes along y direction
    n_per_el    =   9;                                                      % #nodes per element (bilinear/biquadratic)
    n_tot       =	n_x*n_y;                                                % #nodes total
    EQN         =   2;                                                      % #EQN DOF per node (ux,uy in this case)
    
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
    GCOORDS     =   zeros(2,n_x*n_y);                                       % global node coordinates
    GLINDEX     =   zeros(n_y,n_x);                                         % global node index   
    XC          =   zeros(n_x,n_y);                                         % x,y coordinates in 2D form
    YC          =   zeros(n_x,n_y);    

% COMPUTING GLOBAL COORDINATES 
    for i = 1: n_x
        for j = 1 : n_y
            xn = (i-1)*dx;
            yn = (j-1)*dy;
            nid = (i-1)+(j-1)*n_x + 1;
            GCOORDS(1,nid) = xn;
            GCOORDS(2,nid) = yn;                       
            GLINDEX(j,i) = nid;
        end
    end
        
% LOCAL-TO-GLOBAL MAPPING  
    EL_N = zeros(n_per_el,el_tot);
    el = 0;                                                                 % relates local to global node numbers per element       
    for j = 1 : ndof_per_n : n_y - ndof_per_n
        for i = 1 : ndof_per_n : n_x - ndof_per_n      
            el = el + 1;           
            EL_N(1,el)   =   GLINDEX(j,i);
            EL_N(2,el)   =   GLINDEX(j + 1,i);  
            EL_N(3,el)   =   GLINDEX(j + 2,i);
            EL_N(4,el)   =   GLINDEX(j + 2,i + 1);
            EL_N(5,el)   =   GLINDEX(j + 1,i + 1);  
            EL_N(6,el)   =   GLINDEX(j,i + 1);  
            EL_N(7,el)   =   GLINDEX(j,i + 2);
            EL_N(8,el)   =   GLINDEX(j + 1,i + 2);  
            EL_N(9,el)   =   GLINDEX(j + 2,i + 2);               
        end
    end

% NODES-TO-EQUATION NO. MAPPING
    NF = zeros(EQN,n_tot);
    node_eqn = 0;                                                           % relates between nodes and equation numbers
    for j = 1 : n_y 
        for i = 1 : n_x             
            node_eqn            =   node_eqn + 1;
            NF(1,GLINDEX(j,i))  =   node_eqn;                               
            node_eqn            =   node_eqn + 1;
            NF(2,GLINDEX(j,i))  =   node_eqn;
        end
    end
    
% ELEMENTS-TO-EQUATION NO. MAPPING    
    GG = zeros(n_per_el*EQN,el_tot);
    for iel = 1 : el_tot                                                    % relates between elements and equation numbers
    node_eqn    =   0;     
        for j = 1 : n_per_el
            for dof_n = 1 : EQN
                node_eqn = node_eqn + 1;
                GG(node_eqn,iel) = NF(dof_n,EL_N(j,iel));            
            end  
        end
    end
	
% BOUNDARY CONDITIONS
    % in global CS
    bc_dof_b    =   find(GCOORDS(2,:) == 0);
    bc_dof_t    =   find(GCOORDS(2,:) == Ly);
    bc_dof_l    =   find(GCOORDS(1,:) == 0);
    bc_dof_r    =   find(GCOORDS(1,:) == Lx);
    
    % including the x,x...,y,y.. components respectively in order
    bc_dof      =   [NF(1,bc_dof_b) NF(2,bc_dof_b) NF(1,bc_dof_t) NF(2,bc_dof_t)];
    bc_val      =   [zeros(1,length(bc_dof_b)) zeros(1,length(bc_dof_b)) ...
                      0.5*ones(1,length(bc_dof_t)) zeros(1,length(bc_dof_t))];
          
% INITIALIZATION OF ALL KINDS OF STUFF    
    KG          =	zeros(n_tot*EQN,n_tot*EQN);                             % global stiffness matrix   
    FG          =	zeros(n_tot*EQN,1);                                     % global force vector      
    Nloc        =   zeros(n_int,n_per_el);                                  % local shape function matrix
    dNloc       =   zeros(2,n_per_el,n_int);                                % local shape function derivative matrix 
    for i_int = 1 : n_int
        Nloc(i_int,:)      =   Nbiquad2D(int_pts(1,i_int),int_pts(2,i_int));
        dNloc(:,:,i_int)   =   dNbiquad2D(int_pts(1,i_int),int_pts(2,i_int));        
    end       				 
     
    for i = 1: n_x
        for j = 1 : n_y
            nid     =   (i-1)+(j-1)*n_x + 1;
            XC(i,j) =   GCOORDS(1,nid);
            YC(i,j) =   GCOORDS(2,nid);
        end
    end    
    
  figure(1)
  ZC = zeros(n_x,n_y);
  pcolor(XC,YC,ZC);
  xlabel('$X[m]$')
  ylabel('$Y[m]$')
  set(findall(gcf,'type','axes'),'FontSize',18,'LineWidth',2)
  set(findall(gcf,'type','line'),'LineWidth',2)
  set(findall(gcf,'type','text'),'Interpreter','latex','FontSize',20)
  set(findall(gcf,'tag','legend'),'Interpreter','latex','FontSize',18,'LineWidth',1)
  set(findall(gcf,'tag','annotation'),'LineStyle','none','Interpreter','latex','FontSize',20)
  print(gcf,'-dpng','-r300', '2Delastic_init');   

%% LOOPING
  for iel = 1:el_tot % ELEMENT LOOP
    n_now   =   EL_N(:,iel);                                                % which nodes are in the current element?
    eqn_now =   GG(:,iel);                                                  % equation nos. in the current element
    
    Bmatrix =   zeros(3,n_per_el*EQN);
    Kloc    =   zeros(n_per_el*EQN,n_per_el*EQN);
    Floc    =   zeros(n_per_el*EQN,1);

    for i_int = 1:n_int        
        [jacobdN,invjacobdN,det_jacobdN] = jacobDN_compute(dNloc(:,:,i_int), ...
                GCOORDS(:,n_now(:))');               
        dNlocT  =   invjacobdN*dNloc(:,:,i_int);                            % shape functions in physical coordinates           
        Bmatrix(1,1:2:n_per_el*EQN - 1) = dNlocT(1,:);
        Bmatrix(2,2:2:n_per_el*EQN) = dNlocT(2,:);
        Bmatrix(3,1:2:n_per_el*EQN - 1) = dNlocT(2,:);
        Bmatrix(3,2:2:n_per_el*EQN) = dNlocT(1,:);
        Kloc    =   Kloc + Bmatrix'*D*Bmatrix*weight(i_int)*det_jacobdN;    % local stiffness matrix
        Floc1   =   source_x*Nloc(i_int,:)'*weight(i_int)*det_jacobdN;
        Floc2   =   source_y*Nloc(i_int,:)'*weight(i_int)*det_jacobdN;        
        Floc(1:2:n_per_el*EQN-1,1) =   Floc(1:2:n_per_el*EQN-1,1) + Floc1(:);
        Floc(2:2:n_per_el*EQN,1)   =   Floc(2:2:n_per_el*EQN,1) + Floc2(:); % local force vector                
    end

    % add local matrices to global ones
    KG(eqn_now,eqn_now)     =   KG(eqn_now,eqn_now) + Kloc;
    FG(eqn_now)             =   FG(eqn_now) + Floc;
  end % ELEMENT LOOP

% APPLY BOUNDARY CONDITIONS 
  % ordering used in x,y,x,y
  for i = 1:length(bc_dof)
    KG(bc_dof(i),:)        =   0;                                           % set bc-rows to zero
    KG(bc_dof(i),bc_dof(i))=   1;                                           % put 1 onto the main diagonal
    FG(bc_dof(i))          =   bc_val(i);                                   % set bc-value to the rhs
  end

% SOLVER
  U	= KG\FG;           

%% PLOTTING        
  for i = 1: n_x
      for j = 1 : n_y  
          XC(i,j)       =   XC(i,j) + U(NF(1,GLINDEX(j,i)));
          YC(i,j)       =   YC(i,j) + U(NF(2,GLINDEX(j,i)));
      end
  end
  
  figure(2)
  ZC = zeros(n_x,n_y);
  pcolor(XC,YC,ZC);
  xlabel('$X[m]$')
  ylabel('$Y[m]$')
  set(findall(gcf,'type','axes'),'FontSize',18,'LineWidth',2)
  set(findall(gcf,'type','line'),'LineWidth',2)
  set(findall(gcf,'type','text'),'Interpreter','latex','FontSize',20)
  set(findall(gcf,'tag','legend'),'Interpreter','latex','FontSize',18,'LineWidth',1)
  set(findall(gcf,'tag','annotation'),'LineStyle','none','Interpreter','latex','FontSize',20)
  print(gcf,'-dpng','-r300', '2Delastic_final');             
    
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
    function mat = dNbil2D(xi_p,eta_p)                                      
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
        
        % assemble the shape function derivative & dNexp matrix
        mat         =   [val1_xi  val2_xi  val3_xi  val4_xi ; ...
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
    function mat = dNbiquad2D(xi_p,eta_p)                                      
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
        mat         =   [val1_xi  val2_xi  val3_xi  val4_xi val5_xi  ...
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
    
    
    
    
    
    
        