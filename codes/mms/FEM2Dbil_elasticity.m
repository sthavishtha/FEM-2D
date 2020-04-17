% FEM code for order of verification test in 2D elasticity equation
% Modified from the template code as a part of the course "Introduction to
% Finite element methods in geosciences", Spring 2019, ETH Zurich.
% (c) Sthavishtha, 2019
% -------------------------------------------------------------------------
% Features
% - implemented bilinear shape functions
% - 2x2 Gauss-legendre quadrature numerical integration schemes
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
    Dcoeff      =   E/((1 + poisson)*(1 - 2*poisson));
    D           =   D*Dcoeff;                                               % constitutive formulation
    
% NUMERICAL PARAMETERS
    el_x        =   16;                                                     % #elements in x direction
    el_y        =   16;
    el_tot      =	el_x*el_y;                                              % #elements total
    n_x         =	el_x + 1;                                               % total #nodes along x direction
    n_y         =	el_y + 1;                                               % total #nodes along y direction
    n_per_el    =   4;                                                      % #nodes per element (bilinear/biquadratic)
    n_tot       =	n_x*n_y;                                                % #nodes total
    EQN         =   2;                                                      % #EQN DOF per node (ux,uy in this case)
    
% QUADRATURE PARAMETERS
    % 2 point GLg
    n_int       =   4;                                                      % #integration points  
    int_pts   	=   [-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3) ; ...        % integration points in [(-1,-1) ; (1,1)]
                    -1/sqrt(3)  1/sqrt(3) 1/sqrt(3) -1/sqrt(3)];                        
    weight      =   [1  1  1  1];                                           % weights of integration points
                                                                         
% CREATE NUMERICAL GRID
    dx          =   Lx/el_x;                                                % distance between two nodes
    dy          =   Ly/el_y;           
    GCOORDS     =   zeros(2,n_x*n_y);                                       % global node coordinates
    GLINDEX     =   zeros(n_x,n_y);                                         % global node index   
    XC          =   zeros(n_x,n_y);                                         % x,y coordinates in 2D form
    YC          =   zeros(n_x,n_y);    
    XC_ana      =   zeros(n_x,n_y);                                         % x,y coordinates in 2D form
    YC_ana      =   zeros(n_x,n_y);                                         % analytical mesh coordinates

% COMPUTING GLOBAL COORDINATES 
    for i = 1: n_x
        for j = 1 : n_y
            xn = (i-1)*dx;
            yn = (j-1)*dy;
            nid = (i-1)+(j-1)*n_x + 1;
            GCOORDS(1,nid) = xn;
            GCOORDS(2,nid) = yn;          
            GLINDEX(i,j) = nid;
        end
    end
        
% LOCAL-TO-GLOBAL MAPPING  
    EL_N = zeros(n_per_el,el_tot);
    el = 0;                                                                 % relates local to global node numbers per element       
    for j = 1 : n_y - 1
        for i = 1 : n_x - 1        
            el = el + 1;
            EL_N(1,el)   =   GLINDEX(i,j);
            EL_N(2,el)   =   GLINDEX(i,j + 1);
            EL_N(3,el)   =   GLINDEX(i + 1,j + 1);
            EL_N(4,el)   =   GLINDEX(i + 1,j);             
        end
    end

% NODES-TO-EQUATION NO. MAPPING
    NF = zeros(EQN,n_tot);
    node_eqn = 0;                                                           % relates between nodes and equation numbers
    for j = 1 : n_y 
        for i = 1 : n_x             
            node_eqn            =   node_eqn + 1;
            NF(1,GLINDEX(i,j))  =   node_eqn;                               
            node_eqn            =   node_eqn + 1;
            NF(2,GLINDEX(i,j))  =   node_eqn;
        end
    end
    
% ELEMENTS-TO-EQUATION NO. MAPPING    
    GG = zeros(n_per_el*EQN,el_tot);
    for iel = 1 : el_tot                                                    % relates between elements and equation numbers
    node_eqn    =   0;        
        for dof_n = 1 : EQN
            node_eqn = node_eqn + 1;
            GG(node_eqn,iel) = NF(dof_n,EL_N(1,iel));            
        end  
        for dof_n = 1 : EQN
            node_eqn = node_eqn + 1;
            GG(node_eqn,iel) = NF(dof_n,EL_N(2,iel));
        end
        for dof_n = 1 : EQN
            node_eqn = node_eqn + 1;        
            GG(node_eqn,iel) = NF(dof_n,EL_N(3,iel));
        end
        for dof_n = 1 : EQN
            node_eqn = node_eqn + 1;        
            GG(node_eqn,iel) = NF(dof_n,EL_N(4,iel));
        end
    end
	
% BOUNDARY CONDITIONS  
    i_b = 0;
    i_t = 0;
    i_l = 0;
    i_r = 0;
    for i = 1 : n_x
        for j = 1 : n_y
            if GCOORDS(2,GLINDEX(i,j)) == 0 % bottom wall
                i_b = i_b + 1;
                bc_dof_b(i_b)       = NF(1,GLINDEX(i,j));
                bc_dof_b(i_b + 1)   = NF(2,GLINDEX(i,j));  
                bc_val_b(i_b)       = MMS_Solutionx(GCOORDS(1,GLINDEX(i,j)),GCOORDS(2,GLINDEX(i,j)));
                bc_val_b(i_b + 1)   = MMS_Solutiony(GCOORDS(1,GLINDEX(i,j)),GCOORDS(2,GLINDEX(i,j)));
                i_b = i_b + 1;
            elseif GCOORDS(2,GLINDEX(i,j)) == Ly % top wall
                i_t = i_t + 1;
                bc_dof_t(i_t)       = NF(1,GLINDEX(i,j));
                bc_dof_t(i_t + 1)   = NF(2,GLINDEX(i,j));
                bc_val_t(i_t)       = MMS_Solutionx(GCOORDS(1,GLINDEX(i,j)),GCOORDS(2,GLINDEX(i,j))); 
                bc_val_t(i_t + 1)   = MMS_Solutiony(GCOORDS(1,GLINDEX(i,j)),GCOORDS(2,GLINDEX(i,j))); 
                i_t = i_t + 1;
            elseif GCOORDS(1,GLINDEX(i,j)) == 0 % left wall
                i_l = i_l + 1;
                bc_dof_l(i_l)       = NF(1,GLINDEX(i,j));
                bc_dof_l(i_l + 1)   = NF(2,GLINDEX(i,j));
                bc_val_l(i_l)       = MMS_Solutionx(GCOORDS(1,GLINDEX(i,j)),GCOORDS(2,GLINDEX(i,j)));
                bc_val_l(i_l + 1)   = MMS_Solutiony(GCOORDS(1,GLINDEX(i,j)),GCOORDS(2,GLINDEX(i,j)));
                i_l = i_l + 1;
            elseif GCOORDS(1,GLINDEX(i,j)) == Lx % right wall
                i_r = i_r + 1;
                bc_dof_r(i_r)       = NF(1,GLINDEX(i,j));
                bc_dof_r(i_r + 1)   = NF(2,GLINDEX(i,j));
                bc_val_r(i_r)       = MMS_Solutionx(GCOORDS(1,GLINDEX(i,j)),GCOORDS(2,GLINDEX(i,j)));
                bc_val_r(i_r + 1)   = MMS_Solutiony(GCOORDS(1,GLINDEX(i,j)),GCOORDS(2,GLINDEX(i,j)));   
                i_r = i_r + 1;
            end
        end
    end
    
    bc_dof      =     [bc_dof_b bc_dof_t bc_dof_l bc_dof_r];
    bc_val      =     [bc_val_b bc_val_t bc_val_l bc_val_r];     
          
% INITIALIZATION OF ALL KINDS OF STUFF    
    KG          =	zeros(n_tot*EQN,n_tot*EQN);                             % global stiffness matrix   
    FG          =	zeros(n_tot*EQN,1);                                     % global force vector      
    Nloc        =   zeros(n_int,n_per_el);                                  % local shape function matrix
    dNloc       =   zeros(2,n_per_el,n_int);                                % local shape function derivative matrix 
    for i_int = 1 : n_int
        Nloc(i_int,:)      =   Nbil2D(int_pts(1,i_int),int_pts(2,i_int));
        dNloc(:,:,i_int)   =   dNbil2D(int_pts(1,i_int),int_pts(2,i_int));        
    end       				 
     
    for i = 1: n_x
        for j = 1 : n_y
            nid     =   (i-1)+(j-1)*n_x + 1;
            XC(i,j) =   GCOORDS(1,nid);
            YC(i,j) =   GCOORDS(2,nid);
            XC_ana(i,j) =   GCOORDS(1,nid);
            YC_ana(i,j) =   GCOORDS(2,nid);
        end
    end    
    
  figure(1)
  ZC = zeros(n_x,n_y);
  pcolor(XC,YC,ZC);
  title('Initial');
  xlabel('X')
  ylabel('Y')     

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
        x_int   =   Nloc(i_int,:)*GCOORDS(1,n_now)';                        % integration points in global coords.
        y_int   =   Nloc(i_int,:)*GCOORDS(2,n_now)'; 
        source_x=   MMS_RHSx(poisson,x_int,y_int,Dcoeff);                   % sources from MMS_RHS 
        source_y=   MMS_RHSy(poisson,x_int,y_int,Dcoeff);
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
  for i = 1:length(bc_dof)
    KG(bc_dof(i),:)        =   0;                    % set bc-rows to zero
    KG(bc_dof(i),bc_dof(i))=   1;                    % put 1 onto the main diagonal
    FG(bc_dof(i))          =   bc_val(i);            % set bc-value to the rhs
  end

% SOLVER
  U	= KG\FG;
  
  U_x = U(1:2:n_tot*EQN - 1);
  U_y = U(2:2:n_tot*EQN);
  
% COMPUTING DEFORMED GRID      
  for i = 1: n_x
      for j = 1 : n_y 
          XC(i,j)       =   XC(i,j) + U(NF(1,GLINDEX(i,j)));
          YC(i,j)       =   YC(i,j) + U(NF(2,GLINDEX(i,j)));
          XC_ana(i,j)   =   XC_ana(i,j) + MMS_Solutionx(GCOORDS(1,GLINDEX(i,j)),...
              GCOORDS(2,GLINDEX(i,j)));
          YC_ana(i,j)   =   YC_ana(i,j) + MMS_Solutiony(GCOORDS(1,GLINDEX(i,j)),...
              GCOORDS(2,GLINDEX(i,j)));          
      end
  end
  
  %% PLOTTING
  figure(2)
  ZC = zeros(n_x,n_y);
  pcolor(XC,YC,ZC);
  title('FEM soln.');
  xlabel('X')
  ylabel('Y')  
  
  figure(3)
  ZC_ana = zeros(n_x,n_y);
  pcolor(XC_ana,YC_ana,ZC_ana);
  title('MMS soln.');
  xlabel('X')
  ylabel('Y')  
  
  figure(4);
  pcolor(XC,YC,ZC);
  xlabel('$X[m]$')
  ylabel('$Y[m]$')
  set(findall(gcf,'type','axes'),'FontSize',18,'LineWidth',2)
  set(findall(gcf,'type','line'),'LineWidth',2)
  set(findall(gcf,'type','text'),'Interpreter','latex','FontSize',20)
  set(findall(gcf,'tag','legend'),'Interpreter','latex','FontSize',18,'LineWidth',1)
  set(findall(gcf,'tag','annotation'),'LineStyle','none','Interpreter','latex','FontSize',20)
  print(gcf,'-dpng','-r300', '2Delastic_ovs');
  
% COMPUTE THE L2 ERROR
  L2_error = 0.0;
  
  for iel = 1:el_tot
      n_now   = EL_N(:,iel);
      dV      = dx*dy;
            
      for i_int = 1:n_int           
        ux_int  = Nloc(i_int,:)*U_x(n_now);
        uy_int  = Nloc(i_int,:)*U_y(n_now);
        x_int   = Nloc(i_int,:)*GCOORDS(1,n_now)';
        y_int   = Nloc(i_int,:)*GCOORDS(2,n_now)';         
        ux_mms  = MMS_Solutionx(x_int,y_int);
        uy_mms  = MMS_Solutiony(x_int,y_int);
        error_x = ux_int - ux_mms;
        error_y = uy_int - uy_mms;
        L2_error = L2_error + (error_x*error_x + error_y*error_y)*dV;        
      end      
  end
  
  L2_error   = sqrt(L2_error);  
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
    
    % compute the jacobian and its determinant
    function [jacobdN,invjacobdN,det_jacobdN] = jacobDN_compute(dNquad,xy_p)
        jacobdN     =   dNquad*xy_p;
        invjacobdN  =   inv(jacobdN);
        det_jacobdN =   det(jacobdN);
    end    
    
    % Manufactured Soln.
    function val = MMS_Solutionx(x,y)
        val = sin(x).*cos(y);       
    end
    
    function val = MMS_Solutiony(x,y)
        val = cos(x).*sin(y);
    end
    
    % RHS of manufactured Soln. 
    function val = MMS_RHSx(poisson,x,y,Dcoeff)
        val = -(2*poisson - 2)*sin(x).*cos(y)*Dcoeff;
    end
    
    function val = MMS_RHSy(poisson,x,y,Dcoeff)
        val = -(2*poisson - 2)*sin(y).*cos(x)*Dcoeff;
    end    
    
    
    
    
    
    
        