% Takes an array of values and steps once in diffusion using Dirichlet's
function[return_arr] = dirich_diff(patterned_cell, index, square_dim, vis)
    u = cellfun(@(v) v(index), patterned_cell(:,:));
    
    nx = square_dim;            % Number of 
    ny = square_dim;
    nt = 1;            % 42 hrs of culture and 30 minute intervals
    dt = .1;           
    % Representative of a 30 min interval
    dx = square_dim/(nx-1);    
    dy = square_dim/(ny-1);
    x=-square_dim/2:dx:square_dim/2;                       %Range of x(0,2) and specifying the grid points
    y=-square_dim/2:dy:square_dim/2;                       %Range of y(0,2) and specifying the grid points
    u=zeros(nx,ny);                  %Preallocating u
    un=zeros(nx,ny);                 %Preallocating un
    u = cellfun(@(v) v(index), patterned_cell(:,:));
    UW=0;                            %x=0 Dirichlet B.C 
    UE=0;                            %x=L Dirichlet B.C 
    US=0;                            %y=0 Dirichlet B.C 
    UN=0;                            %y=L Dirichlet B.C     
    % Differential Equations
    %B.C vector
    bc=zeros(nx-2,ny-2);
    bc(1,:)=UW/dx^2; 
    bc(nx-2,:)=UE/dx^2;  %Dirichlet B.Cs
    bc(:,1)=US/dy^2; 
    bc(:,ny-2)=UN/dy^2;  %Dirichlet B.Cs

    bc(1,1)=UW/dx^2+US/dy^2; 
    bc(nx-2,1)=UE/dx^2+US/dy^2;
    bc(1,ny-2)=UW/dx^2+UN/dy^2; 
    bc(nx-2,ny-2)=UE/dx^2+UN/dy^2;
    bc=vis*dt*bc;

    %Calculating the coefficient matrix for the implicit scheme
    Ex=sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
    Ax=Ex+Ex'-2*speye(nx-2);        %Dirichlet B.Cs
    Ey=sparse(2:ny-2,1:ny-3,1,ny-2,ny-2);
    Ay=Ey+Ey'-2*speye(ny-2);        %Dirichlet B.Cs
    A=kron(Ay/dy^2,speye(nx-2))+kron(speye(ny-2),Ax/dx^2);
    D=speye((nx-2)*(ny-2))-vis*dt*A;

    un=u;
    
    U=un;U(1,:)=[];U(end,:)=[];U(:,1)=[];U(:,end)=[];
    U=reshape(U+bc,[],1);
    U=D\U;
    U=reshape(U,nx-2,ny-2);
    u(2:nx-1,2:ny-1)=U;
    %Boundary conditions
    %Dirichlet:
    u(1,:)=UW;
    u(nx,:)=UE;
    u(:,1)=US;
    u(:,ny)=UN;
    
    
    return_arr = u;
    %lowestValue = min(U(U(:)>0));
    %highestValue = max(U(:));
    %imagesc(U);
    %cmap = jet(256);
    %caxis(gca, [lowestValue-2/256, highestValue]);
    %cmap(1,:) = [0,0,0];
    %colormap(cmap);
    %colorbar;
    

end