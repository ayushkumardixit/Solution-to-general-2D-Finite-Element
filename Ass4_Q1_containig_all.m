% Assignment 4 Question 1
% Ayush Kumar Dixit 18117021 FEM
Method=input('Enter 1 for three noded triangle,2 for four Noded rectangle,3 for four noded quadrilateral,4 for eight noded quadrilateral, 5 for six noded triangular elements:');
% Start switch for case 1
switch Method
    case 1
        % Three Noded Triangluar Element
% Define Number of Elements
NoElements=input('Enter the number of elements:');
% Define Number of Nodes
NoNodes=input('Enter Total number of nodes:');
% Define Corresponding Nodes
CorresNodes=zeros(NoElements,6);
% Take input of corresponding nodes from element
CorresNodes=input('Enter the Nodes corresponding to global nodes:');
% Take input the value of youngs modulus
E=input('Enter the value of Youngs Modulus:');
% Take input the value of poissons ratio
Nue=input('Enter the value of Poissons ratio:');
% Take input the value of force
force=zeros(2*NoNodes,1);
force=input('Enter the force in form of matrix:');
% Assuming Plane Stress
% Define the D matrix in terms of youngs modulus and poisoons ratio
D=(E/(1-Nue^2))*[1 Nue 0; Nue 1 0; 0 0 (1-Nue^2)/2];
% Define a matrix containing the coordinates
X_Y=zeros(3*NoElements,2);
% Input a matrix containing the coordinates
X_Y=input('Enter the Nodal matrix coordinates 3NoElements*2:');
% Define alpha values in matrix form
Alpha=zeros(3,NoElements);
AlphaC=[0 1 -1; 1 -1 0];
for e=1:NoElements;
for i=1:3;
    Alpha(i,e)=Alpha(i,e)+ X_Y(3*e + AlphaC(1,i)-1,1)*X_Y(3*e+ AlphaC(2,i)-1,2)- X_Y(3*e + AlphaC(2,i)-1,1)*X_Y(3*e + AlphaC(1,i)-1,2);
end 
end
% Define beta values in matrix form
Beta =zeros(3, NoElements);
for e=1:NoElements;
for i=1:3;
  Beta(i,e)= Beta(i,e) + X_Y(3*e + AlphaC(1,i)-1,2)- X_Y(3*e+ AlphaC(2,i)-1,2);
end 
end 
% Define gamma values in matrix form
Gamma =zeros(3, NoElements);
for e=1:NoElements;
for i=1:3;
  Gamma(i,e)= Gamma(i,e) + X_Y(3*e + AlphaC(2,i)-1,1) - X_Y(3*e + AlphaC(1,i)-1,1) ;
end 
end 
% Define the value of determinent as doublearea
DoubleArea=zeros(NoElements,1);
% Calculate the determinent and store the value in doublearea matrix
for e=1:NoElements;
M=[1 X_Y(2*e -1,1) X_Y(2*e-1,2);1 X_Y(2*e,1) X_Y(2*e,2);1 X_Y(2*e+1,1) X_Y(2*e+1,2)]
DoubleArea(e,1)= DoubleArea(e,1) + det(M);
end
% Define global matrix
KGlobal=zeros(2*NoNodes,2*NoNodes);
% Define a matrix for calculating the strain
Bforstrain=zeros(3*NoElements,6);
% Create a loop and calculate the value of local stiffness matrix
for e=1:NoElements;
    % Define B matrix
    B=(1/DoubleArea(e,1))*[Beta(1,e) 0 Beta(2,e) 0 Beta(3,e) 0; 0 Gamma(1,e) 0 Gamma(2,e) 0 Gamma(3,e);Gamma(1,e) Beta(1,e) Gamma(2,e) Beta(2,e) Gamma(3,e) Beta(3,e)];
    % Define local stiffness matrix
    KLocal= B'*D*B*(DoubleArea(e,1)/2);
    Bforstrain(3*e-2:3*e,:)= Bforstrain(3*e-2:3*e,:)+ B;
    %Add and calculate the value of global stiffness matrix
    KGlobal(CorresNodes(e,:),CorresNodes(e,:))=KGlobal(CorresNodes(e,:),CorresNodes(e,:))+KLocal;
end
% Define displacements matrix
displacements=zeros(2*NoNodes,1);
% Input a matrix containing fixed node places
fixedDof=input('Input the nodes that are fixed in row matrix:')
% Define Active Degree of freedom matrix 
activeDof=setdiff([1:2*NoNodes]',fixedDof');
% Calculte the displacements matrix
U(activeDof,1)=KGlobal(activeDof,activeDof)\force(activeDof,1);
% Store the values obtained in already defined displacements matrix
displacements(activeDof,1)= U(activeDof,1);
% Define strain of elements
Strain=zeros(3*NoElements,1);
% Define stress of elements
Stress=zeros(3*NoElements,1);
% Create a loop for calculating strain and stress
for e=1:NoElements;
    % Define local strain for an element
StrainLocal= (Bforstrain(3*e-2:3*e,:))*displacements(CorresNodes(e,:),1);
% Add Strain to already stored matrix
Strain(3*e-2:3*e,1)=Strain(3*e-2:3*e,1)+StrainLocal;
% Define local stress for an element
StressLocal=D*StrainLocal;
% Add stress to already stored matrix value
Stress(3*e-2:3*e,1)=Stress(3*e-2:3*e,1)+StressLocal;
end
% Display the values of displacements at specific nodes
fprintf('Displacements %6.12f\n',double(displacements));
% Display the values of strain at specific nodes
fprintf('Strain %6.12f\n',double(Strain));
% Display the values of stress at specific nodes
fprintf('Stress %6.12f\n',double(Stress));    
% End the loop
end
% Start switch for case 2
switch Method
    case 2
        % Four Noded Rectangular Element
% Enter the number of elements
NoElements=input('Enter the number of elements:');
% Enter the number of nodes
NoNodes=input('Enter Total number of nodes:');
% Define a corresponding nodes matrix
CorresNodes=zeros(NoElements,6);
% Take input values of corrsponding nodes matrix
CorresNodes=input('Enter the Nodes corresponding to global nodes:');
% Input the values of youngs modulus
E=input('Enter the value of Youngs Modulus:');
% Input the values of poisonns ratio
Nue=input('Enter the value of Poissons ratio:');
% Define the force vectors
force=zeros(2*NoNodes,1);
% Input the force vectors
force=input('Enter the force in form of matrix:');
% Assuming Plane Stress
% Define the D in matrix form
D=(E/(1-Nue^2))*[1 Nue 0; Nue 1 0; 0 0 (1-Nue^2)/2];
% Define the matrix containg the coordinates
X_Y=zeros(3*NoElements,2);
% Input the coordinates of all the nodes
X_Y=input('Enter the Nodal matrix: ');
% Define a global stiffness matrix
KGlobal=zeros(2*NoNodes,2*NoNodes);
% Using loop calculate and store values in stiffness matrix
for e=1:NoElements;
    % Define a symbolic variable x
    syms x;
    % Define a symbolic variable y
    syms y;
    % Define the value of length l
    l=X_Y(4*e-2,1)-X_Y(4*e-3,1);
    % Define the value of height h
    h=X_Y(4*e-1,2)-X_Y(4*e-2,2);
    % Define the B matrix
    B=[(1/l)*(y/h -1) 0 (1/l)*(-y/h + 1) 0 y/(l*h) 0 (-y/(l*h)) 0; 0 (1/h)*(x/l -1) 0 (-x/(l*h)) 0 x/(l*h) 0 (1/h)*(-x/l +1);(1/h)*(x/l -1) (1/l)*(y/h -1) (-x/(l*h)) (1/l)*(-y/h + 1) x/(l*h) y/(l*h) (1/h)*(-x/l +1) (-y/(l*h))];
    % Define the variable a1
    a1=X_Y(4*e-3,1); b1=X_Y(4*e-2,1);
    % Define the variable a2
    a2=X_Y(4*e-2,2); b2=X_Y(4*e-1,2);
    % Define the KL matrix for further integration
    KL = (B')*D*B;
    % Define the x and y for further use in integration
    x=(a1+b1)/2 + (1/sqrt(3))*(b1-a1)/2; y=(a2+b2)/2 + (1/sqrt(3))*(b2-a2)/2;
    % Calculate the value at first position
    KL1=subs(KL);
    x=(a1+b1)/2 + (1/sqrt(3))*(b1-a1)/2; y=(a2+b2)/2 - (1/sqrt(3))*(b2-a2)/2;
    % Calculate the value at first position
    KL2=subs(KL);
    % Define the x and y for further use in integration
    x=(a1+b1)/2 - (1/sqrt(3))*(b1-a1)/2; y=(a2+b2)/2 + (1/sqrt(3))*(b2-a2)/2;
    % Calculate the value at first position
    KL3=subs(KL);
    % Define the x and y for further use in integration
    x=(a1+b1)/2 - (1/sqrt(3))*(b1-a1)/2; y=(a2+b2)/2 - (1/sqrt(3))*(b2-a2)/2;
    % Calculate the value at first position
    KL4=subs(KL);
    % Add the value to calculate the local matrix using integration
    KLocal= 0.25*(b2-a2)*(b1-a1)*(KL1 + KL2 + KL3 + KL4);
    % Clear the values of x and y defined
    clear x; clear y;
    % Add and store respective values in global matrix
    KGlobal(CorresNodes(e,:),CorresNodes(e,:))=KGlobal(CorresNodes(e,:),CorresNodes(e,:))+ KLocal;
end
% Define displacements matrix
displacements=zeros(2*NoNodes,1);
% Deifne the matrix containg fixed nodes
fixedDof=input('Input the nodes that are fixed in row matrix:')
% Define the active degree of freedom 
activeDof=setdiff([1:2*NoNodes]',fixedDof');
% Calculate the displacement matrix
U(activeDof,1)=KGlobal(activeDof,activeDof)\force(activeDof,1);
% Store the calculted value in  already defined displacements matrix
displacements(activeDof,1)= U(activeDof,1);
% Define a strain matrix for an elements
Strain=zeros(3,1);
% Define a stress matrix for an elements
Stress=zeros(3,1);
% Input the elemnt on which you want the stress and strain
e=input('The Element on which you want the strain and stress:');
   % Define a symbolic variable x
    syms x;
    % Define a symbolic variable y
    syms y;
    % Define the value of length l
    l=X_Y(4*e-2,1)-X_Y(4*e-3,1);
    % Define the value of height h
    h=X_Y(4*e-1,2)-X_Y(4*e-2,2);
    % Define the B matrix
    B=[(1/l)*(y/h -1) 0 (1/l)*(-y/h + 1) 0 y/(l*h) 0 (-y/(l*h)) 0; 0 (1/h)*(x/l -1) 0 (-x/(l*h)) 0 x/(l*h) 0 (1/h)*(-x/l +1);(1/h)*(x/l -1) (1/l)*(y/h -1) (-x/(l*h)) (1/l)*(-y/h + 1) x/(l*h) y/(l*h) (1/h)*(-x/l +1) (-y/(l*h))];
    % Calculate the strain matrix
    Strain=B*displacements(CorresNodes(e,:),1);  
    % Caculate the stress matrix
Stress=D*Strain;
% Display the values of displacements
fprintf('Displacements %6.12f\n',double(displacements));
% Display the values of strain
disp(Strain);
% Display the values of stress
disp(Stress);
end

switch Method
    case 3
        % Four Noded Quadrilateral Element
% Take input the number of elements
NoElements=input('Enter the number of elements:');
% Take input the number of nodes involved
NoNodes=input('Enter Total number of nodes:');
% Define a matrix of corresponding nodes
CorresNodes=zeros(NoElements,8);
% Take input the corresponding nodes from user
CorresNodes=input('Enter the Nodes corresponding to global nodes in matrix form NoELements*8:');
% Take input the value of youngs modulus
E=input('Enter the value of Youngs Modulus:');
% Take input the value of poissons ratio
Nue=input('Enter the value of Poissons ratio:');
% Define a matrix containg forces at each nodes
force=zeros(2*NoNodes,1);
% Take input the force applied in a matrix
force=input('Enter the force in form of matrix:');
% Assuming Plane Stress
D=(E/(1-Nue^2))*[1 Nue 0; Nue 1 0; 0 0 (1-Nue^2)/2];
% Define the cooordinates of all the nodes
X_Y=zeros(8*NoElements,2);
% Take input the coordinates of all the nodes
X_Y=input('Enter the Nodal matrix 4NoElements*2:');
% Define the B Matrx
B1=[1 0 0 0;0 0 0 1;0 1 1 0];
% Define the global stiffness matrix
KGlobal=zeros(2*NoNodes);
    % Define the Local jacobian matrix element 11
JLocal11=0;
    % Define the Local jacobian matrix element 12
JLocal12=0;
% Define the Local jacobian matrix element 21
JLocal21=0;
% Define the Local jacobian matrix element 22
JLocal22=0;
% Define a loop for storing and calculating global stiffness matrix
 for e=1:NoElements;
     % Define  a symbolic variable eta 
     syms eta;
     % Define a symbolic variable shi
     syms shi;
     % Define a matrix containing shape function
N=[0.25*(1-eta)*(1-shi); 0.25*(1-eta)*(1+shi); 0.25*(1+eta)*(1+shi); 0.25*(1+eta)*(1-shi)];
% Define a loop for creating and storing jacobian matrix
   for p=1:4;
    JLocal11 = JLocal11 + diff(N(p,1),shi)*X_Y(4*e+p-4,1);
    JLocal12 = JLocal12 + diff(N(p,1),shi)*X_Y(4*e+p-4,2);
    JLocal21 = JLocal21 + diff(N(p,1),eta)*X_Y(4*e+p-4,1);
    JLocal22 = JLocal22 + diff(N(p,1),eta)*X_Y(4*e+p-4,2);
   end
   % Storing values again in previosly defined jacobian matrix
   JLocal(1,1)=JLocal11;
   JLocal(1,2)=JLocal12;
   JLocal(2,1)=JLocal21;
   JLocal(2,2)=JLocal22;
   % Defining a B2 matrix 
    B2=(1/det(JLocal))*[JLocal22 (-JLocal12) 0 0;(-JLocal21) JLocal11 0 0; 0 0 JLocal22 (-JLocal12); 0 0 (-JLocal21) JLocal11];
    % Defining a B3 matrix 
    B3=0.25*[eta-1 0 1-eta 0 1+eta 0 (-1-eta) 0; shi-1 0 (-1-shi) 0 1+shi 0 1-shi 0; 0 eta-1 0 1-eta 0 1+eta 0 (-1-eta);0 shi-1 0 (-1-shi) 0 1+shi 0 1-shi];
    % Defining and calculating B matrix as the multiplicaton of B1,B2 & B3
    B=B1*B2*B3;
    % Define a KL matrix
    KL = B'*D*B*det(JLocal);
    % Guass Quadrature Integration
    % Define the first point of integration
    eta= -0.57735; shi= -0.57735;
    % Take the value obatined in KL1 matrix
    KL1=subs(KL);
    % Define the second point of integration
    eta=0.57735; shi= -0.57735;
    % Take the value obatined in KL2 matrix
    KL2=subs(KL);
    % Define the third point of integration
    eta= -0.57735; shi= 0.57735; 
    % Take the value obatined in KL3 matrix
    KL3=subs(KL);
    % Define the fourth point of integration
    eta= 0.57735;shi=0.57735;
    % Take the value obatined in KL4 matrix
    KL4=subs(KL);
    % Calculate the integral and insert value in KLocal matrix
    KLocal= KL1 + KL2 + KL3 + KL4;
    % Calculate and store values in gobal stiffness matrix
    KGlobal(CorresNodes(e,:),CorresNodes(e,:))= KGlobal(CorresNodes(e,:),CorresNodes(e,:))+ KLocal;
    % Clear the symbolic variables
    clear eta; clear shi;
 end
 % Define a displacements matrix
displacements=zeros(2*NoNodes,1);
% Take the fixed number of nodes form the user
fixedDof=input('Input the nodes that are fixed in row matrix:')
% Define the active degree of freedom
activeDof=setdiff([1:2*NoNodes]',fixedDof');
% Calculate the displacement matrix
U(activeDof,1)=KGlobal(activeDof,activeDof)\force(activeDof,1);
% Store the calculated values in already defined displacements matrix
displacements(activeDof,1)= U(activeDof,1);
% Display the values of displacements
fprintf('Displacements %6.12f\n',double(displacements));
% Take input the element on which you want to calculate the strain & stress
e=input('The Element on which you want the strain and stress:');
% Define a symbolic variable eta
     syms eta;
     % Define a symbolic variable shi
     syms shi;
     % Define a shape function matrix in form of eta and shi
N=[0.25*(1-eta)*(1-shi); 0.25*(1-eta)*(1+shi); 0.25*(1+eta)*(1+shi); 0.25*(1+eta)*(1-shi)];
% Create a loop to calulate the jacobian matrix
   for p=1:4;
    JLocal11 = JLocal11 + diff(N(p,1),shi)*X_Y(4*e+p-4,1);
    JLocal12 = JLocal12 + diff(N(p,1),shi)*X_Y(4*e+p-4,2);
    JLocal21 = JLocal21 + diff(N(p,1),eta)*X_Y(4*e+p-4,1);
    JLocal22 = JLocal22 + diff(N(p,1),eta)*X_Y(4*e+p-4,2);
   end
   % Store the calculted value in local jacobian matrix 11
   JLocal(1,1)=JLocal11;
   % Store the calculted value in local jacobian matrix 12
   JLocal(1,2)=JLocal12;
   % Store the calculted value in local jacobian matrix 21
   JLocal(2,1)=JLocal21;
   % Store the calculted value in local jacobian matrix 22
   JLocal(2,2)=JLocal22;
   % Define a matrix B2
    B2=(1/det(JLocal))*[JLocal22 (-JLocal12) 0 0;(-JLocal21) JLocal11 0 0; 0 0 JLocal22 (-JLocal12); 0 0 (-JLocal21) JLocal11];
    % Define a matrix B3 in terms of eta and shi
    B3=0.25*[eta-1 0 1-eta 0 1+eta 0 (-1-eta) 0; shi-1 0 (-1-shi) 0 1+shi 0 1-shi 0; 0 eta-1 0 1-eta 0 1+eta 0 (-1-eta);0 shi-1 0 (-1-shi) 0 1+shi 0 1-shi];
    % Define a matrix B after multiplying B1,B2 and B3
    B=B1*B2*B3;
    % Define a strain matrix
Strain=B*displacements(CorresNodes(e,:),1); 
% Define a stress matrix
Stress=D*Strain;
% Display the displacments values
fprintf('Displacements %6.12f\n',double(displacements));
% Display the strain on the element
disp(Strain);
%  Display the stress on the element
disp(Stress);
% End the switch loop
end

switch Method
    case 4
        % Eight Noded Quadrilateral Element
% Enter the number of elements
NoElements=1;
% Enter the number of nodes
NoNodes=input('Enter Total number of nodes:');
% Take input the value of youngs modulus
E=input('Enter the value of Youngs Modulus:');
% Take input the value of poisonns ratio
Nue=input('Enter the value of Poissons ratio:');
% Define a force matrix
force=zeros(2*NoNodes,1);
% Take input the values of forces applied at each node
force=input('Enter the force in form of matrix:');
% Assuming Plane Stress
% Define the D matrix
D=(E/(1-Nue^2))*[1 Nue 0; Nue 1 0; 0 0 (1-Nue^2)/2];
% Define the Input Coordinates matrix
X_Y=zeros(8*NoElements,2);
% Take input the coordinates
X_Y=input('Enter the Nodal matrix 8NoElements*2:');
% Define B1
B1=[1 0 0 0;0 0 0 1;0 1 1 0];
% Define Global Stiffness Matrix
KGlobal=zeros(2*NoNodes);
% Define Jacobian Elements at 11
JLocal11=0;
% Define Jacobian Elements at 12
JLocal12=0;
% Define Jacobian Elements at 21
JLocal21=0;
% Define Jacobian Elements at 22
JLocal22=0;
% Define a symbolic variable eta
     syms eta;
     % Define a symbolic variable shi
     syms shi;
     % Define a matrix of shape function
  N=[0.25*(1-shi)*(1-eta)*(-1-shi-eta); 0.25*(1+shi)*(1-eta)*(-1+shi-eta); 0.25*(1+shi)*(1+eta)*(-1+shi+eta); 0.25*(1-shi)*(1+eta)*(-1-shi+eta); 0.5*(1-shi*shi)*(1-eta); 0.5*(1+shi)*(1-eta*eta); 0.5*(1-shi*shi)*(1+eta); 0.5*(1-shi)*(1-eta*eta)];
 % Define a for loop to calculate the local jacobian element
  for p=1:8 
    JLocal11 = JLocal11 + diff(N(p,1),shi)*X_Y(p,1);
    JLocal12 = JLocal12 + diff(N(p,1),shi)*X_Y(p,2);
    JLocal21 = JLocal21 + diff(N(p,1),eta)*X_Y(p,1);
    JLocal22 = JLocal22 + diff(N(p,1),eta)*X_Y(p,2);
  end
  % Define the determinent of jacobian matrix
    det=(JLocal11*JLocal22-JLocal12*JLocal21);
    % Define the B2 matrix
    B2=(1/det)*[JLocal22 (-JLocal12) 0 0;(-JLocal21) JLocal11 0 0; 0 0 JLocal22 (-JLocal12); 0 0 (-JLocal21) JLocal11];
    % Define the B3 matrix in the form of eta and shi
    B3=[diff(N(1,1),shi) 0 diff(N(2,1),shi) 0 diff(N(3,1),shi) 0 diff(N(4,1),shi) 0 diff(N(5,1),shi) 0 diff(N(6,1),shi) 0 diff(N(7,1),shi) 0 diff(N(8,1),shi) 0; diff(N(1,1),eta) 0 diff(N(2,1),eta) 0 diff(N(3,1),eta) 0 diff(N(4,1),eta) 0 diff(N(5,1),eta) 0 diff(N(6,1),eta) 0 diff(N(7,1),eta) 0 diff(N(8,1),eta) 0; 0 diff(N(1,1),shi) 0 diff(N(2,1),shi) 0 diff(N(3,1),shi) 0 diff(N(4,1),shi) 0 diff(N(5,1),shi) 0 diff(N(6,1),shi) 0 diff(N(7,1),shi) 0 diff(N(8,1),shi) ;0 diff(N(1,1),eta) 0 diff(N(2,1),eta) 0 diff(N(3,1),eta) 0 diff(N(4,1),eta) 0 diff(N(5,1),eta) 0 diff(N(6,1),eta) 0 diff(N(7,1),eta) 0 diff(N(8,1),eta)];
    % Calculate the B matrix
    B= B1*B2*B3;
    % Define the KL matrix for further use in integration
    KL= B'*D*B*det;
    % Guass Quadrature Integration
    % Define the coordinates at first point
    eta= -(1/sqrt(3)); shi= -(1/sqrt(3));
    % Calculate the value at first point
    KL1=subs(KL);
    % Define the coordinates at first point
    eta=(1/sqrt(3));   shi= -(1/sqrt(3));
    % Calculate the value at first point
    KL2=subs(KL);
    % Define the coordinates at first point
    eta=-(1/sqrt(3));  shi= (1/sqrt(3)); 
    % Calculate the value at first point
    KL3=subs(KL);
    % Define the coordinates at first point
    eta= (1/sqrt(3));  shi=(1/sqrt(3));
    % Calculate the value at first point
    KL4=subs(KL);
    % Add the values obtained at all four points 
    KLocal= KL1 + KL2 + KL3 + KL4;
    % clear the symbolic variables
    clear eta; clear shi;
    % Define global stiffness matrix
    KGlobal(:,:)= KLocal;
    % Define the displacements matrix
displacements=zeros(2*NoNodes,1);
% Define the active degree of freedom matrix
activeDof=[3 4 5 6  8 9 10 11 12 13 14 16];
% Calculate the displacements matrix
U(activeDof,1)=KGlobal(activeDof,activeDof)\force(activeDof,1);
% Store the calculted values in already defined displacements matrix
displacements(activeDof,1)= U(activeDof,1);
% Display the values of displacements at each nodes
fprintf('Displacements %6.12f\n',double(displacements));
% End the loop
end

switch Method
    case 5
        % Six Noded Triangluar Element
% Take input the number of elements
NoElements=input('Enter the number of elements:');
% Take input the number of nodes
NoNodes=input('Enter Total number of nodes:');
% Define a matrix of corresponding nodes
CorresNodes=zeros(NoElements,12);
% Take input the corresponding nodes in matrix form
CorresNodes=input('Enter the Nodes corresponding to global nodes in matrix form NoELements*8:');
% Take input the value of youngs modulus
E=input('Enter the value of Youngs Modulus:');
% Take input the value of poisoons ratio
Nue=input('Enter the value of Poissons ratio:');
% Assuming Plane Stress
% Define the D matrix
D=(E/(1-Nue^2))*[1 Nue 0; Nue 1 0; 0 0 (1-Nue^2)/2];
% Define a matrix containing coordinates
X_Y=zeros(6*NoElements,2);
% Take input the respective coordinates
X_Y=input('Enter the Nodal matrix 4NoElements*2:');
% Take input the value of B1 matrix
B1=[1 0 0 0;0 0 0 1;0 1 1 0];
% Define global stiffness matrix
KGlobal=zeros(2*NoNodes);
% Define the elements of jacobian matrix
JLocal11=0;
JLocal12=0;
JLocal21=0;
JLocal22=0;
% Define a force matrix applied at each nodes 
force=zeros(2*NoNodes,1);
% Define a symbolic varibles xf and yf
syms xf; syms yf;
% Take input the area from user
A=input('Enter the area of triangle where force is applied:');
% Define the values L1f, L2f and L3f
L1f= (1/(2*A))*((X_Y(8,1)*X_Y(9,2) - X_Y(9,1)*X_Y(8,2) + xf*(X_Y(8,2)- X_Y(9,2))+ yf*(X_Y(9,1)-X_Y(8,1))));
L2f= (1/(2*A))*((X_Y(8,1)*X_Y(7,2) - X_Y(7,1)*X_Y(9,2) + xf*(X_Y(9,2)- X_Y(7,2))+ yf*(X_Y(7,1)-X_Y(9,1))));
L3f= (1/(2*A))*((X_Y(7,1)*X_Y(8,2) - X_Y(8,1)*X_Y(7,2) + xf*(X_Y(7,2)- X_Y(8,2))+ yf*(X_Y(8,1)-X_Y(7,1))));
% Define the shape function matrix
Nf=  [L1f*(2*L1f-1); L2f*(2*L2f-1); (1-L1f-L2f)*(1-2*L1f-2*L2f); 4*L1f*L2f; 4*L2f*(1-L1f-L2f); 4*L1f*(1-L1f-L2f)];
% Use the coordinate value xf to be 4
xf=4;
% Define and calculate the value N11
N11=subs(Nf(1,1));
% Define N11 to be a functon of varible y
fun=@(y)N11;
% Integrate and store value in variable F1
F1= int(fun,yf,0,3);
% Define and calculate the value N12
N12=subs(Nf(2,1));
% Define N11 to be a functon of varible y
fun2=@(y)N12;
% Integrate and store value in variable F2
F2=int(fun2,yf,0,3);
% Define and calculate the value N14
N14=subs(Nf(4,1));
% Define N11 to be a functon of varible y
fun4=@(y)N14;
% Integrate and store value in variable F4
F4=int(fun4,yf,0,3);
% Define the force values using above calculated values
force=[-1.5;0;0;0;F1;0;F4;0;F2;0;0;0;-1.5;0;0;0;0;0];
% Clear all the defined symbolic variables
clear xf; clear L1f; clear L2f; clear L3f; clear Nf;
% Define a for loop to calculate the global stiffness matrix
 for e=1:NoElements;
     % Define a symbolic variable L1
     syms L1;
     % Define a symbolic variable L2
     syms L2;
     % Define a shape function in matrix form
N=  [L1*(2*L1-1); L2*(2*L2-1); (1-L1-L2)*(1-2*L1-2*L2); 4*L1*L2; 4*L2*(1-L1-L2); 4*L1*(1-L1-L2)];
% Create a loop to calculate the different elements of jacobian matrix
  for p=1:6 
    JLocal11 = JLocal11 + diff(N(p,1),L1)*X_Y(6*e-6+p,1);
    JLocal12 = JLocal12 + diff(N(p,1),L1)*X_Y(6*e-6+p,2);
    JLocal21 = JLocal21 + diff(N(p,1),L2)*X_Y(6*e-6+p,1);
    JLocal22 = JLocal22 + diff(N(p,1),L2)*X_Y(6*e-6+p,2);
  end
  % Calculte the determinent of jacobian matrix
    det=(JLocal11*JLocal22-JLocal12*JLocal21);
    % Define the B2 matrix
    B2=(1/det)*[JLocal22 (-JLocal12) 0 0;(-JLocal21) JLocal11 0 0; 0 0 JLocal22 (-JLocal12); 0 0 (-JLocal21) JLocal11];
   % Define the B3 matrix
    B3=[diff(N(1,1),L1) 0 diff(N(2,1),L1) 0 diff(N(3,1),L1) 0 diff(N(4,1),L1) 0 diff(N(5,1),L1) 0 diff(N(6,1),L1) 0 ; diff(N(1,1),L2) 0 diff(N(2,1),L2) 0 diff(N(3,1),L2) 0 diff(N(4,1),L2) 0 diff(N(5,1),L2) 0 diff(N(6,1),L2) 0 ; 0 diff(N(1,1),L1) 0 diff(N(2,1),L1) 0 diff(N(3,1),L1) 0 diff(N(4,1),L1) 0 diff(N(5,1),L1) 0 diff(N(6,1),L1) ;0 diff(N(1,1),L2) 0 diff(N(2,1),L2) 0 diff(N(3,1),L2) 0 diff(N(4,1),L2) 0 diff(N(5,1),L2) 0 diff(N(6,1),L2)];
   % Calculate and define the B matrix
    B=B1*B2*B3;
    % Define a varoable Kl for further use in interation
    KL = B'*D*B*det;
    % Guass Quadrature Integration
    % Define the first point for integration
    L1= 1/3; L2= 1/3;
    % Calculate the value at first point
    KL1=(-27/48)*subs(KL);
    % Define the second point for integration
    L1=0.6; L2= 0.2;
    % Calculate the value at second point
    KL2=(25/48)*subs(KL);
    % Define the third point for integration
    L1= 0.2; L2= 0.6; 
    % Calculate the value at third point
    KL3=(25/48)*subs(KL);
    % Define the fourth point for integration
    L1= 0.2; L2= 0.2; 
    % Calculate the value at fourth point
    KL4=(25/48)*subs(KL);
    % Calculate and store the obtained value in Local matrix
    KLocal= KL1+KL2+KL3+KL4;
    % Display the local stiffness matrix
    disp(KLocal);
    % Calculate and add values in global stiffness matrix
    KGlobal(CorresNodes(e,:),CorresNodes(e,:))= KGlobal(CorresNodes(e,:),CorresNodes(e,:))+ KLocal;
    % Clear the symboic variables L1 and L2
    clear L1; clear L2;
 end
 % Define the displacements matrix
displacements=zeros(2*NoNodes,1);
% Define the fixed degree of freedom matrix
fixedDof=input('Input the nodes that are fixed in row matrix:')
% Define the matrix containg active degree of freedom values
activeDof=setdiff([1:2*NoNodes]',fixedDof');
% Calculate the displacements matrix
U(activeDof,1)=KGlobal(activeDof,activeDof)\force(activeDof,1);
% Store the calculted values i already defined dsiaplcements matrix
displacements(activeDof,1)= U(activeDof,1);
% Display the displacements matrix
fprintf('Displacements %6.12f\n',double(displacements));
% End the loop
end