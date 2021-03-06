%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab Code supplementing the paper
% A biomechanical approach to infer size-based functional response in aquatic and terrestrial systems
% by Portalier, Cherif, Fussmann, Loreau 
%
% Frontiers in Ecology and Evolution
%
% August 2021
%
% Matlab version: R2020b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%%%% READ ME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The following code defines the function fnSpecies
% that calculates species traits
% The function takes body mass and a vector of parameters

%%%% Physical parameters taken as arguments:
% 1. Body density (kg.m-3)
% 2. Medium density (kg.m-3)
% 3. Medium dynamic viscosity (N.s.m-2)
% 4. Acceleration due to gravity (m.s-2)
%
%%%% Parameters defined within the function
% see main text or supplementary methods for references
% 1. Ash-free dry mass to wet mass ratio
% 2. Energy (J.kg-1) to ash-free dry mass ratio
% 3. Bite diameter (m) at reference size
% 4. Reference mass (kg) for bite size
% 5. Detection distance (m) at reference size
% 6. Reference mass (kg) for detection distance
%
%%%% The function returns a vector with
% 1. Body volume (m3)
% 2. Body mass (kg)
% 3. Body radius (m)
% 4. Section surface of the body (m2)
% 5. Stroke period (s)
% 6. Energetic content (J)
% 7. Bite size (kg)
% 8. Bite time (s)
% 9. Detection distance (m)
% 10. Species-specific speed (m.s-1)
% 11. Motion cost per time at species-specific speed (J.s-1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Species] = fnSpecies(BodyMass,param)
  %% unwrap physical parameters
  BodyDensity=param(1);
  MediumDensity=param(2);
  DynamicViscosity=param(3);
  Gravity=param(4);
  
  %% Define other parameters
  % Ash-free dry mass to wet mass ratio
  DryMassRatio=0.16;
  % Energy (J.kg-1) to ash-free dry mass ratio
  EnergyDryMassRatio=23000000; 
  % Bite diameter (m) at reference size
  ReferenceBiteSize=0.00026; 
  % Reference mass (kg) for bite size
  BiteSizeReferenceMass=2.9; 
  % Detection distance (m) at reference size
  ReferenceDetectionDistance=0.255; 
  % Reference mass (kg) for detection distance
  ReferenceSizeDetectiondistance=0.0376;
  
  %% response vector
  Species=zeros(1,11);
  
  % Body volume
  Species(1)=BodyMass/BodyDensity;
  % Body mass
  Species(2)=BodyMass;
  % Body radius
  Species(3)=(Species(1)*3.0/4.0/pi)^(1.0/3.0);
  % Section surface
  Species(4)= (Species(3)^2.0)*pi;
  % Stroke period
  Species(5)=0.35*(BodyMass^0.25);
  % Energy content
  Species(6)=BodyMass*DryMassRatio*EnergyDryMassRatio;
  % Bite size
  BiteRadius=(ReferenceBiteSize*((BodyMass/BiteSizeReferenceMass)^0.32))/2.0;
  Species(7)=4.0/3.0*(BiteRadius^3)*pi*BodyDensity;
  % Bite time
  Species(8)= 0.1*(Species(7))^2;
  % Detection distance
  Species(9)=ReferenceDetectionDistance*((Species(1)/ReferenceSizeDetectiondistance)^0.49);
  
  % Species-specific speed calculation
  % Time step for computation
  TimeStep=Species(5)/1000;
  % Switch parameter for fnMotion function
  % i_Switch = 1: fnMotion is called for the calculation of species-specific speed
  i_Switch=1;
  % Vector of parameters for fnMotion function (see corresponding file for details)
  p=[Species(1),BodyMass,Species(3),Species(4),Species(5),Gravity,MediumDensity,DynamicViscosity,TimeStep,i_Switch];
  Result=fnMotion(BodyMass,p);
  % Species-specific speed
  Species(10)=Result(1);
  % Motion cost per time at species-specific speed
  Species(11)=Result(2);
end
