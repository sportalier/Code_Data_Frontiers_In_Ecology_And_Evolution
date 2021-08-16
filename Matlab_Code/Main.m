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
% The following code is used to run the model defined in the paper
% It calls several functions defined in other files
% It defines body masses for n predator species and n prey species
% Then, the code computes for each species: 
% 1. A vector of traits (see fnSpecies function for details)
% The n vectors for the predators are stored into a cell array (PredatorTraits)
% The n vectors for the prey are stored into a cell array (PreyTraits)
% 2. Two arrays of values for predator and prey capture sequence (see fnMotion file for details)
% The n arrays for the predators are stored into a cell array (PredatorArraytot)
% The n arrays for the prey are stored into a cell array (PreyArraytot)
%
% Several vectors of size n are computed to define all aspects of predator-prey relationship
% 1. parameter beta for encounter rate (required for attack rate)
% 2. capture probability (required for attack rate)
% 3. capture time
% 4. Handling time
% 5. attack rate
% 6. total time
%
% Details on calculation and full references for values are provided in 
% the main text or in the supplementary materials of the paper
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Species body mass
% The example is done over 2 predator/prey pairs only
% But any number of species should work
% WARNING: we assume that predator i feeds on prey i 
% Thus, the two vectors have to be sorted accordingly
PredatorMass=[1e-4,1e-2];
PreyMass=[1e-8,1e-5];
NumberofPredators=2; % number of prey = number of predators (see above)

%%%%%%%% Physical parameters 
% Physical parameters (sea water at 20 C) (change values for air)
BodyDensity=1080; % (kg.m-3)
MediumDensity= 1026.95; % 1.247;% (kg.m-3)
DynamicViscosity= 0.00139; % (N.s.m-2) 
Gravity=9.8; % (m.s-2)

param=[BodyDensity,MediumDensity,DynamicViscosity,Gravity];

%%%%%%%% Species traits
% see fnSpecies for details about the different values of the returned vector
PredatorTraits=cell(1,NumberofPredators);

for i=1:NumberofPredators
    PredatorTraits{i}=fnSpecies(PredatorMass(i),param);
end

PreyTraits=cell(1,NumberofPredators);
for i=1:NumberofPredators
    PreyTraits{i}=fnSpecies(PreyMass(i),param);
end

%%%%%%%% Array for capture sequence
PredatorArraytot=cell(1,NumberofPredators);

for i=1:NumberofPredators
    % Retrieve species traits
    Sp=PredatorTraits{i};
    % Time step for computation
    TimeStep=Sp(5)/100;
    % Switch variable for fnMotion function
    % i_Switch = 2: fnMotion is called for the computation of the capture sequence
    i_Switch=2;
    % Vector of parameters for fnMotion function (see fnMotion file for details)
    p=[Sp(1),PredatorMass(i),Sp(3),Sp(4),Sp(5),Gravity,MediumDensity,DynamicViscosity,TimeStep,i_Switch];
    PredatorArraytot{i}=fnMotion(PredatorMass(i),p);
end

PreyArraytot=cell(1,NumberofPredators);

for i=1:NumberofPredators
    % Retrieve species traits
    Sp=PreyTraits{i};
    % Time step for computation
    TimeStep=Sp(5)/100;
    % Switch variable for fnMotion function
    % i_Switch = 2: fnMotion is called for the computation of the capture sequence
    i_Switch=2;
    % Vector of parameters for fnMotion function (see fnMotion file for details)
    p=[Sp(1),PreyMass(i),Sp(3),Sp(4),Sp(5),Gravity,MediumDensity,DynamicViscosity,TimeStep,i_Switch];
    PreyArraytot{i}=fnMotion(PreyMass(i),p);
end

%%%%%%%% Searching sequence
% Defines a vector for encounter rate (parameter "beta", see main text and supplementary material)  
% that is used to compute attack rate
Betaparam=zeros(NumberofPredators,1);

for i=1:NumberofPredators
    Predator=PredatorTraits{i};
    Prey=PreyTraits{i};
    %%% encounter rate (beta parameter)
    % If the predator moves faster than the prey
    if Predator(10)>Prey(10)
        Betaparam(i,1)=(pi*Predator(9)^2.0)*((Prey(10)^2.0)+3.0*(Predator(10)^2.0))/(3.0*Predator(10));
    else
        % if the prey moves faster than the predator
        Betaparam(i,1)=(pi*Predator(9)^2.0)*((Predator(10)^2.0)+3.0*(Prey(10)^2.0))/(3.0*Prey(10));
    end  
end

%%%%%%%% Capture sequence
% Define a vector for capture probability (required to compute attack rate)
% and a vector for capture time
% number of attempts = 1 / capture probability
captureprobability=zeros(NumberofPredators,1);
capturetime=zeros(NumberofPredators,1);

for i=1:NumberofPredators
    PredatorArray=PredatorArraytot{i};
    Predator=PredatorTraits{i};
    
    PreyArray=PreyArraytot{i};
    Prey=PreyTraits{i};
    % time vectors
    predatortime=PredatorArray(4,:);
    preytime=PreyArray(4,:);
    % distance vectors
    predatordistance=PredatorArray(2,:);
    preydistance=PreyArray(2,:);
    % predator detection distance
    DetectiondistancePredator=Predator(9);
    % prey detection distance (initial distance when chase begins)
    DetectiondistancePrey=Prey(9);
    preydistance=preydistance+DetectiondistancePrey;
    % predator work
    predatorwork=PredatorArray(3,:);

    %%% find a crossing point
    % if Predator detection distance > Prey detection distance:
    % the predator can detect the prey before the prey can detect the predator 
    if DetectiondistancePredator>=DetectiondistancePrey
        % fit polynomial functions
        P1=polyfit(predatortime,predatordistance,5);
        P2=polyfit(preytime,preydistance,5);
        % check if coefficients are similar (i.e., parallel curves)
        i_check=abs(P1-P2);
        % parallel curves (no crossing point)
        if max(i_check(1:5))<1e-15
            distanceroot=NaN;
        else
            % if not parallel
            % determine the first root
            options=optimoptions(@fmincon,'Algorithm','active-set','MaxFunctionEvaluations',5000,'FunctionTolerance',1e-9);
            % constraints: min time >= 1e-10, max time <= predatortime(end) 
            A=[-1;1];
            b=[-1e-10;predatortime(end)];
            % starting point
            x0=1e-5;
            % Vector of parameters for fnVelocity function (see corresponding file for details)
            p=cell(1,2);
            p{1}=P1;
            p{2}=P2;
            % optimization
            root=fmincon(@(x)fnCapture(x,p),x0,A,b,[],[],[],[],[],options);
            % check whether or not curves cross (a root exists)
            if root<1e6
                 R1=polyval(P1,root);
                 R2=polyval(P2,root);
                 distanceroot=R2-R1;
            else
                distanceroot=NaN;
            end
        end

        % cross or close enough for capture
        threshold=Predator(3)/100;
        if ~isnan(distanceroot) && distanceroot < threshold
            % the predator can run over this period of time (feasible capture)
            if root<=predatortime(end)
                % the prey can run over this period of time
                if root<=preytime(end)
                    % find position in predator array
                    indextimepredator1=find(predatortime>=root);
                    indextimepredator=indextimepredator1(1);
                    % find position in prey array
                    indextimeprey1=find(preytime>=root);
                    indextimeprey=indextimeprey1(1);
                    %capturecostperattempt=predatorwork(indextimepredator)+fieldmetab(i)*root;
                    % predator speed > 0
                    if PredatorArray(1,indextimepredator)>0
                        % capture probability
                        captureprobability(i,1)=1.0/(1.0+(PreyArray(1,indextimeprey)/PredatorArray(1,indextimepredator)));
                        capturetime(i,1)=PredatorArray(4,indextimepredator);
                    else
                        % predator speed = 0
                        captureprobability(i,1)=0.0;
                        capturetime(i,1)=0.0;
                    end
                else
                    % the prey stops before the predator can reach it
                    % prey speed = 0
                    captureprobability(i,1)=1.0;
                    capturetime(i,1)=PredatorArray(4,indextimepredator);
                end
            else
                % the predator stops before
                captureprobability(i,1)=0.0;
                capturetime(i,1)=PredatorArray(4,end);
            end
        else
            % no crossing point
            captureprobability(i,1)=0.0;
            capturetime(i,1)=0.0;
        end
    else
        % predator detection distance < prey detection distance
        % no crossing point: the prey escapes before detection
        captureprobability(i,1)=0.0;
        capturetime(i,1)=0.0;
    end
end

%%%%%%%% Handling time
% Defines a vector for consumption time, a vector for digestion time
% and a vector for handling time that is equal to consumptiontime + digestiontime
handlingtime=zeros(NumberofPredators,1);
consumptiontime=zeros(NumberofPredators,1);
digestiontime=zeros(NumberofPredators,1);

for i=1:NumberofPredators
    Predator=PredatorTraits{i};    
    Prey=PreyTraits{i};
    %%%% handling time calculation
    %%% consumption time = number of bites * bite time
    % if predator bite size < prey size (multiple bites)
    if Predator(7)<Prey(2)
        consumptiontime(i,1)=(Prey(2)/Predator(7))*Predator(8);
    else
        % prey consumed in one bite
        consumptiontime(i,1)=Predator(8);
    end
    %%% digestion time
    digestiontime(i,1)=2.3e4*(Prey(2)/Predator(2))*(Predator(2)^0.25);
    %%% handling time
    handlingtime(i,1)=consumptiontime(i,1)+digestiontime(i,1);   
end

%%%%%% functional response
% for predator i feeding on prey i
%%% attack rate
attackrate=zeros(NumberofPredators,1);
for i=1:NumberofPredators
    attackrate(i)=Betaparam(i)*captureprobability(i);
end

%%% time
totaltime=zeros(NumberofPredators,1);
for i=1:NumberofPredators
    totaltime(i)=capturetime(i)*handlingtime(i);
end


