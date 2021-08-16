%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
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
% The following code defines 11 functions used in the calculation of motion speeds and costs 
% fnMotion
% fnVelocity 
% fnVelocityStiff
% fnDragCoef
% fnActiveAscending
% fnPassiveAscending
% fnRootPassiveMotion
% fnPassiveDescending
% fnRootPassiveDescending
% fnActiveForwardMotion
% fnPassiveForwardMotion
%
%
%%%% fnMotion is used to called from outside and drives the computation 
% The function takes two arguments (BodyMass and p)
%
%%% BodyMass is used to calculated Maximal Force for motion (MaxForce)
%
%%% p is vector containing biological and physical parameters
%%% It is used by fnVelocity and fnFelocityStiff functions
% 1. Body volume (m3)
% 2. Body mass (kg)
% 3. Body radius (m)
% 4. Body section surface (m2)
% 5. Stroke period (s)
% 6. Acceleration due to gravity (m.s-2)
% 7. Medium density (kg.m-3)
% 8. Medium dynamic viscosity (N.s.m-2)
% 9. Time step for calculation (s)
% 10. Switch parameter (see below)
% 
%
%%% The function returns an object (Result)
% According to the value of the Switch parameter
% If i_Switch=1: the function drives the computation of the species-specific speed  
% Result is a vector containing 2 values:
% 1. Species-specific speed
% 2. Motion cost per time at species-specific speed
%
% If i_Switch=2: the function drives the computation of the capture sequence
% Result is an array with 4 rows:
% 1. Instantaneous horizontal speed (m.s-1)
% 2. Cumulated horizontal distance travelled since the beginning (m)
% 3. Cumulated work (J)
% 4. Cumulated time (s)
%
%
% If i_Switch=5: the function drives the computation of the motion cost during handling time 
% Result is a single value:
% 1. Motion cost per time during handling time
%
%
%%%% fnVelocity is used to compute speed, forces, or motion cost 
%
% The function takes two vectors as arguments (x and p)
%
%%% x is a vector containing vertical and horizontal forces 
% 1. Vertical muscular force (N)
% 2. Horizontal muscular force (N)
%
%%% p is vector containing biological and physical parameters
% 1. Body volume (m3)
% 2. Body mass (kg)
% 3. Body radius (m)
% 4. Body section surface (m2)
% 5. Stroke period (s)
% 6. Acceleration due to gravity (m.s-2)
% 7. Medium density (kg.m-3)
% 8. Medium dynamic viscosity (N.s.m-2)
% 9. Time step for calculation (s)
% 10. Switch parameter (see below)
%
%%% The function returns different outputs according to the value of the switch parameter
%
%%% If i_swith = 1
% fnVelocity returns the ratio between work and horizontal speed: this value has to be minimized for species-specific speed (search sequence)
% This option is used during optimization procedure for species-specific speed
%
%%% If i_swith = 2
% fnVelocity returns 1 / total horizontal distance travelled: this value has to be minimized for predator jump (capture sequence)
% This option is used during optimization procedure of capture sequence
%
%%% If i_swith = 3
% fnVelocity returns a vector with horizontal speed and motion cost 
% This option is used during optimization procedure for capture sequence
%
%%% If i_swith = 4
% fnVelocity returns an array of 4 rows used during capture sequence (used after optimization procedure)
% This array is used for the calculation of capture cost
% The 4 rows give for each time step:
% 1. Instantaneous horizontal speed (m.s-1)
% 2. Cumulated horizontal distance travelled since the beginning (m)
% 3. Cumulated work (J)
% 4. Cumulated time (s)
%
% The ODE solver used for speed compution is ode45 (generic solver)
%
%
%%%% fnVelocityStiff is used to compute speed, forces, or motion cost 
% It is the same function as fnVelocity
% However, the ODE solver for speed compution is ode15s (for stiff functions)
%
%
%%%% fnDragCoef is used to compute the drag coefficient 
%
% The function takes body radius, speed, and a vector of physical parameters (Medium density and Medium dynamic viscosity)
% 
% The function returns a single value: drag coefficient
%
%
%%%% fnActiveAscending is used to compute the derivative of speed with respect to time during active ascending phase of motion (vertical plan)
% This ODE is called by the solver within fnVelocity function
%
% The function takes time, speed, and a vector of parameters
% Parameters: Body mass, Body volume, Body radius, Body section surface, Vertical muscular force, Medium density, Medium dynamic viscosity, and Acceleration due to gravity
%
% The function returns the value of the derivative of speed with respect to time (acceleration) as a list
%
%
%%%% fnPassiveAscending is used to compute the derivative of speed with respect to time during the passive (inertial) ascending phase (vertical plan)
% This ODE is called by the solver within fnVelocity function
%
% The function takes time, speed, and a vector of parameters
% Parameters: Body mass, Body volume, Body radius, Body section surface, Medium density, Medium dynamic viscosity, and Acceleration due to gravity
%
% The function returns the value of the derivative of speed with respect to time (acceleration) as a list
%
%
%%%% fnRootPassiveMotion is used to trigger a stop event to the computation (for passive ascending and passive forward motion phases)
% This root function is called by the solver within fnVelocity function
% 
% The function takes time, speed, and a vector of parameters (same as fnPassiveAscending, fnPassiveDescending and fnPassiveForwardMotion functions)
%
% It checks if speed <= 0: the individual body stops, which triggers the end of the computation
% 
%
%%%% fnPassiveDescending is used to compute the derivative of speed with respect to time during the passive (inertial) descending phase (vertical plan)
% This ODE is called by the solver within fnVelocity function
% The function takes time, speed, and a vector of parameters (same as fnPassiveAscending and fnPassiveForwardMotion functions)
%
% The function returns the value of the derivative of speed with respect to time (acceleration) as a list
%
%
%%%% fnRootPassiveDescending is used to trigger a stop event to the computation (for passive descending motion phase)
% This root function is called by the solver within fnVelocity function
% 
% The function takes time, speed, and a vector of parameters (same as fnPassiveAscending, fnPassiveDescending and fnPassiveForwardMotion functions)
%
% It checks if the derivative of speed (acceleration) <=1e-9 (almost 0): the system is assumed to be at steady state, which triggers the end of the computation
%
%
%%%% fnActiveForwardMotion is used to compute the derivative of speed with respect to time during the active forward motion phase (horizontal plan)
%
% The function takes time, speed, and a vector of parameters
% Parameters: Body mass, Body volume, Body radius, Body section surface, Horizontal muscular force, Medium density, Medium dynamic viscosity, and Acceleration due to gravity
%
% The function returns the value of the derivative of speed with respect to time (acceleration) as a list
%
%
%%%% fnPassiveForwardMotion is used to compute the derivative of speed with respect to time during the passive (inertial) forward motion phase (horizontal plan)
%
% The function takes time, speed, and a vector of parameters (same as fnPassiveAscending and fnPassiveDescending functions)
%
% The function returns the value of the derivative of speed with respect to time (acceleration) as a list
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% fnMotion function
function [Result] = fnMotion(BodyMass,p)
    % Maximal muscular output
    Maxforce=55*BodyMass;
    % Switch parameter for fnVelocity function
    % if i_Switch=1: calculation for species-specific speed
    % if i_Switch=2: calculation for capture sequence
    % if i_Switch=5: calculation for motion cost during handling time
    i_Switch=p(10);
    
    % species-specific speed
    if i_Switch==1
        % i_Switch = 1: the output of fnVelocity returns a vector with horizontal and vertical optimized forces
        
        % options for optimization
        options=optimoptions(@fmincon,'Algorithm','active-set','MaxFunctionEvaluations',5000);
        % constraints: VerticalForce>=0, HorizontalForce>=0, VerticalForce+HorizontalForce<=MaxForce
        A=[-1 0;0 -1;1 1];
        b=[0;0;Maxforce];
        % starting point
        x0=[Maxforce/10; Maxforce/10];
        % optimization
        % Clear last warning message
        lastwarn(''); 
        % switch variable to handle stiff functions
        flag=0;
        try
            x=fmincon(@(x)fnVelocity(x,p),x0,A,b,[],[],[],[],[],options);
            [warnMsg, warnId] = lastwarn;
            % if there is any integration issue (probably due to stiff functions)
            if ~isempty(warnMsg)
                flag=1;
            end
        catch
            % handle exceptions for stiff functions during optimization procedure
            flag=1;
        end

        if flag==0
            % i_Switch = 3: fnVelocity will return a vector with speed and motion cost, for optimized force values
            i_Switch=3;  
            %p=[Species(1),BodyMass,Species(3),Species(4),Species(5),Gravity,MediumDensity,DynamicViscosity,TimeStep,i_Switch];
            p(10)=i_Switch;
            Res=fnVelocity(x,p);
            % Species-specific speed
            speed=Res(1);
            % Motion cost per time
            cost=Res(2);
        else
            % handle stiff functions
            x=fmincon(@(x)fnVelocityStiff(x,p),x0,A,b,[],[],[],[],[],options);
            % i_Switch = 3: fnVelocityStiff will return a vector with speed and motion cost, for optimized force values
            i_Switch=3;  
            p(10)=i_Switch;
            Res=fnVelocityStiff(x,p);
            % Species-specific speed
            speed=Res(1);
            % Motion cost per time
            cost=Res(2);
        end
        % returns species-specific speed and cost per time
        Result=[speed,cost];
    end
    
    % capture sequence
    if i_Switch==2
        % options for optimization
        options=optimoptions(@fmincon,'Algorithm','active-set','MaxFunctionEvaluations',5000);
        % constraints: VerticalForce>=0, HorizontalForce>=0, VerticalForce+HorizontalForce<=MaxForce
        A=[-1 0;0 -1;1 1];
        b=[0;0;Maxforce];
        % starting point for vertical and horizontal forces
        x0=[Maxforce/5; Maxforce/5];
        % optimization
        x=fmincon(@(x)fnVelocity(x,p),x0,A,b,[],[],[],[],[],options);
        % i_Switch = 4: fnVelocity will retun a 4 rows array with speed, distance, work , and time, for optimized force values
        i_Switch=4;
        p(10)=i_Switch;
        Result=fnVelocity(x,p);
    end
    
    % motion cost during handling time
    if i_Switch==5
        Result=fnHandlingMotion(Maxforce,p);
    end
end

%% fnVelocity function
function [Output] = fnVelocity(x,p)
    
    % unwrap p
    BodyVolume=p(1);
    BodyMass=p(2);
    BodyRadius=p(3);
    BodySectionSurface=p(4);
    StrokePeriod=p(5);
    Gravity=p(6);
    MediumDensity=p(7);
    DynamicViscosity=p(8);
    TimeStep=p(9);
    i_Switch=p(10);
    
    % unwrap x
    VerticalForce=x(1);
    if i_Switch == 2
        % capture sequence optimization: 
        % HorizontalForce = MaxForce - VerticalForce
        HorizontalForce=(BodyMass*55)-x(1);
    else
        HorizontalForce=x(2);
    end
    

    % Define a flag to avoid non-biological conditions (e.g., negative speed)
    % 0: values are realistic
    % 1: values out of bounds
    flag=0;

    %%%%% Vertical component
    %%% Active phase
    % initial speed = 0
    y0=0; 
    % time vector for computation
    time=0:TimeStep:StrokePeriod;
    % vector of parameters for computation
    paramactiveascending=[BodyMass,BodyVolume,BodyRadius,BodySectionSurface,VerticalForce,MediumDensity,DynamicViscosity,Gravity];
    % solve speed with respect to time during ascending active phase (stroke period)
    options=odeset('RelTol',1e-8);
    [t1,res1] = ode45(@(t,y) fnActiveAscending(t,y,paramactiveascending), time, y0, options);
    % vertical speed at the end of the active phase
    verticalspeed=res1(end);
    % Instaneous distance travelled
    instantdistancevertical=res1*TimeStep;
    % Cumulated vertical distance
    totaldistance=cumsum(instantdistancevertical);
    % Total vertical distance travelled during active phase
    activedistancetravelled=totaldistance(end);
    % Check if vertical force allows effective motion (or if it is too weak)
    if (activedistancetravelled<=0)
      flag=1;
    end
    
    if flag==0
        %%% Inertial ascending phase
        % Speed = vertical speed at the end of the active phase
        y0=verticalspeed;
        time=0:TimeStep:TimeStep*10000;
        % Vector of parameters for active phases (without forces)
        parampassivephase=[BodyMass,BodyVolume,BodyRadius,BodySectionSurface,MediumDensity,DynamicViscosity,Gravity];
        % solve speed with respect to time during ascending passive (inertial) phase
        % A stop event is triggered when speed = 0 (individual body stops)
        options=odeset('RelTol',1e-8,'Events',@(t,y) fnRootPassiveMotion(t,y));
        [t2,res2,te,ye,ie] = ode45(@(t,y) fnPassiveAscending(t,y,parampassivephase), time, y0, options);
        % Total vertical distance travelled during intertial ascending phase 
        instantdistancevertical=res2*TimeStep;
        totaldistance=cumsum(instantdistancevertical);
        passivedistanceupwards=totaldistance(end-1);
        % Duration of inertial ascending phase
        passivetimeupwards=t2(end-1); 
    else
        passivedistanceupwards=0;
        passivetimeupwards=0;
    end
    
    % Total distance travelled during ascending (active + passive) phases 
    totalascendingdistance=activedistancetravelled+passivedistanceupwards;

    %%% Inertial descending phase
    y0=0;
    time=0:TimeStep:TimeStep*1e5;
    % solve speed with respect to time during descending passive (inertial) phase
    % A stop event is triggered when the derivative of speed = 0 (steady state: terminal sinking velocity)
    parampassivephase=[BodyMass,BodyVolume,BodyRadius,BodySectionSurface,MediumDensity,DynamicViscosity,Gravity];
    options=odeset('RelTol',1e-8,'Events',@(t,y) fnRootPassiveDescending(t,y,parampassivephase));
    [t3,res3,te,ye,ie] = ode45(@(t,y) fnPassiveDescending(t,y,parampassivephase), time, y0, options);
    % Total vertical distance travelled during intertial descending phase 
    instantdistancevertical=res3*TimeStep;
    totaldistance=cumsum(instantdistancevertical);
    % Check if distance travelled > total ascending distance travelled
    if totaldistance(end) > totalascendingdistance
        % Search when distance travelled = total ascending distance travelled
        index=find(totaldistance >= totalascendingdistance);
        passivetimedownwards=t3(index(1));
    else
        % Extrapolate when distance travelled = total ascending distance travelled (knowing that the individual body has reached its terminal velocity)
        remainingdistance=totalascendingdistance-totaldistance(end);
        remainingtime=remainingdistance/res3(end); 
        passivetimedownwards=t3(end)+remainingtime; 
    end

    % Total duration of vertical motion
    totaltime=StrokePeriod+passivetimeupwards+passivetimedownwards;
    
    %%%%% Horizontal component
    %%% Active phase
    y0=0;
    time=0:TimeStep:StrokePeriod;
    % vector of parameters (same as active ascending phase, except Horizontal muscular force)
    paramactiveforward=[BodyMass,BodyRadius,BodySectionSurface,HorizontalForce,MediumDensity,DynamicViscosity];
    % solve speed with respect to time during active phase (stroke period)
    options=odeset('RelTol',1e-8);
    [t4,res4] = ode45(@(t,y) fnActiveForwardMotion(t,y,paramactiveforward), time, y0, options);
    % horizontal speed at the end of the active phase
    horizontalspeed=res4(end); 
    % Total horizontal distance travelled during active motion phase forwards
    instantdistancehorizontal=res4*TimeStep;
    totaldistance=cumsum(instantdistancehorizontal);
    activedistancetravelledforwards=totaldistance(end);
    % Check if horizontal force allows effective motion (or if it is too weak)
    if (activedistancetravelledforwards<0)
      flag=1;
    end
    
   if flag==0
       %%% Inertial phase
       y0=horizontalspeed;
       time=0:TimeStep:totaltime;
       if y0 >0
           % solve speed with respect to time during passive (inertial) motion phase forwards
           % A stop event is triggered when speed = 0 (individual body stops)
           parampassivephaseforward=[BodyMass,BodyRadius,BodySectionSurface,MediumDensity,DynamicViscosity];
           options=odeset('RelTol',1e-8,'Events',@(t,y) fnRootPassiveMotion(t,y));
           [t5,res5,te,ye,ie] = ode45(@(t,y) fnPassiveForwardMotion(t,y,parampassivephaseforward), time, y0, options);
           % Total horizontal distance travelled during active motion phase forwards
           instantdistancehorizontal=res5*TimeStep;
           totaldistance=cumsum(instantdistancehorizontal);
           passivedistancetravelledforwards=totaldistance(end);
       else
           passivedistancetravelledforwards=0;
       end
       % Total distance travelled forwards
       totaldistanceforwards=activedistancetravelledforwards+passivedistancetravelledforwards;
   end
   
   %%%% Outputs
   
    %%% Case 1: i_Switch = 1
    % Optimization procedure for species-specific speed
    % The function returns the ratio between work and speed
    if i_Switch==1
      if flag==0
          Workpt=(activedistancetravelled*VerticalForce+activedistancetravelledforwards*HorizontalForce)/(totaltime);
          %HorizontalSpeed=totaldistanceforwards/totaltime;
          Output=Workpt/totaldistanceforwards;
      else
          % Values out of bounds lead to huge penalty
          Output=1e9;
      end
    end
    
  
    %%% Case 2: i_Switch = 2
    % Optimization procedure for capture sequence
    % The function returns - total horizontal distance (i.e., the maximal total horizontal distance)
    if i_Switch==2
        if flag==0
            Output=-totaldistanceforwards;
        else
            % Values out of bounds lead to huge penalty
            Output=1e9;
        end
    end
    
    %%% Case 3: i_Switch = 3
    % After optimization, function returns species-specific speed and motion cost per time
    if (i_Switch==3)
        Workpt=(activedistancetravelled*VerticalForce+activedistancetravelledforwards*HorizontalForce)/(totaltime);
        HorizontalSpeed=totaldistanceforwards/totaltime;
        Output=[HorizontalSpeed, Workpt];
    end

    %%% Case 4: i_Switch = 4
    % After optimization, function returns an array for capture sequence
    if i_Switch==4
        if flag==0
            %%% Instantneous speed
            instantspeed=[res4(2:end).',res5.'];
            instantspeed(end)=0;
            %%% cumulated horizontal distance
            % active phase
            instantdistancehorizontal=res4(2:end)*TimeStep;
            totaldistanceactive=cumsum(instantdistancehorizontal);
            finaldistanceactive=totaldistanceactive(end);
            % passive phase
            instantdistancehorizontal=res5*TimeStep;
            totaldistancepassive=cumsum(instantdistancehorizontal);
            totaldistancepassive=totaldistancepassive+finaldistanceactive;
            % cumulated distance
            cumulatedhorizontaldistance=[totaldistanceactive.',totaldistancepassive.'];

            %%% cumulated work
            cumulatedworkactivephase=totaldistanceactive*HorizontalForce;
            finalactivework=cumulatedworkactivephase(end);
            % passive phase: instaneous work = 0, thus cumulated work = final work at the end of active phase
            cumulatedworkpassivephase=ones(1,length(totaldistancepassive))*finalactivework;
            totalcumulatedwork=[cumulatedworkactivephase.',cumulatedworkpassivephase];

            %%% cumulated time
            finaltimeactive=t4(end);
            cumulatedtimepassive=t5;
            cumulatedtimepassive=cumulatedtimepassive+finaltimeactive;
            cumulatedtime=[t4(2:end).',cumulatedtimepassive.'];

            %%% adjust size of speed vector
            maxsize=max(length(instantspeed),length(cumulatedhorizontaldistance));
            if length(instantspeed)<maxsize
                missingdata=maxsize-length(instantspeed);
                instantspeed=[instantspeed,zeros(1,missingdata)];
            end
      
            %%% data frame
            % row 1: instantaneous speed (at any time step)
            % row 2: cumulated horizontal distance covered
            % row 3: cumulated work
            % row 4: cumulated time
            Output=vertcat(instantspeed,cumulatedhorizontaldistance,totalcumulatedwork,cumulatedtime);
        else
            Output=NaN;
        end
    end
  
end

%% fnVelocityStiff function
function [Output] = fnVelocityStiff(x,p)
  
    % unwrap p
    BodyVolume=p(1);
    BodyMass=p(2);
    BodyRadius=p(3);
    BodySectionSurface=p(4);
    StrokePeriod=p(5);
    Gravity=p(6);
    MediumDensity=p(7);
    DynamicViscosity=p(8);
    TimeStep=p(9);
    i_Switch=p(10);
    
    % unwrap x
    VerticalForce=x(1);
    if i_Switch == 2
        % capture sequence optimization: 
        % HorizontalForce = MaxForce - VerticalForce
        HorizontalForce=(BodyMass*55)-x(1);
    else
        HorizontalForce=x(2);
    end

    % Define a flag to avoid non-biological conditions (e.g., negative speed)
    % 0: values are realistic
    % 1: values out of bounds
    flag=0;

    %%%%% Vertical component
    %%% Active phase
    % initial speed = 0
    y0=0; 
    % time vector for computation
    time=0:TimeStep:StrokePeriod;
    % vector of parameters for computation
    paramactiveascending=[BodyMass,BodyVolume,BodyRadius,BodySectionSurface,VerticalForce,MediumDensity,DynamicViscosity,Gravity];
    % solve speed with respect to time during ascending active phase (stroke period)
    options=odeset('RelTol',1e-8);
    [t1,res1] = ode15s(@(t,y) fnActiveAscending(t,y,paramactiveascending), time, y0, options);
    % vertical speed at the end of the active phase
    verticalspeed=res1(end);
    % Instaneous distance travelled
    instantdistancevertical=res1*TimeStep;
    % Cumulated vertical distance
    totaldistance=cumsum(instantdistancevertical);
    % Total vertical distance travelled during active phase
    activedistancetravelled=totaldistance(end);
    % Check if vertical force allows effective motion (or if it is too weak)
    if (activedistancetravelled<=0)
      flag=1;
    end
    
    if flag==0
        %%% Inertial ascending phase
        % Speed = vertical speed at the end of the active phase
        y0=verticalspeed;
        time=0:TimeStep:TimeStep*10000;
        % Vector of parameters for active phases (without forces)
        parampassivephase=[BodyMass,BodyVolume,BodyRadius,BodySectionSurface,MediumDensity,DynamicViscosity,Gravity];
        % solve speed with respect to time during ascending passive (inertial) phase
        % A stop event is triggered when speed = 0 (individual body stops)
        options=odeset('RelTol',1e-8,'Events',@(t,y) fnRootPassiveMotion(t,y));
        [t2,res2,te,ye,ie] = ode15s(@(t,y) fnPassiveAscending(t,y,parampassivephase), time, y0, options);
        % Total vertical distance travelled during intertial ascending phase 
        instantdistancevertical=res2*TimeStep;
        totaldistance=cumsum(instantdistancevertical);
        passivedistanceupwards=totaldistance(end);
        % Duration of inertial ascending phase
        passivetimeupwards=t2(end); 
    else
        passivedistanceupwards=0;
        passivetimeupwards=0;
    end
    
    % Total distance travelled during ascending (active + passive) phases 
    totalascendingdistance=activedistancetravelled+passivedistanceupwards;

    %%% Inertial descending phase
    y0=0;
    time=0:TimeStep:TimeStep*1e5;
    % solve speed with respect to time during descending passive (inertial) phase
    % A stop event is triggered when the derivative of speed = 0 (steady state: terminal sinking velocity)
    parampassivephase=[BodyMass,BodyVolume,BodyRadius,BodySectionSurface,MediumDensity,DynamicViscosity,Gravity];
    options=odeset('RelTol',1e-8,'Events',@(t,y) fnRootPassiveDescending(t,y,parampassivephase));
    [t3,res3,te,ye,ie] = ode15s(@(t,y) fnPassiveDescending(t,y,parampassivephase), time, y0, options);
    % Total vertical distance travelled during intertial descending phase 
    instantdistancevertical=res3*TimeStep;
    totaldistance=cumsum(instantdistancevertical);
    % Check if distance travelled > total ascending distance travelled
    if totaldistance(end) > totalascendingdistance
        % Search when distance travelled = total ascending distance travelled
        index=find(totaldistance >= totalascendingdistance);
        passivetimedownwards=t3(index(1));
    else
        % Extrapolate when distance travelled = total ascending distance travelled (knowing that the individual body has reached its terminal velocity)
        remainingdistance=totalascendingdistance-totaldistance(end);
        remainingtime=remainingdistance/res3(end); 
        passivetimedownwards=t3(end)+remainingtime; 
    end

    % Total duration of vertical motion
    totaltime=StrokePeriod+passivetimeupwards+passivetimedownwards;
    
    %%%%% Horizontal component
    %%% Active phase
    y0=0;
    time=0:TimeStep:StrokePeriod;
    % vector of parameters (same as active ascending phase, except Horizontal muscular force)
    paramactiveforward=[BodyMass,BodyRadius,BodySectionSurface,HorizontalForce,MediumDensity,DynamicViscosity];
    % solve speed with respect to time during active phase (stroke period)
    options=odeset('RelTol',1e-8);
    [t4,res4] = ode15s(@(t,y) fnActiveForwardMotion(t,y,paramactiveforward), time, y0, options);
    % horizontal speed at the end of the active phase
    horizontalspeed=res4(end); 
    % Total horizontal distance travelled during active motion phase forwards
    instantdistancehorizontal=res4*TimeStep;
    totaldistance=cumsum(instantdistancehorizontal);
    activedistancetravelledforwards=totaldistance(end);
    % Check if horizontal force allows effective motion (or if it is too weak)
    if (activedistancetravelledforwards<0)
      flag=1;
    end
    
   if flag==0
       %%% Inertial phase
       y0=horizontalspeed;
       size1=totaltime/TimeStep;
       if size1 < 1e6
           time=0:TimeStep:totaltime;
       else
           TimeStep2=totaltime/1e6;
           time=0:TimeStep2:totaltime;
       end
       if (y0 >0)
           % solve speed with respect to time during passive (inertial) motion phase forwards
           % A stop event is triggered when speed = 0 (individual body stops)
           parampassivephaseforward=[BodyMass,BodyRadius,BodySectionSurface,MediumDensity,DynamicViscosity];
           options=odeset('RelTol',1e-8,'Events',@(t,y) fnRootPassiveMotion(t,y));
           [t5,res5,te,ye,ie] = ode15s(@(t,y) fnPassiveForwardMotion(t,y,parampassivephaseforward), time, y0, options);
           % Total horizontal distance travelled during active motion phase forwards
           instantdistancehorizontal=res5*TimeStep;
           totaldistance=cumsum(instantdistancehorizontal);
           passivedistancetravelledforwards=totaldistance(end);
       else
           passivedistancetravelledforwards=0;
       end
       % Total distance travelled forwards
       totaldistanceforwards=activedistancetravelledforwards+passivedistancetravelledforwards;
   end
   
   %%%% Outputs
   
    %%% Case 1: i_Switch = 1
    % Optimization procedure for species-specific speed
    % The function returns the ratio between work and speed
    if i_Switch==1
      if flag==0
          Workpt=(activedistancetravelled*VerticalForce+activedistancetravelledforwards*HorizontalForce)/(totaltime);
          %HorizontalSpeed=totaldistanceforwards/totaltime;
          Output=Workpt/totaldistanceforwards;
      else
          % Values out of bounds lead to huge penalty
          Output=1e9;
      end
    end
    
  
    %%% Case 2: i_Switch = 2
    % Optimization procedure for capture sequence
    % The function returns - total horizontal distance (i.e., the maximal total horizontal distance)
    if i_Switch==2
        if flag==0
            Output=-totaldistanceforwards;
        else
            % Values out of bounds lead to huge penalty
            Output=1e9;
        end
    end
    
    %%% Case 3: i_Switch = 3
    % After optimization, function returns species-specific speed and motion cost per time
    if (i_Switch==3)
        Workpt=(activedistancetravelled*VerticalForce+activedistancetravelledforwards*HorizontalForce)/(totaltime);
        HorizontalSpeed=totaldistanceforwards/totaltime;
        Output=[HorizontalSpeed, Workpt];
    end

    %%% Case 4: i_Switch = 4
    % After optimization, function returns an array for capture sequence
    if i_Switch==4
        if flag==0
            %%% Instantneous speed
            instantspeed=[res4(2:end).',res5.'];
            instantspeed(end)=0;
            %%% cumulated horizontal distance
            % active phase
            instantdistancehorizontal=res4(2:end)*TimeStep;
            totaldistanceactive=cumsum(instantdistancehorizontal);
            finaldistanceactive=totaldistanceactive(end);
            % passive phase
            instantdistancehorizontal=res5*TimeStep;
            totaldistancepassive=cumsum(instantdistancehorizontal);
            totaldistancepassive=totaldistancepassive+finaldistanceactive;
            % cumulated distance
            cumulatedhorizontaldistance=[totaldistanceactive.',totaldistancepassive.'];

            %%% cumulated work
            cumulatedworkactivephase=totaldistanceactive*HorizontalForce;
            finalactivework=cumulatedworkactivephase(end);
            % passive phase: instaneous work = 0, thus cumulated work = final work at the end of active phase
            cumulatedworkpassivephase=ones(1,length(totaldistancepassive))*finalactivework;
            totalcumulatedwork=[cumulatedworkactivephase.',cumulatedworkpassivephase];

            %%% cumulated time
            finaltimeactive=t4(end);
            cumulatedtimepassive=t5;
            cumulatedtimepassive=cumulatedtimepassive+finaltimeactive;
            cumulatedtime=[t4(2:end).',cumulatedtimepassive.'];

            %%% adjust size of speed vector
            maxsize=max(length(instantspeed),length(cumulatedhorizontaldistance));
            if length(instantspeed)<maxsize
                missingdata=maxsize-length(instantspeed);
                instantspeed=[instantspeed,zeros(1,missingdata)];
            end
      
            %%% data frame
            % row 1: instantaneous speed (at any time step)
            % row 2: cumulated horizontal distance covered
            % row 3: cumulated work
            % row 4: cumulated time
            Output=vertcat(instantspeed,cumulatedhorizontaldistance,totalcumulatedwork,cumulatedtime);
        else
            Output=NaN;
        end
    end
  
end

%% fnDragCoef function
function[DragCoef] = fnDragCoef(Radius, Speed, param)
  
  % unwrap physical parameters
  MediumDensity=param(1);
  DynamicViscosity=param(2);
  
  % Avoid negative values
  if (Speed<0.0)
    Speed=-Speed;
  end
  
  % Reynolds' number
  Reynolds = MediumDensity*Speed*Radius*2/DynamicViscosity;
  % Drag coefficient
  if (Reynolds>0.0)
    DragCoef=(0.352+(0.124+24.0/Reynolds)^0.5)^2.0;
  else % in case of speed = 0
    DragCoef=0.0;
  end
end
  

%% fnActiveAscending function
function [y] = fnActiveAscending(t,x,parms)
  
  % unwrap parameters
  BodyMass=parms(1);
  BodyVolume=parms(2);
  BodyRadius=parms(3);
  BodySectionSurface=parms(4);
  VerticalForce=parms(5);
  MediumDensity=parms(6);
  DynamicViscosity=parms(7);
  Gravity=parms(8);
  
  % compute drag coefficient
  paramDragCoef=[MediumDensity,DynamicViscosity];
  DragCoef=fnDragCoef(BodyRadius,x,paramDragCoef);
  
  % compute derivative
  y=VerticalForce/BodyMass-Gravity+MediumDensity*BodyVolume*Gravity/BodyMass-0.5*MediumDensity*BodySectionSurface*DragCoef*(x^2.0)/BodyMass;
end


%% fnPassiveAscending function
function [y] = fnPassiveAscending(t,x,parms)
  
  % unwrap parameters
  BodyMass=parms(1);
  BodyVolume=parms(2);
  BodyRadius=parms(3);
  BodySectionSurface=parms(4);
  MediumDensity=parms(5);
  DynamicViscosity=parms(6);
  Gravity=parms(7);
  
  % compute drag coefficient
  paramDragCoef=[MediumDensity,DynamicViscosity];
  DragCoef=fnDragCoef(BodyRadius,x,paramDragCoef);
  
  % compute derivative
  y=-Gravity+MediumDensity*BodyVolume*Gravity/BodyMass-0.5*MediumDensity*BodySectionSurface*DragCoef*(x^2.0)/BodyMass; 
end

%% fnRootPassiveMotion function
function [value,isterminal,direction] = fnRootPassiveMotion(t,y)
  % check if y<=0: individual body stops
  value=y;
  isterminal=1;
  direction=-1;
end
 


%% fnPassiveDescending function
function [y] = fnPassiveDescending(t,x,parms)
  
  % unwrap parameters
  BodyMass=parms(1);
  BodyVolume=parms(2);
  BodyRadius=parms(3);
  BodySectionSurface=parms(4);
  MediumDensity=parms(5);
  DynamicViscosity=parms(6);
  Gravity=parms(7);
  
  % compute drag coefficient
  paramDragCoef=[MediumDensity,DynamicViscosity];
  DragCoef=fnDragCoef(BodyRadius,x,paramDragCoef);
  
  % compute derivative
  y=Gravity-MediumDensity*BodyVolume*Gravity/BodyMass-0.5*MediumDensity*BodySectionSurface*DragCoef*(x^2.0)/BodyMass;
end

%% fnRootPassiveDescending function
function [value,isterminal,direction] = fnRootPassiveDescending(t,y,parms)
    % check value of the derivative
    % if derivative < 1e-9: triggers end of computation (the system is assumed to be at steady state)
    y1=fnPassiveDescending(t,y,parms);
    value=y1-1e-9;
    isterminal=1;
    direction=-1;
end

%% fnActiveForwardMotion function
function [y] = fnActiveForwardMotion(t,x,parms)
  
  % unwrap parameters
  BodyMass=parms(1);
  BodyRadius=parms(2);
  BodySectionSurface=parms(3);
  HorizontalForce=parms(4);
  MediumDensity=parms(5);
  DynamicViscosity=parms(6);
  
  % compute drag coefficient
  paramDragCoef=[MediumDensity,DynamicViscosity];
  DragCoef=fnDragCoef(BodyRadius,x,paramDragCoef);
  
  % compute derivative
  y=HorizontalForce/BodyMass-0.5*MediumDensity*BodySectionSurface*DragCoef*(x^2.0)/BodyMass;
end

%% fnPassiveForwardMotion function
function [y] = fnPassiveForwardMotion(t,x,parms)
  
  % unwrap parameters
  BodyMass=parms(1);
  BodyRadius=parms(2);
  BodySectionSurface=parms(3);
  MediumDensity=parms(4);
  DynamicViscosity=parms(5);
  
  % compute drag coefficient
  paramDragCoef=[MediumDensity,DynamicViscosity];
  DragCoef=fnDragCoef(BodyRadius,x,paramDragCoef);
  
  % compute derivative
  y=-0.5*MediumDensity*BodySectionSurface*DragCoef*(x^2.0)/BodyMass;
end
