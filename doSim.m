%% doSim.m
% runs HH myelin model using global variables from main
% returns nothing
% parameters: none
% ~~~~~~~~~ BEGIN: ~~~~~~~~~ 
function doSim()
%% global variables

global pulseStart myelin_curve_freq timeDesired c_m g_i g_out g_m Vrest N n0 m0 h0 axon_radius axon_length Cm Cm_myelin Ri Rout Rm_myelin GK_max GNa_max GL myelination myelin_thickness;

nodes = (1:N)';

    %% do square modulus myelin curves (as in, every n nodes, make one
    %myelinated) or inverse (every n nodes, make one unmyelinated).
    
    square_myelin_curve = zeros(size(nodes));
    inverted_square_myelin_curve = zeros(size(nodes));
    
    square_myelin_curve = mod(nodes,myelin_curve_freq)==0; %every freq, make one myelinated section
    inverted_square_myelin_curve = square_myelin_curve==0; %every freq, make one unmyelinated section
    
    % now multiply by actual thickness
    square_myelin_curve = square_myelin_curve * myelin_thickness;
    inverted_square_myelin_curve = inverted_square_myelin_curve * myelin_thickness;

    % use this specific curve as our myelination curve
    myelination = inverted_square_myelin_curve;
    %myelination = ones(size(nodes)) .* myelin_thickness;
    %myelination = zeros(size(nodes));
    %myelination(length(nodes)/2) = 0;
    myelination(end) = 0;%last node should be unmyelinated
    myelination(1)=0;%make sure first node is not myelinated so that AP propagates
    myelination

    %% calculating conductance (to build G matrix)
    % NOTE: all distance units are in cm
    node_length = axon_length / N; %length of each node in cm
    surface_area_per_node = 2 * pi * axon_radius * node_length;
    surface_area_with_myelin = 2 .* pi .* node_length .* (axon_radius + myelination); %vector for all nodes
    cross_sectional_area = pi * axon_radius^2;


    % calculate capacitance and conductances of cable from specific ones above
    % these values are now in ohms & Farads (geometry-independent)
    r_i = (Ri / cross_sectional_area) * node_length;
    r_out = (Rout / cross_sectional_area) * node_length; % make sure Rout units work out (same as Ri)
    r_m = Rm_myelin ./ surface_area_with_myelin; %for myelinated parts (this is a vector for all nodes)
    c_m = ones(size(nodes));
    c_m(myelination>0) = Cm_myelin .* surface_area_with_myelin(myelination>0); %for both myelinated and unmyelinated parts (this is a vector for all nodes)
    c_m(myelination==0) = Cm * surface_area_per_node;
    %fix capacitance for myelinated parts (unmyelinated parts will stay the
    %same)
    %c_m(myelination>0) = c_m(myelination>0) / surface_area_per_node .* surface_area_with_myelin(myelination>0); % for myelinated sections, correct surface area to the one using myelin surface area
    %c_m = c_m .* (membrane_thickness) ./ (myelination + membrane_thickness);% divide capacitance by new thickness (with myelin)


    % fix the units for max conductances
    GK_max = GK_max * surface_area_per_node;
    GNa_max = GNa_max * surface_area_per_node;
    GL = GL * surface_area_per_node;

    % convert resistance to conductance
    g_i = 1 / r_i; 
    g_out = 1 / r_out; 
    g_m = 1 ./ r_m; %for myelinated parts (this is a vector for all nodes)

    %% get the initialStates
    initialStates = ones(4*N,1); %set up initial state for CableDyn as a column vector (made up of [Vm; n; m; h] )
    initialStates(1:N) = Vrest;
    initialStates(N+1:2*N) = n0;
    initialStates(2*N+1:3*N) = m0;
    initialStates(3*N+1:4*N) = h0;

    %% loop through input current range, and use one at a time
    %spike_frequencies = zeros(size(input_range));
    %for i = 1:length(input_range)
        %% running the simulation
        %Io = input_range(1); %changing input current
        [time,s]=ode23s('MyelinDyn', timeDesired, initialStates);

        Vm = s(:,1:N);
        n = s(:, N+1:2*N);
        m = s(:, 2*N+1:3*N);
        h = s(:, 3*N+1:4*N);

        %% find average conduction velocity
    % first, find at what time (t_final) the highest peak for for the last node is.
    % Then, delta_t = t_final - t_impulse_start
    % delta_x = length of axon - 1 node length
    %velocity = delta_x / delta_t

    %find all local maxima
    [peaks, peak_indices] = findpeaks(Vm(:,end), 'MINPEAKHEIGHT', -0.03); %returns peak and index of highest peak in the last node (highest spike)
            
    %if there are any maxima, calculate conduction velocity and spiking
    %frequency.
    if(~isempty(peak_indices))
        
        if(timeDesired(peak_indices(1)) < pulseStart)
            if(length(peak_indices)==1)% if there's only one peak and it's before stimulus, then there are no real peaks.
                disp('cant calculate velocity since there are no peaks2!');
                return;
            else % if there is a peak (spontaneous) before stimulus, use only the peaks from the 2nd one onwards
                t_final = timeDesired(peak_indices(2:end));
            end
        else
                t_final = timeDesired(peak_indices(:));       
        end
        
        %now that the times of the true peaks are found,
        
        %calculate conduction velocity
        delta_t = t_final(1) - pulseStart;
        delta_x = axon_length - node_length;
        conduction_velocity = delta_x / delta_t;
        
        %calculate spike_frequency
        time_between_first_and_last_peak = t_final(end) - t_final(1);
        spike_frequency = (length(t_final)-1) / time_between_first_and_last_peak;
        
        %display them in console (DEBUG)
        disp(strcat('the conduction velocity is: ', num2str(conduction_velocity(1)), ' cm/sec and the spike frequency is: ', num2str(spike_frequency)));
    else % no peaks found
        disp('cant calculate velocity since there are no peaks1!');
    end
    
    %save the states up in this weird matrix and save it to file named
    %test.mat (so that data can be plotted later using:
    %replayStatesAnimation('test');
        states = [Vm n m h time];
        save('test','states');

    

end