%% replayStatesAnimation.m
% plots out and animates results from HH myelin simulation
% parameter: filename_string name of file containing data
function replayStatesAnimation(filename_string)

load(filename_string);

global myelination pulseStart myelin_thickness;
%extract data from saved file (data is in states matrix)
N = (length( (states(1,:)))-1) / 4;
Vm = states(:,1:N);
n = states(:, N+1:2*N);
m = states(:, 2*N+1:3*N);
h = states(:, 3*N+1:4*N);
timeDesired = states(:,end);


%Vm(i,:) gives voltage for all nodes at time=i
%Vm(:,i) gives voltage for node i over all time

    %% calculate frequency of spikes
    % grab all local maxima greater than 0 Volts
    %[peaks, peak_indices] = findpeaks(Vm(:,end), 'MINPEAKHEIGHT', -0.03);
    %number_of_peaks = length(peaks);

    % now calculate frequency by dividing the number of peaks with the time
    % between the first and last peak
    %if(number_of_peaks>0)
        %time_between_first_and_last_peak = timeDesired(peak_indices(number_of_peaks)) - timeDesired(peak_indices(1));
    %    spike_frequencies(i) = (number_of_peaks-1) / time_between_first_and_last_peak;
    %else
    %    disp('no peaks found!');
    %end
    
    
%% plot all the voltage and gate probabilities at all nodes over time
    


figure(1);
voltage_offset = 0.15;
gate_offset = 1.1;

%plot the last node
%    plot(timeDesired, Vm(:,1));
    subplot(211);
    plot(myelination');
    xlabel('Node', 'FontSize', 16, 'FontName', 'Helvetica');
    ylabel('Myelin Thickness in cm', 'FontSize', 16, 'FontName', 'Helvetica');
    title('Myelination Pattern of Axon', 'FontSize', 16, 'FontName', 'Helvetica');
    axis([0 N -0.1*myelin_thickness max(max(myelination)*1.1, myelin_thickness)] ); %sets the axes for myelin vs nodes

    subplot(212);
    plot(timeDesired, Vm(:,end));
    hold on;
    plot([pulseStart pulseStart], [min(Vm(:,end)) max(Vm(:,end))], 'r');
    text(pulseStart/2, mean([min(Vm(:,end)) max(Vm(:,end))]),'Stimulus', 'FontSize', 16, 'FontName', 'Helvetica');
    xlabel('Time in seconds', 'FontSize', 16, 'FontName', 'Helvetica');
    ylabel('Voltage in Volts', 'FontSize', 16, 'FontName', 'Helvetica');
    title('Voltage over time at last axon section', 'FontSize', 16, 'FontName', 'Helvetica');  
    hold off;
%for i=1:N
%for i=[1 20]

    %if(mod(i,15)~=0)
    %    continue;
    %end
    %subplot(211);
%    hold on;
%    plot(timeDesired, Vm(:,i) + voltage_offset*i);
%    text(max(timeDesired),Vm(end,i)+ voltage_offset*i,strcat('node ', num2str(i)) );
    
%   if(i==1 && number_of_peaks>0)
%        plot([timeDesired(peak_indices(1)) timeDesired(peak_indices(1))], [min(Vm(:,1)) Vm(peak_indices(1),end)+(voltage_offset*N) ], 'r');% plot a red line along the the first peak of the first node
%    end
    
%     subplot(212);
%     hold on;
%     plot(timeDesired, n(:,i) +gate_offset*i, 'blue');
%     plot(timeDesired, m(:,i) +gate_offset*i, 'green');
%     plot(timeDesired, h(:,i) +gate_offset*i, 'red');
%     text(max(timeDesired),n(end,i) +gate_offset*i,num2str(i));
%     text(max(timeDesired),m(end,i) +gate_offset*i,num2str(i));
%     text(max(timeDesired),h(end,i) +gate_offset*i,num2str(i)); 
%end


%% plot animation of all nodes over time (animated)
figure(2);

voltagePlotHandle = plot(1:length(Vm(1,:)), Vm(1,:));
axis([0 length(Vm(1,:)) min(min(Vm(:,:))) max(max(Vm(:,:)))] ); %sets the axes for Vm vs nodes
xlabel('Node', 'FontSize', 16, 'FontName', 'Helvetica');
ylabel('Voltage in Volts', 'FontSize', 16, 'FontName', 'Helvetica');

pause(3);  %wait 3 seconds before starting animation

i=1;
while(1)
%for i = 1:length(Vm(:,1))
  if(i==length(timeDesired))
      i=1;
      continue;
  end
  
  %if(mod(i,2)~=0)
  %    i = i + 1;
  %    continue;
  %end
  
  % redraw the voltage for this time 
  figure(2);
  set(voltagePlotHandle, 'XData', 1:length(Vm(i,:)));
  set(voltagePlotHandle, 'YData', Vm(i,:));
  title(strcat('t = ', num2str(timeDesired(i)) ) );  
  
  pause( 1 / length(Vm(:,1)) ); %slight pause to give animation illusion
  i = i + 1;
  

end


end