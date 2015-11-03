%This code is a simulation of a quadrature transmitter/receiver. Currently
%only works for an even number of bits.

%Created on 11/3/2015 by RWB

clear all
close all

t = linspace(0,.001,1000);%sampling at 1000x the frequency
%t = 0:1/1000000:.001

A = 1;
f = 1000;
y1 = A*sin(2*pi*f*t); %+sin
y2 = -A*sin(2*pi*f*t);%-sin
y3 = A*cos(2*pi*f*t); %+cos
y4 = -A*cos(2*pi*f*t);%-cos

stream = [1,0,1,0,0,1,0,1,1,1,1,0,0,0,0,1]; %bit stream

%transmitter initials
len = size(stream,2);
ang = zeros(1,len/2);
count = 1;

%Receiver initials
sz = size(ang,2);
bits = zeros(1,len);
count2 = 1;

%Mimics the transmitter/Receiver
graph = [];
graphi = [];
graphq = [];
for i = 1:2:len
    %Mimics the Transmitter
    
    %Assigns the I and Q from the incoming bit stream (2 bits)
    msb = num2str(stream(i));%msb as a string
    lsb = num2str(stream(i+1));%lsb as a string
    x = [msb,lsb];%concatenates the string
    
    switch x
        case '11'
            I = y3;
            Q = y1;
        case '01'
            I = y4;
            Q = y1;
        case '00'
            I = y4;
            Q = y2;
        case '10'
            I = y3;
            Q = y2;
    end
    
    %Mixes the two signals by adding them together
    out = I+Q;%This is what is transmitted
    
    graph = [graph , out]; %phase angles
    graphi = [graphi, I]; %I
    graphq = [graphq, Q]; %Q
    
    %<Channel> --No delay, no noise--
    
    %Mimics the Receiver
    
    %Finds the phase angle from a cosine reference
    [M,I] = max(out);
    [Mref,Iref] = max(y3);
    tau = t(I) - t(Iref);%time delay between peaks
    ang(count) = 2*pi*f*tau;%stores phase angle into an array

    %Based on the angle, assigns bit values (see constellation)
    if ang(count) >= 0 && ang(count) < pi/2          %First Quadrant
        bits(count2:count2+1) = [1,1];
    elseif ang(count) >= pi/2 && ang(count) < pi     %Second Quadrant
        bits(count2:count2+1) = [0,1];
    elseif ang(count) >=pi && ang(count) < 3*pi/2    %Third Quadrant
        bits(count2:count2+1) = [0,0];
    elseif ang(count) >= 3*pi/2 && ang(count) < 2*pi %Fourth Quadrant
        bits(count2:count2+1) = [1,0];
    end
    
    count2 = count2+2;
    count = count+1;
end

%display(bits)

ti = linspace(0,.001*len/2,1000*len/2);%time for i
tq = linspace(0,.001*len/2,1000*len/2);%time for q
figure
plot(ti,graphi);
hold on; grid on
plot(tq,graphq,'LineWidth',2);
legend('I','Q')
xlabel('Time [ms]')
ylabel('Amplitude')

tnow = linspace(0,.001*len/2,1000*len/2);
% figure
plot(tnow,graph);

