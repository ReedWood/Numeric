%% Simulated annealing and the traveling sailesman problem

%% Define coordinates of cities

location = [0  2 10  7 3 12 20 13 15 1;
            0 10 10 20 3  4 12 12 3  10]';
        

%% Find shortest trip using simulated annealing
T = 100;
trip = 1 : 10;
dist = 10000;

while T > 0
    [trip, dist] = anneal( location, trip, dist, T );
    T = T - .01;
end
