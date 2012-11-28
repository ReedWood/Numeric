function[trip, dist] = anneal( location, trip, dist, T )
% Special function to get shortest trip length using simulated annealing.

N = length( location );
numberPermutations = ceil( .1 * N );
if numberPermutations == 0; numberPermutations = 1; end
perFrom = round( 1 + rand( numberPermutations, 1 ) * ( N - 1 ) );
perTo = round( 1 + rand( numberPermutations, 1 ) * ( N - 1 ) );

perTrip = trip;
for i = 1 : numberPermutations
    perTrip( perTrip == perTo  ) = -1;
    perTrip( perTrip == perFrom ) = perTo;
    perTrip( perTrip == -1 ) = perFrom;
end

perTrip = [perTrip, perTrip(1)];
perDist = 0;
for i = 1 : N
    perDist = perDist + sqrt( sum( ( location( perTrip( i ), : ) - ...
                                     location( perTrip( i + 1 ), : ) ).^2 ) );
end
perTrip = perTrip( 1 : end-1 );

if( perDist < dist )
    trip = perTrip;
    dist = perDist;
else
    thres = exp( -( perDist - dist ) / T );
    if rand( 1 ) < thres
        trip = perTrip;
        dist = perDist;
    end
end
    