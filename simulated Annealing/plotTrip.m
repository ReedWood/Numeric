function plotTrip( location, trip )

figure();
hold all;
xlim( [min( location( :, 1 ) ), max( location( :, 1 ) )] );
ylim( [min( location( :, 2 ) ), max( location( :, 2 ) )] );
for i = 1 : length( trip )
    plot( location( trip( i ), 1 ), location( trip( i ), 2 ), '*' );
    waitforbuttonpress;    
end