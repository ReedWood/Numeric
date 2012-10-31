%% Numerical methods in physics
% Exercises on distributions and Central Limit Theorem
% Wolfgang.Mader@fdm.uni-freiburg.de

%% Draw random numbers from different distributions
numberReal = 10000;
% Draw realizations from a uniformly distributed random Variable, U(0,1)
uni = rand( numberReal, 1 );

% Draw realizations from a normaly distributed random Variablem, N(0,1)
gauss = randn( numberReal, 1 );

% Draw realizations from a Cauchy(0,1) distributed random Variable
% Here, we have to generate the distribution using uniformly distributed
% random numbers. How this can be done was shown in the lecture, but
% wikipida is your friend, too. Check the article of the Cauchy
% distribution.
%
% In short, the pdf of the Cauchy distribution is
%   F(x) = 1/pi * arctan( ( x - x_0 ) / g ) + 0.5,
% and hence, inversion leads to
%   F(y) = tan( pi * ( y - 0,5 ) ),
% with y ~ N(0,1).
cauchy = tan( pi * ( rand( numberReal, 1 ) - .5 ) );

%% Plot all realizations
figure();
plot( uni, '*' );
title( 'Realization of a uniformly distributed random variable, U(0,1)' );
xlabel( 'Samples' );
ylabel( 'Realization' );

figure();
plot( gauss, '*' );
title( 'Realization of a normaly distributed random variable, N(0,1)' );
xlabel( 'Samples' );
ylabel( 'Realization' );

figure();
plot( cauchy, '*' );
title( 'Realization of a Cauchy distributed random variable, Cauchy(0,1)' );
xlabel( 'Samples' );
ylabel( 'Realization' );


%% Historgramms for all realizations
figure();
hist( uni );
title( 'Uniformly distributed random variables' );

figure();
hist( gauss );
title( 'Normaly distributed random variables' );

figure();
hist( cauchy );
title( 'Cauchy distributed random variables' );




%% Calculate normalized sums of random variables
numberVar = 1000;

% Uniform distribution
% Here, we are asked to draw random numbers from a uniform distribution for
% which the mean is 0 and the variance is 1. The mean of a random variable
% is defined as
%   mean(x) = int_a^b x * p(x) * dx
%           = p * int_a^b x * dx
%           = p * 1/2 * (b^2 - a^2)
% which gives as a possible solution b = -a.
%
% For the variance we have
%   var(x) = int_a^b x^2 * p(x) * dx
%          = p * int_a^b x^2 * dx
%          = p * 1/3 * ( b^3 - a^3 )
% with p = 1 / (2 * a) we get
%   var(x) = 1 / (2 * a) * 2/3 * a^3
%          = 1/3 * a^2
% which has to be 1, and hence the limits must be
%   a = sqrt( 3 ).
% Mind the fact, that we have to shift the variables form [0,1] to
% [-.5,.5], and therefore have to multiply with the factor
%   f = 2*a = 2*sqrt( 3 ).
mUni = ( rand( numberVar, numberReal ) - .5 ) * 2 * sqrt( 3 );
sUni = sum( mUni ) / sqrt( numberVar );

% Normaly distribution
mGauss = randn( numberVar, numberReal );
sGauss = sum( mGauss ) / sqrt( numberVar );

% Cauchy distribution
% Since we want the "mean" and "variance" of the realization to be 0 and 1,
% and this can not be calculatey for the Cauchy distribution, because mean
% and variance have no meaning in this context, we have to force the
% realization post-hoc, by subtracting the mean and dividing by the
% standard devation for every random variable.
mCauchy = tan( pi * ( rand( numberVar, numberReal ) - .5 ) );
for i = 1 : numberVar
    mCauchy( i, : ) = mCauchy( i, : ) - mean( mCauchy( i, : ) );
    mCauchy( i, : ) = mCauchy( i, : ) / std( mCauchy( i, : ) );
end
sCauchy = sum( mCauchy ) / sqrt( numberVar );

%% Historgramms for all summs
figure();
hist( sUni );
title( 'Sum of 1000 uniformly distributed random variables' );

figure();
hist( sGauss );
title( 'Sum of 1000 normaly distributed random variables' );

figure();
hist( sCauchy );
title( 'Sum of 1000 Cauchy distributed random variables' );



%% QQ-Plots
% Get sorted normaly distributed random numbers for x-axes
xAxes = sort( randn( numberReal, 1 ) );

% QQ-Plot: Uniform distribution
figure();
plot( xAxes, sort( uni ) );
title( 'QQ-Plot of uniformly and normaly distributed random numbers' );
xlabel( 'Normal distribution' );
ylabel( 'Uniform distribution' );

% QQ-Plot: Normal distribution
figure();
plot( xAxes, sort( gauss ) );
title( 'QQ-Plot of normaly and normaly distributed random numbers' );
xlabel( 'Normal distribution' );
ylabel( 'Normal distribution' );

% QQ-Plot: Uniform distribution
figure();
plot( xAxes, sort( cauchy ) );
title( 'QQ-Plot of Cauchy and normaly distributed random numbers' );
xlabel( 'Normal distribution' );
ylabel( 'Cauchy distribution' );
