function methods = generateMethods(parameters,m)


methods(1).CallName = 'M1FO';
methods(1).Name = 'M1';
methods(1).Space = @M1;
methods(1).Initial = parameters.Conditions.IC;
methods(1).BC = parameters.Conditions.BC;
methods(1).RealLimiterFlag = 42;
methods(1).FluxLimiterFlag = 0;
methods(1).Model = 'FullMoment';
methods(1).HighestMoment = 1;

methods(end+1).CallName = 'M1NoLimiter';
methods(end).Name = 'M1';
methods(end).Space = @M1;
methods(end).Initial = parameters.Conditions.IC;
methods(end).BC = parameters.Conditions.BC;
methods(end).RealLimiterFlag = 0;
methods(end).FluxLimiterFlag = 0;
methods(end).Model = 'FullMoment';
methods(end).HighestMoment = 1;


methods(end+1).CallName = 'M1TVBMs0';
methods(end).Name = 'M1';
methods(end).Space = @M1;
methods(end).Initial = parameters.Conditions.IC;
methods(end).BC = parameters.Conditions.BC;
methods(end).RealLimiterFlag = 0;
methods(end).FluxLimiterFlag = 1;
methods(end).Model = 'FullMoment';
methods(end).HighestMoment = 1;

methods(end+1).CallName = 'M1TVBMs50';
methods(end).Name = 'M1';
methods(end).Space = @M1;
methods(end).Initial = parameters.Conditions.IC;
methods(end).BC = parameters.Conditions.BC;
methods(end).RealLimiterFlag = 0;
methods(end).FluxLimiterFlag = 1+50*1i;
methods(end).Model = 'FullMoment';
methods(end).HighestMoment = 1;

methods(end+1).CallName = 'M1Real';
methods(end).Name = 'M1';
methods(end).Space = @M1;
methods(end).Initial = parameters.Conditions.IC;
methods(end).BC = parameters.Conditions.BC;
methods(end).RealLimiterFlag = 1;
methods(end).FluxLimiterFlag = 0;
methods(end).Model = 'FullMoment';
methods(end).HighestMoment = 1;


methods(end+1).CallName = 'M1CharLimiter';
methods(end).Name = 'M1';
methods(end).Space = @M1;
methods(end).Initial = parameters.Conditions.IC;
methods(end).BC = parameters.Conditions.BC;
methods(end).RealLimiterFlag = 0;
methods(end).FluxLimiterFlag = 2+0*1i;
methods(end).Model = 'FullMoment';
methods(end).HighestMoment = 1;

methods(end+1).CallName = 'M1CharLimiter50';
methods(end).Name = 'M1';
methods(end).Space = @M1;
methods(end).Initial = parameters.Conditions.IC;
methods(end).BC = parameters.Conditions.BC;
methods(end).RealLimiterFlag = 0;
methods(end).FluxLimiterFlag = 2+50*1i;
methods(end).Model = 'FullMoment';
methods(end).HighestMoment = 1;

methods(end+1).CallName = 'M1CharRealLimiter';
methods(end).Name = 'M1';
methods(end).Space = @M1;
methods(end).Initial = parameters.Conditions.IC;
methods(end).BC = parameters.Conditions.BC;
methods(end).RealLimiterFlag = 1;
methods(end).FluxLimiterFlag = 2+0*1i;
methods(end).Model = 'FullMoment';
methods(end).HighestMoment = 1;

methods(end+1).CallName = 'M1CharRealLimiter50';
methods(end).Name = 'M1';
methods(end).Space = @M1;
methods(end).Initial = parameters.Conditions.IC;
methods(end).BC = parameters.Conditions.BC;
methods(end).RealLimiterFlag = 1;
methods(end).FluxLimiterFlag = 2+50*1i;
methods(end).Model = 'FullMoment';
methods(end).HighestMoment = 1;


methods(end+1).CallName = 'M1TVBMs0Real';
methods(end).Name = 'M1';
methods(end).Space = @M1;
methods(end).Initial = parameters.Conditions.IC;
methods(end).BC = parameters.Conditions.BC;
methods(end).RealLimiterFlag = 1;
methods(end).FluxLimiterFlag = 1;
methods(end).Model = 'FullMoment';
methods(end).HighestMoment = 1;

methods(end+1).CallName = 'M1TVBMs50Real';
methods(end).Name = 'M1';
methods(end).Space = @M1;
methods(end).Initial = parameters.Conditions.IC;
methods(end).BC = parameters.Conditions.BC;
methods(end).RealLimiterFlag = 1;
methods(end).FluxLimiterFlag = 1+50*1i;
methods(end).Model = 'FullMoment';
methods(end).HighestMoment = 1;

cd('../..');
if nargin<2
    m = 1:length(methods);
end
methods = methods(m);
end