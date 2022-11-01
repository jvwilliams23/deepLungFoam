    scalar inhalationPeriod = 0.5 * breathingPeriod;

    scalar volumeWithTimePrevious = tidalVolume * (offsetTime-dt)/inhalationPeriod;
    volumeWithTime = tidalVolume * offsetTime/inhalationPeriod;

    flowRateWithTime = (volumeWithTime - volumeWithTimePrevious) / dt;
    Info << "flowRate " << flowRateWithTime << endl;

