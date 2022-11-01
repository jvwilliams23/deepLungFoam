      scalar volumeWithTimePrevious = -0.5 * (
              tidalVolume * Foam::cos(2.0*pi*(offsetTime-dt)/breathingPeriod)
              - tidalVolume
              );
      volumeWithTime = -0.5 * (
              tidalVolume * Foam::cos(2.0*pi*offsetTime/breathingPeriod)
              - tidalVolume
              );
      /// maybe should use Next instead of previous (forward integration)
      flowRateWithTime = (volumeWithTime - volumeWithTimePrevious) / dt;

