
double RHO_0;
double dt;
double t;
double offsetTime (0.0);
int t_step;
scalar tidalVolume;
scalar breathingPeriod;
scalar inhalationDuration;
scalar breathHoldDuration;
int numBreathsTracked (0.0);
bool breathHoldFound (false);
scalar nextBreathHoldStart (0.0);
scalar nextBreathHoldEnd (0.0);
scalar accumulatedBreathHoldTime (0.0);
int numLobes;
#include <vector>
std::vector<float> lobe_area;
std::vector<float> lobe_vol_fraction; //how much of total volume is represented by lobe j
double drivingPressure;
double drivingPressure_previous;
double drivingPressure_previous2;
double flowRateWithTime;
double R_globalCmH20;
double C_globalCmH20;
int N_OUTLETS;
std::vector<string> patch_names;

/* Windkessel Structure Definition */
typedef struct {
  double Q_current;     /* Current  time stepproximal flow rate */
  double Q_previous;    /* Previous time step proximal flow rate */
  double Q_previous2;   /* 2 Previous time step proximal flow rate */
  double P_current;     /* Current time stepproximal pressure */
  double P_previous;      /* Previous time step proximal pressure */
  double P_previous2;   /* 2 Previous time step proximal pressure */
  double volumeCurrent;
  double volumePrevious;
    double volumePrevious2;
  int lobeIndex;
  double outletArea; /* outlet area */
  double areaRatio; /* ratio of outlet area to sum of all area in lobe */
  int id;     /* Windkessel element id */
  double R;   /* Resistance */
  double C;   /* Compliance */
  double Z;   /* Impedance */
} WindKessel;


WindKessel *wk;

void initialise(const dictionary& windkesselProperties)
{


 /* Initialising Windkessel object */

  int i;
  wk = (WindKessel*)malloc(N_OUTLETS*sizeof(WindKessel));
  for (i=0;i<N_OUTLETS;i++)
  {
    /* retrieve from Scheme instantiation : First instantiation variables */
    wk[i].Q_current   = 0;
    wk[i].Q_previous  = 0;
    wk[i].Q_previous2   = 0;
    wk[i].P_current   = 0; /* Proximal pressure */
    wk[i].P_previous  = 0;
    wk[i].P_previous2   = 0;
    wk[i].volumeCurrent = 0.0;
    wk[i].volumePrevious = 0.0;
        wk[i].volumePrevious2 = 0.0;
  } 

  /* Retrieving values from windkessel dictionary*/

 const wordList outletNames(windkesselProperties.toc());

 forAll(outletNames, item)
   {
     const word& outletName = outletNames[item];
     Info << "Evaluating properties for " << outletName << endl;

     const dictionary& subDict = windkesselProperties.subDict(outletName);

     scalar C = readScalar(subDict.lookup("C"));
     scalar R = readScalar(subDict.lookup("R"));
     scalar Z = readScalar(subDict.lookup("Z"));
     scalar real_index = readScalar(subDict.lookup("outIndex"));
     label lobeIndex = readScalar(subDict.lookup("lobeIndex"));


     int out_index = real_index;
   
     Info << "C, R, Z and index are " << C << ", " << R << ", " << Z << ", " << out_index << "." <<endl; 

     wk[out_index].R      = R;
     wk[out_index].C      = C;
     wk[out_index].Z      = Z;
     wk[out_index].id       = out_index;
     wk[out_index].lobeIndex = lobeIndex;
     wk[out_index].outletArea = 0.0; 
     wk[out_index].areaRatio = 0.0;
   }


 /* End of file*/

}

double calculate_flow_rate(int i, fvMesh & mesh, surfaceScalarField & phi, volVectorField& U)
{
  scalar outflow = 0.0;

  label outletPatch = mesh.boundaryMesh().findPatchID(patch_names[i]);

  if (outletPatch >=0)
  {
    outflow = sum(phi.boundaryField()[outletPatch]);
  }

  reduce(outflow, sumOp<scalar>());

  return outflow;

}

double integrate_pressure_euler(int i, double R_outlet, double C_outlet)
{
    // calculate pressure at outlet
    scalar P_current = R_outlet * wk[i].Q_current 
        + wk[i].volumeCurrent / C_outlet
        + drivingPressure;
  if (debugChecks)
  {
    //Info << "R_outlet " << R_outlet << tab << "C_outlet " << C_outlet << endl;
    Info << "Term1 " << R_outlet * wk[i].Q_current << tab
         << "Term2 " << wk[i].volumeCurrent / C_outlet << tab
         << "Term3 " << drivingPressure << endl;
  }

    return P_current;
}

double integrate_pressure_AB2(int i, double R_outlet, double C_outlet)
{
    // calculate pressure at outlet w/ second order AB
    scalar C1 = 1.0;
    scalar C2 = 0.0;
    //scalar C3 = 0.0; // leave here for later third order test
    if (t_step >= 1)
    {
        C1 = 3.0 / 2.0;
        C2 = -1.0 / 2.0;
    }
    
    scalar P_current =  
        C1 * (R_outlet * wk[i].Q_current
        + wk[i].volumeCurrent / C_outlet
        + drivingPressure) 
        + C2 * (R_outlet * wk[i].Q_previous
        + wk[i].volumePrevious / C_outlet
        + drivingPressure_previous);
    return P_current;
}

void Wk_pressure_update(int i, double rho, fvMesh & mesh, surfaceScalarField & phi, scalarIOList & store, volVectorField& U)
{
  // scalar p,dpc,dpq;
  scalar cmH20_to_pa = 98.0665;

  // define some global parameters
  scalar R_global = R_globalCmH20 * cmH20_to_pa / 1.0e-6;
  scalar C_global = C_globalCmH20 * 1.0e-6 / cmH20_to_pa;

  // get flow rate at outlet, and amount of volume exited the outlet
  wk[i].Q_current = calculate_flow_rate(i,mesh,phi,U); //max(calculate_flow_rate(i,mesh,phi), 0.0);

  /*
    When inhaling, cutoff negative flowrate at outlets
    During exhalation, cutoff positive flowrate at outlets
  */
  if (flowRateWithTime >= 0.0)
  {
    if (wk[i].Q_current < 0.0)
    {
      wk[i].Q_current = 0.0;
    }
  }
  else
  {
    if (wk[i].Q_current > 0.0)
    {
      wk[i].Q_current = 0.0;
    }
  }
 
  if (breathHoldFound) wk[i].Q_current = 0.0;

  if (debugChecks)
  {
    Info << "Flow rate at outlet " << tab << wk[i].Q_current << tab << "vol" << tab << wk[i].volumeCurrent << endl;
  }
  scalar R_outlet = R_global / wk[i].areaRatio;
  scalar C_outlet = C_global * wk[i].areaRatio;

  scalar P_outlet_current = integrate_pressure_euler(i, R_outlet, C_outlet); 

  if (breathHoldFound) wk[i].Q_current = 0.0;

  if (flowRateWithTime >= 0.0)
  {
    wk[i].P_current = min(P_outlet_current, 0.0);  
  }
  else
  {
    wk[i].P_current = max(P_outlet_current, 0.0); 
  }

  /*Saving the pressure in a scalar array*/
  store[i] = wk[i].P_current;

  // debug check
  if (debugChecks)
  {
    Info << "Driving pressure: " << drivingPressure << tab 
         << "outlet pressure: " << wk[i].P_current << endl; 
  }

}

void get_area_ratios(fvMesh & mesh, const dictionary& windkesselProperties, const dictionary& lungProperties)
{
  int i;
    //const dictionary& subDictLobeVol = lungProperties.subDict("lobeVols").lookup("lobeVols");
    const List<scalar> lobeVols = lungProperties.lookup("lobeVols");

    // check number of entries is correct
    if (lobeVols.size() != numLobes)
    {
        FatalErrorInFunction
            << "Incorrect num entries for lobeVols in lungProperties" << endl
            << "Num entries is " << lobeVols.size() << " but numLobes is " << numLobes << endl
            << abort(FatalError);
    }

    //const List<scalarField> lobeVolNames(subDictLobeVol);

    // initialise scalar to check total volume fraction equals 1
    scalar totalVolFraction = 0.0;

  // initialise vector of total areas per lobe
  for (i=0;i<numLobes;i++)
  {
    lobe_area.push_back(0.0);
        lobe_vol_fraction.push_back(lobeVols[i]);

        totalVolFraction += lobeVols[i];
        //const scalar lobeVol = lobeVolNames[i];
        Info << "lobeVol is " << lobeVols[i] << endl;
  }

    if (totalVolFraction < 0.95 or totalVolFraction > 1.05)
    {
        FatalErrorInFunction
            << "total lobar volume fraction does not equal one: " << totalVolFraction << endl
            << abort(FatalError);
    }

  const wordList outletNames(windkesselProperties.toc());

    // TODO: add check to see that there is a lobe for each numLobe and each lobeIndex in windkesselProperties

  // get outlet areas for all outlets
  Info << "Get outlet areas" << endl;
  forAll(outletNames, item)
  {
    if (debugChecks) Info << "item " << item << endl;
    const word& outletName = outletNames[item];
    if (debugChecks) Info << "outlet name " << outletName << endl;

    const dictionary& subDict = windkesselProperties.subDict(outletName);
    scalar real_index = readScalar(subDict.lookup("outIndex"));
    int out_index = real_index;
    int lobeIndex = wk[out_index].lobeIndex;

    if (debugChecks) Info << "get patch label" << endl;
    label patchID = mesh.boundaryMesh().findPatchID(outletName); 
    if (debugChecks) Info << "find patch area" << endl;
    wk[out_index].outletArea = gSum(mesh.magSf().boundaryField()[patchID]);
    Info << "patch " << patch_names[out_index] << tab 
         << "area " << wk[out_index].outletArea << endl;
    lobe_area[lobeIndex] += wk[out_index].outletArea;
  }

  // get ratio of area at outlet i compared to total outlet area at lobe j
  Info << "Get area ratio" << endl;
  forAll(outletNames, item)
  {
    const word& outletName = outletNames[item];

    const dictionary& subDict = windkesselProperties.subDict(outletName);

    scalar real_index = readScalar(subDict.lookup("outIndex"));
    int out_index = real_index;


    int lobeIndex = wk[out_index].lobeIndex;
    
    wk[out_index].areaRatio = wk[out_index].outletArea * lobeVols[lobeIndex] / lobe_area[lobeIndex];
    Info << "outlet " << out_index <<" area ratio is " << wk[out_index].areaRatio << tab << "lobe " << lobeIndex << endl;
  }

}


void execute_pressure_update(fvMesh & mesh, surfaceScalarField & phi, scalarIOList & store, volVectorField& U)
{
  int i;

  for (i=0;i<N_OUTLETS;i++)
    {
      
      // empty patch with no faces, so skip to avoid floating point errors
      if (wk[i].outletArea == 0)
      {
        continue;
      }
      if (debugChecks)
      {
        Info << "patch" << patch_names[i]<< endl;
      }

      Wk_pressure_update(i, RHO_0, mesh, phi, store, U);

    }
}


void execute_at_end(fvMesh & mesh, surfaceScalarField & phi, scalarIOList & store)
{
  // scalar pa_to_cmH20 = 0.010197162129779282;
  scalar cmH20_to_pa = 98.0665;

  scalar R_global = R_globalCmH20 * cmH20_to_pa / 1.0e-6;
  scalar C_global = C_globalCmH20 * 1.0e-6 / cmH20_to_pa;

  scalar pi = 3.141591;
  int i;

  // check if new breath hold is scheduled, and turn on bool switch
  if (not breathHoldFound and t >= nextBreathHoldStart and  breathHoldDuration > 0.0)
  {
      breathHoldFound = true;
      Info << "BREATH HOLD FOUND at t = " << t << tab << nextBreathHoldStart << endl;
  }

  // if no breath-hold, perform normal calculations to get flowrate etc
  if (not breathHoldFound)
  {
      scalar volumeWithTimePrevious = -0.5 * (
              tidalVolume * Foam::cos(2.0*pi*(offsetTime-dt)/breathingPeriod)
              - tidalVolume
              );
      scalar volumeWithTime = -0.5 * (
              tidalVolume * Foam::cos(2.0*pi*offsetTime/breathingPeriod)
              - tidalVolume
              );
      /// maybe should use Next instead of previous (forward integration)
      flowRateWithTime = (volumeWithTime - volumeWithTimePrevious) / dt;

      // update drivingPressure
      drivingPressure_previous2 = drivingPressure_previous;
      drivingPressure_previous = drivingPressure;
      drivingPressure = -1.0*R_global * flowRateWithTime - volumeWithTime / C_global;

      if (debugChecks)
      {
          Info << "volumeWithTime " << volumeWithTime << tab
              << "flowRateWithTime " << flowRateWithTime << endl;
        Info << "next breath hold starts and ends at : " << nextBreathHoldStart 
            << tab << nextBreathHoldEnd << endl;
      }

  }
  else
  {
      flowRateWithTime = 0.0;

      drivingPressure_previous2 = 0.0;
      drivingPressure_previous = 0.0;
      drivingPressure = 0.0;

      if (breathHoldFound and t >= nextBreathHoldEnd)
      {
          breathHoldFound = false;
          accumulatedBreathHoldTime += breathHoldDuration;
          numBreathsTracked ++;
          nextBreathHoldStart = (breathingPeriod * numBreathsTracked) 
                                + inhalationDuration + accumulatedBreathHoldTime;
          nextBreathHoldEnd = nextBreathHoldStart + breathHoldDuration;
          Info << "BREATH HOLD FINISHED. Num breaths " << numBreathsTracked 
          << tab << "accumulatedBreathHoldTime " << accumulatedBreathHoldTime 
          << tab << "breathHoldDuration " << breathHoldDuration << endl;

      }
  }

  for (i=0;i<N_OUTLETS;i++) 
    {
      
      // empty patch with no faces, so skip to avoid floating point errors
      if (wk[i].outletArea == 0)
      {
        continue;
      }

      /* Save previous states */
      wk[i].P_previous2 = wk[i].P_previous;
      wk[i].Q_previous2 = wk[i].Q_previous;

      wk[i].P_previous = wk[i].P_current;
      wk[i].Q_previous = wk[i].Q_current;

      wk[i].volumePrevious2 = wk[i].volumePrevious;
      wk[i].volumePrevious = wk[i].volumeCurrent;

      // wk[i].Q_current = max(calculate_flow_rate(i,mesh,phi,U), 0.0);
      wk[i].volumeCurrent = wk[i].volumePrevious + wk[i].Q_current*dt;
    }
    // track number of time-steps (for AB scheme coeffs)
    t_step ++;
}
