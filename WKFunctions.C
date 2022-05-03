
double RHO_0;
double dt;
double t;
int t_step;
scalar tidalVolume;
scalar breathingPeriod;
int numLobes;
// float lobe_area;
// float lobe_area[0];
#include <vector>
std::vector<float> lobe_area;
std::vector<float> lobe_vol_fraction; //how much of total volume is represented by lobe j
// double time;
double drivingPressure;
double drivingPressure_previous;
double drivingPressure_previous2;
double R_globalCmH20;
double C_globalCmH20;
int N_OUTLETS;
//char* patch_names[] = {"OUTLET_ACA","OUTLET_MCA"};
// DynamicList<string> patch_names(10); // 10 has been set as the maximum limit of outlets that are expected
std::vector<string> patch_names;
// DynamicList<double> lobe_area

/* Windkessel Structure Definition */
typedef struct {
	double Q_current;    	/* Current  time stepproximal flow rate */
	double Q_previous; 	  /* Previous time step proximal flow rate */
	double Q_previous2;  	/* 2 Previous time step proximal flow rate */
	double P_current;		  /* Current time stepproximal pressure */
	double P_previous;		  /* Previous time step proximal pressure */
	double P_previous2;	  /* 2 Previous time step proximal pressure */
	double Pout_current;		/* Current back-pressure */
	double Pout_previous;	/* Previous back-pressure */
	double Pout_previous2;  /* 2 Previous back-pressure */
	double Pc_current;		/* Current back-pressure */
	double Pc_previous;	/* Previous back-pressure */
	double Pc_previous2;  /* 2 Previous back-pressure */
	double volumeCurrent;
	double volumePrevious;
    double volumePrevious2;
	int lobeIndex;
	double outletArea; /* outlet area */
	double areaRatio; /* ratio of outlet area to sum of all area in lobe */
	// double drivingPressure;
	int id;			/* Windkessel element id */
	double R;		/* Resistance */
	double C;		/* Compliance */
	double Z;		/* Impedance */
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
		wk[i].Q_current 	= 0;
		wk[i].Q_previous 	= 0;
		wk[i].Q_previous2 	= 0;
		wk[i].P_current 	= 0; /* Proximal pressure */
		wk[i].P_previous 	= 0;
		wk[i].P_previous2 	= 0;
		wk[i].Pout_current 	= 0;  /* Venous pressure */
		wk[i].Pout_previous = 0;
		wk[i].Pout_previous2 = 0;
		wk[i].Pc_current 	= 0;	/* Intramural pressure */
		wk[i].Pc_previous 	= 0;
		wk[i].Pc_previous2 	= 0;
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

     wk[out_index].R 			= R;
     wk[out_index].C 			= C;
     wk[out_index].Z 			= Z;
     wk[out_index].id 			= out_index;
     wk[out_index].lobeIndex = lobeIndex;
     wk[out_index].outletArea = 0.0; //mesh.magSf().boundaryField()[outletName];
     wk[out_index].areaRatio = 0.0;
   }


 /* End of file*/

}


double first_order_derivative(double x, double xp)
{
	double derivative;

	derivative = (x - xp)/dt;

	return derivative;
}

double front_first_order_derivative()
{
	double derivative;

	derivative = 1.0/dt;

	return derivative;
}

double back_first_order_derivative(double xp)
{
	double derivative;

	derivative = -xp/dt;

	return derivative;
}

double second_order_derivative(double x, double xp,double xp2)
{
	double derivative;

	derivative = (1.5*x - 2.0*xp + 0.5*xp2)/dt;

	return derivative;
}

double front_second_order_derivative()
{
	double derivative;

	derivative = 1.5/dt;

	return derivative;
}

double back_second_order_derivative(double xp, double xp2)
{
	double derivative;

	derivative = (-2.0*xp + 0.5*xp2)/dt;

	return derivative;
}

double derivative(double x, double xp, double xp2)
{
	double d;

	// d = first_order_derivative(x,xp);

	d = second_order_derivative(x,xp,xp2);

	return d;
}

double front_derivative()
{
	double d;
	//d = front_first_order_derivative();

	d = front_second_order_derivative();

	return d;
}

double back_derivative(double xp, double xp2)
{
		double d;
		//d = back_first_order_derivative(xp);
		d = back_second_order_derivative(xp,xp2);

		return d;
}

/*double calculate_flow_rate_additional(int i, fvMesh & mesh, volVectorField & U)
{
	scalar flux = 0.0;

	word outlet_name = patch_names[i];

	//access boundary elements 
	const fvPatchList& patches = mesh.boundary(); 

	//loop over boundaries 
	forAll(patches, patchI) 
	{ 
	    const fvPatch& cPatch = patches[patchI]; 

	    //check boundary name 
	    if(cPatch.name() == outlet_name) 
	      { 

		//Access boundary face area vectors 
		const vectorField& inletAreaVectors = cPatch.Sf(); 

		//loop over all faces of selected 
		//boundary 
		forAll(cPatch, faceI) 
		{ 

		    //calculate boundary face mass flux 
                    scalar cFaceFlux = U.boundaryField()[patchI][faceI] & inletAreaVectors[faceI]; 

                     //sum face mass fluxes 
                     flux += cFaceFlux; 
		} 

	      } 

	  } 

	Info<< "FlowRate: " << flux << " m^3/s" << endl; 


	return flux;

}*/

double calculate_flow_rate(int i, fvMesh & mesh, surfaceScalarField & phi)
{
	scalar outflow = 0.0;

	label outletPatch = mesh.boundaryMesh().findPatchID(patch_names[i]);
    labelList patchCells = mesh.boundaryMesh()[outletPatch].faceCells();

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

void Wk_pressure_update(int i, double rho, fvMesh & mesh, surfaceScalarField & phi, scalarIOList & store)
{
	// scalar p,dpc,dpq;
	scalar cmH20_to_pa = 98.0665;

	// define some global parameters
	scalar R_global = R_globalCmH20 * cmH20_to_pa / 1.0e-6;
	scalar C_global = C_globalCmH20 * 1.0e-6 / cmH20_to_pa;

	// get flow rate at outlet, and amount of volume exited the outlet
	wk[i].Q_current = max(calculate_flow_rate(i,mesh,phi), 0.0);
    if (wk[i].Q_current < 0)
    {
        FatalErrorInFunction
            << "Reversed flow at outlet" << endl
            << abort(FatalError);
    }
	wk[i].volumeCurrent = wk[i].volumePrevious + wk[i].Q_current*dt;
	scalar R_outlet = R_global / wk[i].areaRatio;
	scalar C_outlet = C_global * wk[i].areaRatio;

    scalar P_outlet_current = integrate_pressure_euler(i, R_outlet, C_outlet);
    wk[i].P_current = min(P_outlet_current, 0.0);

	/*Saving the pressure in a scalar array*/
	store[i] = wk[i].P_current;

	// debug check
	if (debugChecks)
	{
        Info << "Flow rate at outlet " << tab << wk[i].Q_current << tab << "vol" << tab << wk[i].volumeCurrent << endl;
		Info << "R_outlet " << R_outlet << tab << "C_outlet " << C_outlet << endl;
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
 	forAll(outletNames, item)
  {
		const word& outletName = outletNames[item];

		const dictionary& subDict = windkesselProperties.subDict(outletName);
		scalar real_index = readScalar(subDict.lookup("outIndex"));
		int out_index = real_index;
		int lobeIndex = wk[out_index].lobeIndex;

 		label patchID = mesh.boundaryMesh().findPatchID(outletName); 
		wk[out_index].outletArea = gSum(mesh.magSf().boundaryField()[patchID]);
		Info << "patch " << patch_names[out_index] << tab 
				 << "area " << wk[out_index].outletArea << endl;
		lobe_area[lobeIndex] += wk[out_index].outletArea;
  }

  // get ratio of area at outlet i compared to total outlet area at lobe j
 	forAll(outletNames, item)
  {
		const word& outletName = outletNames[item];

		const dictionary& subDict = windkesselProperties.subDict(outletName);

		scalar real_index = readScalar(subDict.lookup("outIndex"));
		int out_index = real_index;


		int lobeIndex = wk[out_index].lobeIndex;
		wk[out_index].areaRatio = wk[out_index].outletArea * lobeVols[lobeIndex] / lobe_area[lobeIndex];
		Info << "outlet area ratio is " << wk[out_index].areaRatio << endl;
	}

}


void execute_at_end(fvMesh & mesh, surfaceScalarField & phi, scalarIOList & store)
{
	// Update variables to new "previous" and "previous2", then update pressure at outlet

	// scalar pa_to_cmH20 = 0.010197162129779282;
	scalar cmH20_to_pa = 98.0665;

	scalar R_global = R_globalCmH20 * cmH20_to_pa / 1.0e-6;
	scalar C_global = C_globalCmH20 * 1.0e-6 / cmH20_to_pa;

  scalar pi = 3.141591;

	scalar volumeWithTimePrevious = -0.5 * (
		tidalVolume * Foam::cos(2.0*pi*(t-dt)/breathingPeriod)
		- tidalVolume
	);
	scalar volumeWithTime = -0.5 * (
		tidalVolume * Foam::cos(2.0*pi*t/breathingPeriod)
		- tidalVolume
	);


	int i;
	/// maybe should use Next instead of previous (forward integration)
	scalar flowRateWithTime = (volumeWithTime - volumeWithTimePrevious) / dt;

    // update drivingPressure
    drivingPressure_previous2 = drivingPressure_previous;
    drivingPressure_previous = drivingPressure;
	drivingPressure = -1.0*R_global * flowRateWithTime - volumeWithTime / C_global;
	if (debugChecks)
	{
		Info << "volumeWithTime " << volumeWithTime << tab 
				 << "flowRateWithTime " << flowRateWithTime << endl;
	}

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

      /* Save previous states */
      wk[i].P_previous2 = wk[i].P_previous;
      wk[i].Q_previous2 = wk[i].Q_previous;
      wk[i].Pout_previous2 = wk[i].Pout_previous;
      wk[i].Pc_previous2 = wk[i].Pc_previous;

      wk[i].P_previous = wk[i].P_current;
      wk[i].Q_previous = wk[i].Q_current;
      wk[i].Pout_previous = wk[i].Pout_current;
      wk[i].Pc_previous = wk[i].Pc_current;

      wk[i].volumePrevious2 = wk[i].volumePrevious;
      wk[i].volumePrevious = wk[i].volumeCurrent;

      /*Update WindKessel values*/

      Wk_pressure_update(i, RHO_0, mesh, phi,store);

    }
    // track number of time-steps (for AB scheme coeffs)
    t_step ++;
}
