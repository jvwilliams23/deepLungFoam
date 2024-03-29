/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    flow_solver_1.1

Description
    Large time-step transient solver for incompressible, turbulent flow, using
    the PIMPLE (merged PISO-SIMPLE) algorithm.
    Windkessel boundary condition that is calculated at the end of every timesetp
    and is applied to the boundary condition with a scalarIOList.
    Windkessel parameters are provided in a dictionary known as windkesselProperties
    that is located in the constant folder of the case files. 
    

\*---------------------------------------------------------------------------*/

bool debugChecks = false;
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "scalarIOList.H"
#include "WKFunctions.C"   // Windkessel function file
#include "WKBCFvPatchScalarField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    timeSelector::addOptions();
    argList::addBoolOption
    (
        "debug",
        "print debug checks"
    );

    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createWindkessel.H"  //Windkessel header file (has to be placed in this order since p depends on scalar list storage)
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    turbulence->validate();
    get_area_ratios(mesh, windkesselProperties, lungProperties);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //          SOLVING FOR FLUID FLOW
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    

Info<< "\nStarting time loop\n" << endl;

    debugChecks = args.optionFound("debug");

    scalar sineWaveCoCutoff = 0.2;
    scalar sineWaveFreq = 0.0;

    while (runTime.run())
    {
        t = runTime.value();
        sineWaveFreq = mag(Foam::sin(2.0*3.141*t/breathingPeriod));

        //Info << "maxCo " << maxCo << endl;

        #include "readTimeControls.H"
        #include "CourantNo.H"
        if (sineWaveFreq < sineWaveCoCutoff)
        {
            maxCo = earlyStageMaxCo;
        }
        else
        {
            maxCo = developedFlowMaxCo;
        }
        #include "setDeltaT.H"

        Info << "maxCo " << maxCo << tab << "sineWaveFreqCriterion " << sineWaveFreq << endl;

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        dt = runTime.deltaTValue();

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            execute_pressure_update(mesh,phi,store,U);

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
            // execute_pressure_update(mesh,phi,store);
            // p.correctBoundaryConditions();
            // U.correctBoundaryConditions();
        }

    execute_at_end(mesh,phi,store);

    /* Updating the Windkessel struct data structure*/
    //execute_at_end(mesh,phi,store);
    //p.correctBoundaryConditions();
    //U.correctBoundaryConditions();
    //p.relax();

    runTime.write();

    if ((dt<1.e-6) and debugChecks)
    {
        FatalErrorInFunction
            << " Timestep too small. Something is wrong with the setup or solver" 
            << nl << exit(FatalError);

            
    }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
