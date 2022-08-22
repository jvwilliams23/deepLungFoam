/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    DPMFoam

Description
    Transient solver for the coupled transport of a single kinematic particle
    cloud including the effect of the volume fraction of particles on the
    continuous phase.

\*---------------------------------------------------------------------------*/

bool debugChecks = false;
bool oneWayCoupling = false;
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "PhaseIncompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "fvOptions.H" //JW 07/12/2020
#include "scalarIOList.H"
#include "WKFunctions.C"   // Windkessel function file
#include "WKBCFvPatchScalarField.H"

bool newBreathHoldFound = false;

#include "basicKinematicMPPICCloud.H"
#define basicKinematicTypeCloud basicKinematicMPPICCloud

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "cloudName",
        "name",
        "specify alternative cloud name. default is 'kinematicCloud'"
    );

    argList::addBoolOption
    (
        "debug",
        "print debug checks"
    );

    argList::addBoolOption
    (
        "uncoupled",
        "make flow one-way coupled"
    );


    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createWindkessel.H"  //Windkessel header file (has to be placed in this order since p depends on scalar list storage)
    #include "createFields.H"
    #include "initContinuityErrs.H"

    Info<< "\nStarting time loop\n" << endl;

    if (args.optionFound("debug"))
    {
        debugChecks = true;
    }

    get_area_ratios(mesh, windkesselProperties, lungProperties);

    scalar sineWaveCoCutoff = 0.2;
    scalar sineWaveFreq = 0.0;

    offsetTime = 0.0;

    Info << "initial breath hold starts and ends at : " << nextBreathHoldStart
        << tab << nextBreathHoldEnd << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"


/// NEED TO INCLUDE ACCUMULATED BREAHT HOLD TIME HERE
        offsetTime = t - accumulatedBreathHoldTime;


        runTime++;
        t = runTime.value();
        sineWaveFreq = mag(Foam::sin(2.0*3.141*offsetTime/breathingPeriod));
        if (sineWaveFreq < sineWaveCoCutoff)
        {
            maxCo = earlyStageMaxCo;
        }
        else
        {
            maxCo = developedFlowMaxCo;
        }

        #include "setDeltaT.H"
        dt = runTime.deltaTValue();
    
        if (breathHoldDuration > 0)
        {
            Info<< "Time = " << runTime.timeName() << tab 
            << "offsetTime = " << offsetTime << tab 
            << "accumulatedBreathHoldTime = " << accumulatedBreathHoldTime 
            << nl << endl;
        }
        else
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;
        }

        continuousPhaseTransport.correct();
        muc = rhoc*continuousPhaseTransport.nu();

        Info<< "Evolving " << kinematicCloud.name() << endl;
        kinematicCloud.evolve();

        // Update continuous phase volume fraction field

        alphac = max(1.0 - kinematicCloud.theta(), alphacMin);
        alphac.correctBoundaryConditions();
        alphacf = fvc::interpolate(alphac);
        alphaPhic = alphacf*phic;


        
        fvVectorMatrix cloudSU(kinematicCloud.SU(Uc));
        volVectorField cloudVolSUSu
        (
            IOobject
            (
                "cloudVolSUSu",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedVector
            (
                "0",
                cloudSU.dimensions()/dimVolume,
                Zero
            ),
            zeroGradientFvPatchVectorField::typeName
        );

        cloudVolSUSu.primitiveFieldRef() = -cloudSU.source()/mesh.V();
        cloudVolSUSu.correctBoundaryConditions();
        cloudSU.source() = Zero;
        

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {

            #include "UcEqn.H"



            execute_pressure_update(mesh,phic,store,Uc);

            // --- PISO loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                continuousPhaseTurbulence->correct();
            }
        }

        execute_at_end(mesh,phic,store);


        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;




    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
