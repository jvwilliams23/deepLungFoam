/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Application
    Any custom solver with Windkessel boundary conditions

Description
    Header file for the initialisation of Windkessel properties.
    Loads the WindkesselProperties dictionary
    Extracts the number of outlets in the system
    Updates the global variable N_OUTLETS
    Assigns a scalarList for storage of pressure values
    Finds the timestep and updates the global variable dt
    Runs the function intiialiseWindkessel.C; creates the WK struct

Dependencies
    WKFunctions.C
    windkesselProperties dictionary (when running)

\*---------------------------------------------------------------------------*/


/* Load dictionary of Windkessel properties */

IOdictionary windkesselProperties
(
    IOobject
    (
        "windkesselProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary lungProperties
(
    IOobject
    (
        "lungProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

numLobes = readScalar(lungProperties.lookup("numLobes"));
breathingPeriod = readScalar(lungProperties.lookup("breathingPeriod"));
inhalationDuration = breathingPeriod / 2.0;
tidalVolume = readScalar(lungProperties.lookup("tidalVolume"));
breathHoldDuration = lungProperties.lookupOrDefault("breathHoldDuration", 0.0);
C_globalCmH20 = readScalar(lungProperties.lookup("C_globalCmH20"));
R_globalCmH20 = readScalar(lungProperties.lookup("R_globalCmH20"));

scalar accumulatedBreathHoldTime = 0.0;


/*
 simulation is sensitive to timestep when flow is developing or transitioning to
 exhalation. Therefore, as a simple fix we provide two maximum Courant numbers.
 One will be used when the flow is developed (between 0.15 to 0.85*t_inhale).
 The other will be used at early (and late) stages of inhalation.
*/

scalar earlyStageMaxCo = readScalar(lungProperties.lookup("earlyStageMaxCo"));
scalar developedFlowMaxCo = readScalar(lungProperties.lookup("developedFlowMaxCo"));
//scalar developedFlowMaxCo = 0.5;
//scalar earlyStageMaxCo = 0.1; // use this value as maxCo 

/* Finding the number of outlets */

 const wordList outletNames(windkesselProperties.toc());
 int sizeOf = outletNames.size();
 N_OUTLETS = sizeOf;
 Info << "Number of outlets are :" << N_OUTLETS << endl;

/* Define density for use with Windkessel */ 

IOdictionary transportProperties_extra
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rho
(
    "rho",
    transportProperties_extra.lookup("rho")
);


RHO_0 = rho.value();

/* Updating global array patch_names */

forAll(outletNames, item)
{
 // initialise list of length = number of outlets
 patch_names.push_back("-1");
}

Info << "assigning outIndex" << endl;
forAll(outletNames, item)
{
 const word& outletName = outletNames[item];

 const dictionary& subDict = windkesselProperties.subDict(outletName);
 scalar outlet_index = readScalar(subDict.lookup("outIndex"));

 patch_names[outlet_index] = outletName;
}

for (int it_check = 0; it_check < N_OUTLETS; it_check++)
{
  Info << "Patch array " << it_check << " is " << patch_names[it_check] << endl;
}

/* Creating scalar list for storage of pressure values*/

Info<< "Reading scalar list store \n" << endl;
scalarIOList store
(
    IOobject
    (
        "store",
        runTime.timeName(),
        mesh
    ),
    N_OUTLETS
);


/*Initialising the storage array*/

for (int it = 0; it <N_OUTLETS; it++)
{
  store[it] = 0;
}

/* Initialise the WK properties*/
initialise(windkesselProperties);
    

