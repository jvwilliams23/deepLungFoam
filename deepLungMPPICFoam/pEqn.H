{
    volVectorField HbyA(constrainHbyA(rAUc*UcEqn.H(), Uc, p));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
           fvc::flux(HbyA)
         + rAUcf*fvc::ddtCorr(Uc, phic)
        )
    );

    if (p.needReference())
    {
        adjustPhi(phiHbyA, Uc, p);
    }

    phiHbyA += phicForces;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, Uc, phiHbyA, rAUcf);

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUcf, p)
         ==
            fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phic = phiHbyA - pEqn.flux();

            p.relax();

            Uc = HbyA
              + rAUc*fvc::reconstruct((phicForces - pEqn.flux())/rAUcf);
            Uc.correctBoundaryConditions();
        }
    }
}

#include "continuityErrs.H"

fvOptions.correct(Uc); //JW 07/12/2020
