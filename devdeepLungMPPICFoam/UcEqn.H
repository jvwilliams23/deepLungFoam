fvVectorMatrix UcEqn
(
    fvm::ddt(Uc) + fvm::div(phic, Uc)
  + continuousPhaseTurbulence->divDevRhoReff(Uc)
 ==
    (0.0/rhoc)*cloudSU
  + fvOptions(Uc) //JW 07/12/2020
);

UcEqn.relax();

fvOptions.constrain(UcEqn);  //JW 07/12/2020

volScalarField rAUc(1.0/UcEqn.A());
surfaceScalarField rAUcf("Dp", fvc::interpolate(rAUc));

surfaceScalarField phicForces
(
   fvc::flux(rAUc*0.0*cloudVolSUSu/rhoc) + rAUcf*(g & mesh.Sf())
);

if (pimple.momentumPredictor())
{
    solve
    (
        UcEqn
     ==
        fvc::reconstruct
        (
            phicForces/rAUcf - fvc::snGrad(p)*mesh.magSf()
        )
    );
    fvOptions.correct(Uc); //JW 07/12/2020
}
